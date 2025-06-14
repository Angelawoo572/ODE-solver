/*
cism.c - A cache simulator

This program simulates a cache memory using:
- number of set index bits ('s'), number of sets = 2^s
- number of lines per set ('E'), associativity
- number of block offset bits ('b'), block size = 2 ^ b

It processes a memory tracefile consisting of load/store operations and reports:
- number of cache hits, misses, and evictions
- total dirty bytes currenty in cache
- total dirty bytes that have been evicted to memory

It supports LRU replacement policy and a write-back, write-allocate policy.

It accepts the following command-line arguments (same as csim-ref):
-s <s>     : Number of set index bits (2^s sets)
-E <E>     : Number of lines per set (associativity)
-b <b>     : Number of block offset bits (2^b bytes per block)
-t <tracefile> : File name of the memory trace to process
-v         : (optional) Verbose mode (prints each access and its result)
-h         : Print help message and exit

Usage: ./csim [-v] -s <s> -E <E> -b <b> -t <tracefile>

Key data structures:
- A cache is a 2D array: an array of sets, each containing E cache lines.
- Each cache line contains:
- valid bit
- tag
- dirty bit (1 if the block has been modified)
- LRU timestamp (used for eviction policy)

Dirty bit behavior:
- Any **store** operation marks the corresponding line as dirty.
- If a dirty line is evicted, its entire block is written back to memory.
The number of dirty bytes evicted = block size.
At simulation end, remaining dirty blocks are counted toward - dirty bytes in
cache.

LRU is implemented with a global access counter ('time'), and each line stores
the timestamp of its most recent access. On replacement, the line with the
oldest timestamp is evicted.

Important constraints:
- All memory is dynamically allocated to support large cache sizes
- Trace parsing strictly follows the "Op Addr,Size" format
- The 'size' field in traces is ignored since memory accesses are block-aligned

When modifying:
- Be careful not to incorrectly clear or ignore the dirty bit on replacement
- Memory must be freed and file handles closed to avoid leaks
- Verbose and summary output must match csim-ref exactly

This simulator uses a counter-based LRU implementation for simplicity and
performance. It's easier to manage than pointer-based LRU (like linked lists).

The dirty bit tracking is required to model realistic write-back behavior.
Instead of writing data immediately to memory on a store, the block is marked
dirty and only written back on eviction. This mirrors actual hardware caching
behavior and influences performance.

Using dynamically allocated memory (via malloc/calloc) allows scalability and
avoids stack overflows. Struct-based cache lines help keep code modular and
clear.

Author: Angela Wu
Andrew ID: qianruw
 */
#include "cachelab.h"
#include <errno.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/**
 * Maximum length of a line read from the trace file.
 * This ensures fgets does not overflow the buffer when reading a line.
 */
#define LINELEN 128

/**
 * Represents a single cache line.
 * Each line stores:
 *   - valid: whether this line currently holds a valid block
 *   - tag: the tag bits from the memory address
 *   - last_used: global timestamp to implement LRU
 *   - dirty: whether this line has been modified (used in write-back policy)
 */
typedef struct {
    int valid;
    unsigned long tag;
    unsigned long last_used;
    int dirty;
} line;
/**
 * A block is a pointer to a cache line.
 * A set contains E blocks (E lines).
 */
typedef line *block;
/**
 * Cache: an array of 2^s sets,
 * where each set is an array of E blocks.
 */
block *cache;

/**
 * Cache parameters:
 *   - s : number of set index bits (number of sets = 2^s)
 *   - E : number of lines per set (associativity)
 *   - b : number of block offset bits (block size = 2^b bytes)
 */
int s, E, b;
/**
 * Global timestamp used for LRU.
 * Incremented after each memory access, and stored in last_used.
 */
unsigned long time = 0;
/**
 * Verbose mode flag. Set to true if the '-v' option is on.
 * If enabled, each memory operation will be printed with its result.
 */
bool v = false;
/**
 * Name of the trace file to process.
 * This is set by the '-t' command-line argument.
 */
char *tracefile = NULL;
/**
 * Structure for storing cache simulation statistics:
 *   - number of hits, misses, evictions
 *   - number of dirty bytes remaining in cache
 *   - number of dirty bytes evicted to memory
 */
csim_stats_t stats;

/**
 * @brief Allocate and initialize the cache data structure.
 *
 * This function dynamically allocates memory for the entire cache.
 * The cache is represented as an array of 2^s sets,
 * where each set contains E cache lines.
 *
 * Each line is initialized to:
 *  - valid = false
 *  - tag = 0
 *  - last_used = 0
 *  - dirty = false
 *
 * Memory layout:
 *  cache (block*) points to an array of sets (each is an array of lines).
 *
 * Warning:
 *  - This function must be called before any simulation begins.
 *  - Must call free_cache() after simulation to prevent memory leaks.
 */
void init_cache(void) {
    int S = 1 << s;
    cache = malloc((size_t)S * sizeof(block));
    if (!cache) {
        fprintf(stderr, "Allocation failed: cache\n");
        exit(1);
    }
    for (int i = 0; i < S; i++) {
        cache[i] = malloc((size_t)E * sizeof(line));
        if (!cache[i]) {
            fprintf(stderr, "Allocation failed: set %d\n", i);
            exit(1);
        }
        for (int j = 0; j < E; j++) {
            // Initialize valid bit; all lines are invalid at simulation start
            cache[i][j].valid = false;
            // Tag defaults to 0, irrelevant when valid = 0
            cache[i][j].tag = 0;
            // Zero last_used; LRU timestamps are based on global time counter
            cache[i][j].last_used = 0;
            // Reset dirty bit to ensure no false dirty evictions
            cache[i][j].dirty = false;
        }
    }
}

/**
 * @brief Free all allocated memory for the cache.
 *
 * Iterates over each set and frees its array of cache lines.
 * Finally, frees the top-level array of sets.
 *
 * Assumes the cache was previously allocated by init_cache().
 * Does nothing if cache is NULL.
 *
 * Must be called exactly once after simulation ends to avoid memory leaks.
 */
void free_cache(void) {
    int S = 1 << s;
    for (int i = 0; i < S; i++) {
        free(cache[i]);
    }
    free(cache);
}

/**
 * @brief Simulate a single memory access (load or store) to the cache.
 * This function performs hit/miss/eviction logic, updates LRU,
 * and maintains correct dirty byte tracking for write-back behavior.
 *
 * @param address   The memory address being accessed.
 * @param is_store  True if the operation is a store (write), false for load.
 *
 * Global state affected:
 *  - stats: hit/miss/eviction/dirty counters updated
 *  - cache: line metadata (valid, tag, dirty, last_used) updated
 *  - time: incremented to reflect LRU clock
 */
void accessCache(size_t address, bool is_store) {
    size_t S = 1UL << s; // 2^s sets, Number of sets
    size_t B = 1UL << b; // 2^b blocks, Block size in bytes
    // Extract set index and tag from the address
    size_t set_idx = (address >> b) & (S - 1);
    size_t tag = address >> (s + b);
    block lines = cache[set_idx];
    // Search for a cache hit by checkin valid lines with matching tag
    int hit_idx = -1;
    for (int i = 0; i < E; i++) {
        if (lines[i].valid && lines[i].tag == tag) {
            hit_idx = i;
            break;
        }
    }
    if (hit_idx >= 0) {
        // cache hit;
        stats.hits++;
        /* Increment time before updating last_used,
           so that newer accesses have strictly higher timestamps, enforces a
           correct LRU ordering. */
        time++;
        lines[hit_idx].last_used = time;
        /* If this is a store and the line was previously clean,
           mark it dirty and track one full block as dirty bytes.
           This avoids double-counting repeated writes. */
        if (is_store && !lines[hit_idx].dirty) {
            lines[hit_idx].dirty = 1;
            stats.dirty_bytes += B;
        }
        return;
    }
    // cache miss
    stats.misses++;
    // Find empty line or choose LRU victim to place new block
    int empty_idx = -1;
    int lru_idx = 0;
    unsigned long min_time = lines[0].last_used;
    /* Search for an unused line first.
       If none found, track the least-recently used line for eviction. */
    for (int i = 0; i < E; i++) {
        if (!lines[i].valid) {
            empty_idx = i;
            break;
        }
        if (lines[i].last_used < min_time) {
            min_time = lines[i].last_used;
            lru_idx = i;
        }
    }

    int use_idx = (empty_idx >= 0 ? empty_idx : lru_idx);

    /* If we are replacing a valid line, this is an eviction.
       Dirty blocks must be written back to memory. */
    if (lines[use_idx].valid) {
        stats.evictions++;
        /* If the evicted line was dirty, increment dirty_evictions
           and remove its bytes from dirty_bytes still in cache. */
        if (lines[use_idx].dirty) {
            stats.dirty_evictions += B;
            stats.dirty_bytes -= B;
        }
    }

    // fill the line with new block
    lines[use_idx].valid = 1;
    lines[use_idx].tag = tag;
    /* As with hits, update last_used timestamp after bumping time,
       to reflect most recent access. */
    time++;
    lines[use_idx].last_used = time;
    lines[use_idx].dirty = 0; // initially clean

    // if store, mark dirty
    /* If this was a store, mark new line dirty and track dirty bytes. */
    if (is_store) {
        lines[use_idx].dirty = 1;
        stats.dirty_bytes += B;
    }
}

/**
 * @brief Parse and process a memory access trace file.
 *
 * Reads each memory access from the trace and simulates the corresponding
 * cache behavior by calling accessCache().
 *
 * Lines starting with 'L' (load) or 'S' (store) are processed.
 * lines are ignored or flagged as errors.
 *
 * @param trace Path to the trace file
 * @return 0 on full success; 1 if any parse errors were encountered
 */
int process_trace_file(const char *trace) {
    FILE *tfp = fopen(trace, "rt");
    /* Check file open success immediately so we don't pass NULL to fgets,
       which would result in undefined behavior or a crash. */
    if (!tfp) {
        fprintf(stderr, "Error opening '%s': %s\n", trace, strerror(errno));
        return 1;
    }
    char linebuf[LINELEN]; // How big should LINELEN be?
    // buffer for each line of the trace
    int parse_error = 0;
    while (fgets(linebuf, LINELEN, tfp)) {
        /* Skip lines that don't represent memory accesses:
           - empty lines (newline or null terminator) */
        if (linebuf[0] == '\n' || linebuf[0] == '\0')
            continue;
        char operation;
        unsigned long addr;
        int size;

        /* Parse each line strictly using sscanf.
           Reject lines that don't match the required format.
           This prevents undefined behavior due to uninitialized values. */
        if (sscanf(linebuf, " %c %lx,%d", &operation, &addr, &size) != 3) {
            fprintf(stderr, "Line: parse error: %s", linebuf);
            parse_error = 1;
            continue;
        }

        /* The simulator only supports 'L' (load) and 'S' (store) ops.
           Reject unsupported op codes early to prevent silent bugs. */
        if (operation != 'L' && operation != 'S') {
            fprintf(stderr, "Line: invalid op '%c'\n", operation);
            parse_error = 1;
            continue;
        }

        /* Simulate the access.
           Pass 'is_store = true' if operation is 'S'. */
        accessCache(addr, operation == 'S');
    }
    fclose(tfp);
    return parse_error;
}

/**
 * @brief Print usage information for the csim simulator.
 *
 * This function is called when the user provides '-h' or invalid arguments.
 * It prints to stdout; if used in an error path, stderr may be more
 * appropriate.
 *
 * Note: It prints to stdout using unformatted output functions to avoid printf
 * failure.
 */
void printUsage(void) {
    fputs("Usage: ./csim [-hv] -s <s> -E <E> -b <b> -t <tracefile>\n", stdout);
    fputs("Options:\n", stdout);
    fputs("  -h         Print this help message.\n", stdout);
    fputs(
        "  -v         Verbose mode: report effects of each memory operation\n",
        stdout);
    fputs("  -s <s>     Number of set index bits.(there are 2**s sets)\n",
          stdout);
    fputs("  -E <E>     Number of lines per set. (associativity)\n", stdout);
    fputs("  -b <b>     Number of block bits.(there are 2**b blocks)\n",
          stdout);
    fputs("  -t <tracefile>  File name of the memory trace to process\n",
          stdout);
}

/**
 * @brief Entry point for the cache simulator.
 *
 * Parses command line arguments and processes the memory trace file.
 * Validates all required arguments and initializes the cache simulator state.
 * Any failure (invalid args, trace errors) results in immediate termination
 * with usage info.
 */
int main(int argc, char **argv) {
    int opt;
    /* If no arguments are given, immediately print usage.
       This prevents undefined behavior due to missing required args. */
    if (argc == 1) {
        fprintf(stderr, "Error: No arguments provided.\n");
        printUsage();
        exit(1);
    }

    /* Use getopt to parse command-line arguments.
       All required arguments (-s, -E, -b, -t) are captured by their respective
       cases. If any unrecognized flag is passed, fall through to default case.
     */
    while ((opt = getopt(argc, argv, "s:b:E:t:hv")) != -1) {
        switch (opt) {
        case 's':
            s = atoi(optarg);
            break;
        case 'E':
            E = atoi(optarg);
            break;
        case 'b':
            b = atoi(optarg);
            break;
        case 't':
            tracefile = optarg;
            break;
        case 'h':
            printUsage();
            exit(0);
        case 'v':
            v = true;
            break;
        default:
            /* Unknown option: usage string will guide user. */
            printUsage();
            exit(1); // false
        }
    }

    /* If there are unexpected extra arguments after all options,
      we reject them to avoid misinterpreting the input structure. */
    if (optind < argc) {
        fprintf(stderr, "Error: unexpected argument: %s\n", argv[optind]);
        printUsage();
        exit(1);
    }

    /* Validate required arguments:
       - s, E, b must be non-negative
       - s + b must not exceed 64 (bit width limit for tag computation)
       - tracefile must be provided
       These checks prevent silent memory corruption or segmentation faults. */
    if (s < 0 || E <= 0 || b < 0 || s + b > 64 || tracefile == NULL) {
        if (s < 0)
            fprintf(
                stderr,
                "Error: s (number of set index bits) must be non-negative.\n");
        if (E <= 0)
            fprintf(stderr,
                    "Error: E (number of lines per set) must be positive.\n");
        if (b < 0)
            fprintf(stderr,
                    "Error: b (number of block bits) must be non-negative.\n");
        if (s + b > 64)
            fprintf(stderr,
                    "Error: s + b must not exceed 64 (total address width).\n");
        if (tracefile == NULL)
            fprintf(stderr,
                    "Error: missing required trace file (-t <tracefile>).\n");
        printUsage();
        exit(1);
    }

    /* Initialize the cache structure.
       If init_cache ever fails (e.g., malloc fails), it must exit internally.
     */
    init_cache();

    /* Process the trace file. If it returns nonzero parse error,
       abort immediately with usage info to ensure we do not proceed on invalid
       state. */
    if (process_trace_file(tracefile)) {
        fprintf(stderr, "Error: trace file '%s' is malformed or unreadable.\n",
                tracefile);
        printUsage();
        exit(1);
    }

    /* Clean up memory before exiting. */
    free_cache();
    /* Report final simulation statistics. */
    printSummary(&stats);

    return 0;
}
