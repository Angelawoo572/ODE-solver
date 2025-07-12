# ODE-Solver with SUNDIALS

This project utilizes the SUNDIALS library ([https://github.com/LLNL/sundials](https://github.com/LLNL/sundials)) to numerically solve ordinary differential equations (ODEs).

## Installation

Please first download and install the SUNDIALS library from:

```
https://github.com/LLNL/sundials
```

Make sure SUNDIALS is correctly set up in your environment.

## Compilation

### ode\_try

To compile `ode_try`, run:

```
make -f Makefile
```

### ode\_formal

To compile and run `ode_formal`, execute:

```
make
make run
```

## Examples

### ode\_try

`ode_try` demonstrates solving a simple ODE of the form:

```
dy/dt = -y
```

This serves as a basic test to ensure your SUNDIALS setup is functioning properly.

### ode\_formal

`ode_formal` solves a more complex system of ODEs describing magnetization dynamics, using the Landau–Lifshitz–Gilbert (LLG) equation with constant parameters:

* `chk`, `che`, `alpha`, `chg`, `cha`, and `chb`.

The general form of the equation solved is:

LLG Equation

where parameters like `α` (damping constant) significantly affect system dynamics.

### Results

The results below illustrate the impact of varying the parameter `cha`:

* **Case 1:** `cha = 0` shows oscillatory behavior without damping.

* **Case 2:** `cha = 1.5` introduces significant damping and results in a distinct elliptical trajectory.

The visualizations highlight how parameter tuning alters the system's dynamic behaviors and trajectories.

## Usage

To run simulations, adjust parameters as needed within the source files, recompile using the commands provided above, and execute the resulting binaries to observe the output.

---

For further details, refer to the comments within each source file or contact the repository owner.
