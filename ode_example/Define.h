#pragma once

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NON_CONFORMING_WCSTOK

#define _CRTDBG_MAP_ALLOC

///////////////////////////////////////////////////////////////////////////////////////////////////

// 1Grainの計算ログを出力
//#define GPU_LOG

// CPUで計算する
//#define CPU

// アニメーション用ログの全部出し
//#define ALL_ANIMATION

///////////////////////////////////////////////////////////////////////////////////////////////////
#define	FILE_READ_BUFFER				256

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算用
#define PI								3.14159f
#define GAMA							(1.76e7f / 1.0e9f)

#define MF								10
#define METH							(MF / 10)
#define MITER							(MF - (10 * METH))
#define MAXORD							12					// MORD[METH - 1];
//	const int MORD[2] {12, 5};

#define MAXCOR							3
#define MXNCF							10

#define MXSTP0							500
#define MXHNL0							10

#define RTO								0.5e-4f
#define ATO								0.5e-4f

#define RAND_ARRAY						97

#define SIZE_EL							13
#define SIZE_ELCO_D1					13
#define SIZE_ELCO_D2					12
#define SIZE_TESCO_D1					3
#define SIZE_TESCO_D2					12

#define DUMMY_DOUBLE					-9876.5

///////////////////////////////////////////////////////////////////////////////////////////////////
// 集計用
#define NM_2_A							10.0
#define SNR_START_TRACK					100

#define FOLDER_SIMON					"simon"
#define FOLDER_RESULT					"result"

#define EXTENSION_CSV					".csv"
#define EXTENSION_PAR					".par"
#define EXTENSION_MED					".med"
#define EXTENSION_TPF					".tpf"
#define EXTENSION_HFD					".hfd"
#define EXTENSION_GA_DAT				"_ga.dat"
#define EXTENSION_GP_DAT				"_gp.dat"
#define EXTENSION_VEC_MZ				"_mz.csv"
#define EXTENSION_END					".end"
#define EXTENSION_BPOL					".bpol"
#define EXTENSION_MS					".ms"
#define EXTENSION_KIN					".kin"
#define EXTENSION_CON					".con"
#define EXTENSION_LOG					".log"
#define EXTENSION_OZ					".oz"
#define EXTENSION_SNR					".snr"
#define EXTENSION_MAG					".mag"

#define FILE_VEC						"vec_"
#define FILE_DYM						"dym_"
#define FILE_TMP						"tmp_"

#define FILE_PARAMETER					"parameter.par"
#define FILE_HEAD						"head.csv"
#define FILE_OFFSET						"offset.csv"
#define FILE_PATTERN					"pattern.csv"

#define FILE_RESULT_WAVEFORM			"Waveform.csv"
#define FILE_RESULT_WAVEFORM_1CYCLE		"WaveformCycle.csv"
#define FILE_RESULT_TRACK				"Track.csv"
#define FILE_RESULT_AVERAGE				"Average.csv"
#define FILE_RESULT_SNR					"SNR.csv"

#define PATTERN_VEC_MZ					"_mz.csv"



