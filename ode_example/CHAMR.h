#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include <memory>
#include <vector>

#include <string.h>

#include "Define.h"


class CHAMR
{
public:
	CHAMR();
	virtual ~CHAMR();

	enum
	{
	};

	struct TRACK
	{
		TRACK()
		{
			memset(szName, 0, 64);
			fOffset = 0.0f;
			fBitlength = 0.0f;
		}

		char		szName[64];
		float		fOffset;
		float		fBitlength;
	};
	using VectorTrack = std::vector<TRACK>;

public:
	bool Run(const char* pszPath, const char* pszPattern, const char* pszParameter, const int nLogGrain);

	void SetParFile(const char* psz) {strcpy(m_fpa_in, psz);}
	void SetEndFile(const char* psz) {strcpy(m_fiend, psz);}

	const char* GetParFile() {return m_fpa_in;}
	const char* GetEndFile() {return m_fiend;}

protected:

	void myavg(float &amx, float &amy, float &amz);
	float ran1(int &IDUM);
	float gasdev(int &idum);

	float rumach();

#ifdef CPU

	void caltemp(float asx, float asy);
	void fun(float YDx[], float YDy[], float YDz[]);
	bool slsode(float T, float TOUT, int &ISTATE);
	void scfode();
	int sintdy (float T);
	bool sstode();
	void sewset();
	float svnorm(float Vx[], float Vy[], float Vz[], float Wx[], float Wy[], float Wz[]);


#endif

protected:
	std::unique_ptr<float[]>		m_atpx;
	std::unique_ptr<float[]>		m_atpy;
	std::unique_ptr<float[]>		m_atpz;
	std::unique_ptr<float[]>		m_atpt;

	// å˙Ç≥
	std::unique_ptr<float[]>		m_Thickness;
	std::unique_ptr<float[]>		m_tpf;

	// í∏ì_
	std::unique_ptr<float[]>		m_vx;
	std::unique_ptr<float[]>		m_vy;

	// _é•â◊
	std::unique_ptr<int[]>			m_nnv;
	std::unique_ptr<float[]>		m_ymx;
	std::unique_ptr<float[]>		m_ymy;
	std::unique_ptr<float[]>		m_ymz;

	// Exchange Field
	std::unique_ptr<float[]>		m_heb;
	std::unique_ptr<float[]>		m_hko;

	// .med
	std::unique_ptr<float[]>		m_tempc;
	std::unique_ptr<float[]>		m_sigmahk;
	std::unique_ptr<float[]>		m_ahko;
	std::unique_ptr<float[]>		m_sigmatc;
	std::unique_ptr<int[]>			m_amsrt;
	std::unique_ptr<float[]>		m_aexo;
	std::unique_ptr<float[]>		m_ept;
	std::unique_ptr<float[]>		m_eptAex;
	std::unique_ptr<float[]>		m_alf;
	std::unique_ptr<float[]>		m_alfg;
	std::unique_ptr<float[]>		m_alft;
	std::unique_ptr<float[]>		m_amso;

	// 
	std::unique_ptr<float[]>		m_tcbulk;
	std::unique_ptr<float[]>		m_sigmatc_gs;
	std::unique_ptr<float[]>		m_sigmatc_ots;

	// TC distribution
	std::unique_ptr<float[]>		m_tc;
	std::unique_ptr<float[]>		m_hks;

	// _ga
	std::unique_ptr<float[]>		m_gpx;
	std::unique_ptr<float[]>		m_gpy;
	std::unique_ptr<float[]>		m_gps;
	std::unique_ptr<float[]>		m_gvmi;
	std::unique_ptr<float[]>		m_gdia;
	std::unique_ptr<float[]>		m_gv;

	// head field
	std::unique_ptr<float[]>		m_ttw;
	std::unique_ptr<float[]>		m_hhw;

	// Thermal Field
	std::unique_ptr<float[]>		m_htho;

	// Start Simulation
	std::unique_ptr<float[]>		m_Yx;
	std::unique_ptr<float[]>		m_Yy;
	std::unique_ptr<float[]>		m_Yz;
	std::unique_ptr<float[]>		m_YHx;
	std::unique_ptr<float[]>		m_YHy;
	std::unique_ptr<float[]>		m_YHz;
	std::unique_ptr<float[]>		m_EWTx;
	std::unique_ptr<float[]>		m_EWTy;
	std::unique_ptr<float[]>		m_EWTz;
	std::unique_ptr<float[]>		m_SAVFx;
	std::unique_ptr<float[]>		m_SAVFy;
	std::unique_ptr<float[]>		m_SAVFz;
	std::unique_ptr<float[]>		m_ACORx;
	std::unique_ptr<float[]>		m_ACORy;
	std::unique_ptr<float[]>		m_ACORz;
	std::unique_ptr<float[]>		m_tmpt;
	std::unique_ptr<float[]>		m_bms;
	std::unique_ptr<float[]>		m_hthi;
	std::unique_ptr<float[]>		m_htx;
	std::unique_ptr<float[]>		m_hty;
	std::unique_ptr<float[]>		m_htz;
	std::unique_ptr<float[]>		m_hax;
	std::unique_ptr<float[]>		m_hay;
	std::unique_ptr<float[]>		m_haz;
	std::unique_ptr<float[]>		m_huk;
	std::unique_ptr<float[]>		m_hea;
	std::unique_ptr<float[]>		m_amk;
	std::unique_ptr<float[]>		m_amgx;
	std::unique_ptr<float[]>		m_amgy;
	std::unique_ptr<float[]>		m_amgz;
	std::unique_ptr<float[]>		m_amdz;
	std::unique_ptr<float[]>		m_ztmp;
	std::unique_ptr<float[]>		m_zhed;

	// .hfd
	std::unique_ptr<float[]>		m_hfix;
	std::unique_ptr<float[]>		m_hfiy;
	std::unique_ptr<float[]>		m_hfiz;
	std::unique_ptr<float[]>		m_zhfd;
	std::unique_ptr<float[]>		m_hftx;
	std::unique_ptr<float[]>		m_hfty;
	std::unique_ptr<float[]>		m_hftz;

	// Head Field Interpolation
	std::unique_ptr<float[]>		m_hfx;
	std::unique_ptr<float[]>		m_hfy;
	std::unique_ptr<float[]>		m_hfz;

	// Grain Coordinate
	std::unique_ptr<float[]>		m_gsx;
	std::unique_ptr<float[]>		m_gsy;

	// Random Seed
	std::unique_ptr<int[]>			m_rs;

	// 
	std::unique_ptr<float[]>		m_ELCO;
	std::unique_ptr<float[]>		m_TESCO;

	// 

//	int								m_nqt;				// m_nqs3

	int								m_ZCount;			// ëwêî				nz
	int								m_Grains;			// ÉZÉãêî
	int								m_Vertexs;			// í∏ì_ÇÃç≈ëÂêî
	int								m_nqs3;				// m_nz * 3
	int								m_nht;				// 

	int								m_ntpx;
	int								m_ntpy;
	int								m_ntpz;
	int								m_nhfx;
	int								m_nhfy;
	int								m_nhfz;
	float							m_MP_edge;
	float							m_tpf_edge;
	float							m_gymax;
	int								m_kxofs;
	int								m_nzhfd;

	int								m_ntmpw;
	int								m_nwr;

	float							m_vhksigma;			// 
	float							m_alftscaling;		// 
	float							m_dtc;				// 
	float							m_tambient;			// 
	float							m_tpk;				// 
	float							m_tgk;				// 
	float							m_spd;				// 
	float							m_pmove;			// 
	float							m_pstart;			// 
	float							m_offset;			// 
	float							m_gsz;				// 
	float							m_asz;				// 
	float							m_dtsimi;			// 
	float							m_antf;				// 
	float							m_hapi;				// 
	float							m_trise;			// 
	float							m_nps;				// 
	float							m_sbit;				// 
	float							m_sdelay;			// 
	int								m_npwrite;			// 
	int								m_nmwrite;			// 
	int								m_ntwrite;			// 
	int								m_animation;		// 
	int								m_ichki;			// 
	float							m_sigma_gs;			// 
	float							m_gstc_d;			// 
	float							m_gstc_v;			// 

	char							m_fpfx[256];
	char							m_fpa_in[256];

	char							m_fimed[256];
	char							m_fihfd[256];
	char							m_fitpf[256];
	char							m_fitpf_re[256];
	char							m_fiend[256];

	// slsode
	int								m_JSTART;
	int								m_KFLAG;
	int								m_L;
	int								m_METH;				// 1
	int								m_MITER;			// 0
	int								m_NQ;
	int								m_NST;
	int								m_NQU;

	float							m_EL0;
	float							m_H;
	float							m_HU;
	float							m_RC;
	float							m_TN;
	float							m_UROUND;

	int								s_NHNIL;
	int								s_NSLAST;
//	int								s_NYH;				// m_N nqs3

	// sstode
	int								s_IALTH;
//	int								s_NQNYH;			// m_NQ * nps3

	float							s_CONIT;
	float							s_CRATE;
	float							s_RMAX;

#ifdef CPU
	float							s_EL[13];
	float							s_ELCO[12][13];
	float							s_TESCO[12][3];
#endif

	size_t							m_nLLG;

};

