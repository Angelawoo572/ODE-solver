#include "CUDAFunction.cuh"

// 乱数用
#define M1				259200
#define IA1				7141
#define IC1				54773
#define M2				134456
#define IA2				8121
#define IC2				28411
#define M3				243000
#define IA3				4561
#define IC3				51349
#define RM1				3.8580247E-6f
#define RM2				7.4373773E-6f

// 
enum
{
	// 共通項目

	// 整数
	SHARE_ISTATE,				// 
	SHARE_INT_MAX,

	// 少数
	SHARE_GSX = 0,				// m_gsx[i]
	SHARE_GSY,					// m_gsy[i]
	SHARE_ASY,					// bsy + tsim * m_spd + m_pstart;
	SHARE_TMPT,					// m_tmpt[0]
	SHARE_GV,					// sqrt(gvm / gvmi[indx]);
	SHARE_HHW,					// hhw[iht]
	SHARE_FLOAT_MAX,
};



///////////////////////////////////////////////////////////////////////////////////////////////
// LLG
__device__ void fun(const int ZCount, float &YDx, float &YDy, float &YDz, const float* hea, const float* Yx, const float* Yy, const float* Yz, const float &huk, const float &hax, const float &hay, const float &haz, const float &alfg, const float &htx, const float &hty, const float &htz, const float &amk)
{
	__syncthreads();

	// +1しているのはif文を減らすための工夫。
	const int z = threadIdx.x + 1;
	const int zn = z - 1;
	const int zp = z + 1;

	const float h = hea[z];
	const float hn = hea[zn];

	const float yx = Yx[z];
	const float yy = Yy[z];
	const float yz = Yz[z];

	const float yxn = Yx[zn];
	const float yyn = Yy[zn];
	const float yzn = Yz[zn];

	const float yxp = Yx[zp];
	const float yyp = Yy[zp];
	const float yzp = Yz[zp];

	// 本来ならここで一番上下の層はif文分けされる。
	// 配列の前後に余分な0の領域を持たしif文を減らしている。
	const float hex = (h * yxp) + (hn * yxn);
	const float hey = (h * yyp) + (hn * yyn);
	const float hez = (h * yzp) + (hn * yzn);

	const float hanz = huk * yz;

	const float ahx = hax + hex;
	const float ahy = hay + hey;
	const float ahz = hanz + haz + hez;

	const float agx = alfg * (ahx + htx);
	const float agy = alfg * (ahy + hty);
	const float agz = alfg * (ahz + htz);

	const float spx = (ahy * yz) - (ahz * yy);
	const float spy = (ahz * yx) - (ahx * yz);
	const float spz = (ahx * yy) - (ahy * yx);

	const float agm = (agx * yx) + (agy * yy) + (agz * yz);

	YDx = amk * (spx + agx - (agm * yx));
	YDy = amk * (spy + agy - (agm * yy));
	YDz = amk * (spz + agz - (agm * yz));
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 0～1までの乱数
__device__ float ran1(float* R, int &IX1, int &IX2, int &IX3)
{
	IX1 = (IA1 * IX1 + IC1) % M1;
	IX2 = (IA2 * IX2 + IC2) % M2;
	IX3 = (IA3 * IX3 + IC3) % M3;
	int J = ((RAND_ARRAY - 1) * IX3) / M3;

	float ran1 = R[J];
	R[J] = ((float)IX1 + (float)IX2 * RM2) * RM1;

	return ran1;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 大体-1から1までに収まる乱数
__device__ float gasdev(int &iset, float &gset, float* R, int &IX1, int &IX2, int &IX3)
{
	float gasdev = gset;

	if (iset == 0)
	{
		float v1 = 0.0, v2 = 0.0, rsq = 0.0;
		do
		{
			v1 = 2.0f * ran1(R, IX1, IX2, IX3) - 1.0f;
			v2 = 2.0f * ran1(R, IX1, IX2, IX3) - 1.0f;
			rsq = (v1 * v1) + (v2 * v2);
		} while ((rsq >= 1.0f) || (rsq == 0.0f));

		float fac = sqrt(-2.0f * log(rsq) / rsq);
		gset = v1 * fac;
		gasdev = v2 * fac;
	}

	iset = !iset;

	return gasdev;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 乱数の初期設定
__device__ void random(const int seed, float* R, int &IX1, int &IX2, int &IX3)
{
	IX1 = (IC1 - (-255 + seed)) % M1;
	IX1 = (IA1 * IX1 + IC1) % M1;
	IX2 = IX1 % M2;
	IX1 = (IA1 * IX1 + IC1) % M1;
	IX3 = IX1 % M3;

	for (int J = 0; J < RAND_ARRAY; J++)
	{
		IX1 = (IA1 * IX1 + IC1) % M1;
		IX2 = (IA2 * IX2 + IC2) % M2;
		R[J] = ((float)IX1 + (float)IX2 * RM2) * RM1;
	}

	ran1(R, IX1, IX2, IX3);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
// SINTDY computes interpolated values of the K-th derivative of the dependent variable vector y, and stores it in DKY.
// This routine is called within the package with K = 0 and T = TOUT, but may also be called by the user for any K up to the current order.
// (See detailed instructions in the usage documentation.)
// 
// The computed values in DKY are gotten by interpolation using the Nordsieck history array YH.
// This array corresponds uniquely to a vector-valued polynomial of degree NQCUR or less, and DKY is set to the K-th derivative of this polynomial at T.
// The formula for DKY is:
//             q
// DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
//             j=K
// where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
// 
// The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are communicated by COMMON.
// The above sum is done in reverse order.
// IFLAG is returned negative if either K or T is out of bounds.
// 
// SINTDY は従属変数ベクトル y の K 番目の導関数の補間値を計算し、DKY に格納する。
// このルーチンはパッケージ内で K = 0、T = TOUT で呼び出されますが、ユーザが現在の次数までの任意の K に対して呼び出すこともできます。
// (使用説明書の詳細な説明を参照)。
// 
// DKYの計算値は、Nordsieckの履歴配列YHを使った補間によって得られる。
// この配列は次数NQCUR以下のベクトル値多項式に一意に対応し、DKYはTにおけるこの多項式のK番目の導関数に設定される。
// DKYの式は
//             q
// DKY(i) = 和 c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
//             j=K
// ここで、c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCURである。
// 
// 量nq = NQCUR、l = nq+1、N = NEQ、tn、hはCOMMONで通信される。
// 上記の合計は逆順に行われる。
// IFLAGは、KまたはTのいずれかが範囲外の場合、負を返す。
// 
__device__ int sintdy(const int &ZCount, const int &L, const int &NQ, const float &T, const float &TN, const float &HU, const float &UROUND, const float &H, float &Yx, float &Yy, float &Yz, const float* YHx, const float* YHy, const float* YHz, const int &nYHIndex)
{
	// 元のコードから、Cは1固定のため計算式から削除している。

	float TP = TN - HU - 100.0f * UROUND * (TN + HU);

	if (((T - TP) * (T - TN)) > 0.0f)
	{
// 90
		return -2;
	}

	float S = (T - TN) / H;

// 15
	int nOffset = ZCount * L + nYHIndex;
	Yx = YHx[nOffset];
	Yy = YHy[nOffset];
	Yz = YHz[nOffset];

	int i;
	for (i = 0; i <= NQ; i++)
	{
// 35
		nOffset = ZCount * (NQ - i) + nYHIndex;
		Yx = YHx[nOffset] + (S * Yx);
		Yy = YHy[nOffset] + (S * Yy);
		Yz = YHz[nOffset] + (S * Yz);
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 重み付き二乗平均平方根ベクトル
__device__ float svnorm(const int &ZCount, float* Add, const float &Vx, const float &Vy, const float &Vz, const float &Wx, const float &Wy, const float &Wz)
{
	Add[threadIdx.x] = pow((Vx * Wx), 2.0f) + pow((Vy * Wy), 2.0f) + pow((Vz * Wz), 2.0f);

	// Warp内で値をまとめる
	int i;
	for (i = (blockDim.x >> 1); i > 0; i >>= 1)
	{
		__syncthreads();

		if (i > threadIdx.x)
			Add[threadIdx.x] += Add[threadIdx.x + i];
	}
	__syncthreads();

	// xyz分なのでZ*3している。
	if (0 == threadIdx.x)
		Add[0] = sqrt(Add[0] / (float)(ZCount * 3));
	__syncthreads();

	return Add[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// アニメーション用ログ
__device__ void animation_log(const int &Grains, const int &Grain, const int &ZCount, float* Add, const int &iht, const int &nmwskip, const int &npwskip, const float &Yx, const float &Yy, const float &Yz, const float &tmpt, const float &hfz, float *amdx, float *amdy, float *amdz, float *ztmp, float *d_zhed)
{
	int i;

#ifdef ALL_ANIMATION
	// Warp内で値をまとめる

	if (0 == threadIdx.x)
	{
		for (i = 0; i < MAX_Z; i++)
			Add[i] = 0.0f;
	}

	// 磁化 Z
	Add[threadIdx.x] = Yz;
	__syncthreads();
	for (i = (blockDim.x >> 1); i > 0; i >>= 1)
	{
		__syncthreads();

		if (i > threadIdx.x)
			Add[threadIdx.x] += Add[threadIdx.x + i];
	}
	__syncthreads();

	if (0 == threadIdx.x)
	{
		int nIndex = iht * Grains + Grain;
		amdz[nIndex] = Add[0] / (float)ZCount;

		// 熱
		ztmp[nIndex] = tmpt;

		// ヘッド
		d_zhed[nIndex] = hfz;
	}

#else
	int nIndex;
	if (0 == (iht % nmwskip))
	{
		// Warp内で値をまとめる
		nIndex = (iht / nmwskip) * Grains + Grain;

		if (0 == threadIdx.x)
		{
			for (i = 0; i < MAX_Z; i++)
				Add[i] = 0.0f;
		}

		// 磁化 Z
		Add[threadIdx.x] = Yz;
		__syncthreads();
		for (i = (blockDim.x >> 1); i > 0; i >>= 1)
		{
			__syncthreads();

			if (i > threadIdx.x)
				Add[threadIdx.x] += Add[threadIdx.x + i];
		}
		__syncthreads();
		if (0 == threadIdx.x)
			amdz[nIndex] = Add[0] / (float)ZCount;
	}

	// 熱
	if (0 == threadIdx.x)
	{
		if (0 == (iht % npwskip))
		{
			ztmp[iht / npwskip * Grains + Grain] = tmpt;
			d_zhed[iht / npwskip * Grains + Grain] = hfz;
		}
	}

#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// 計算のメイン
// 
// Fortranのgotoをフラグに置き変え処理している
// フラグの番号はgotoの飛び先の番号としている。
// 
// 計算失敗時は強制終了させる
// さらに磁化zに場所を表す異常値を入れている
// x、yにはその時の状況を見れる値を入れている
// 
// 高速化のため
// 共有メモリ多数使用してる (やりすぎると計算がおかしくなったので一部はグローバルメモリ)
// 関数はコールするのではなく、なるべく展開している
//
__global__ void hamrcd(const int ZCount, const int nht, float* d_ymx, float* d_ymy, float* d_ymz,
					   float* d_gsx, float* d_gsy, const int kxofs, float* d_ttw,
					   const int ntpx, const int ntpy, float* d_tpf,
					   const int nhfx, const int nhfy, const float nps, const float MP_edge, const float tpf_edge,
					   float* d_tc, const float alftscaling, const float dtc, const float tambient, float* d_gv,
					   float* d_htho, float* d_hhw, float* d_hfx, float* d_hfy, float* d_hfz, float* d_hko, float* d_hks, float* d_ept, float* d_eptAex, float* d_heb,
					   const float tf, const float dtf, const float rumach,
					   float* d_alfg, int* d_rs,
					   const float fOffset,
					   float* d_rand, float* a_YHx, float* a_YHy, float* a_YHz,
					   const int animation, const int nmwskip, const int npwskip, float* d_amdx, float* d_amdy, float* d_amdz, float* d_ztmp, float* d_zhed,
					   float* s_ELCO, float* s_TESCO,
#ifdef GPU_LOG
					   const int nLog, float* d_asx, float* d_asy, float* d_bsy, float* d_alft, float* d_amk, float* d_hhwl, float* d_tmpt, float* d_htz, float* d_hax, float* d_hay, float* d_haz, float* d_huk, float* d_hea, float* d_bms, float* d_yy, float* d_yz, float* d_amy, float* d_amz, float* d_tmpc, float* d_hthi,
#endif
					   const int Grains)
{
	// s_Yx、s_Yy、s_Yz、s_heaは挿入時にthreadIdx.xを+1している。
	// これはfunでの計算時にif分をなくすためである。
	// その事もあり、Zの最大数は「MAX_Z - 2」となる。

	__shared__ float d_YHx[MAX_Z * (MAXORD + 1)];
	__shared__ float d_YHy[MAX_Z * (MAXORD + 1)];
	__shared__ float d_YHz[MAX_Z * (MAXORD + 1)];

	__shared__ float s_Yx[MAX_Z];
	__shared__ float s_Yy[MAX_Z];
	__shared__ float s_Yz[MAX_Z];
	__shared__ float s_hea[MAX_Z];
	__shared__ int s_amk[MAX_Z];
	__shared__ float s_Add[MAX_Z];
	__shared__ float s_EL[SIZE_EL];
//	__shared__ float s_ELCO[SIZE_ELCO_D2 * SIZE_ELCO_D1];
//	__shared__ float s_TESCO[SIZE_TESCO_D2 * SIZE_TESCO_D1];

	__shared__ int n_Temp[SHARE_INT_MAX];
	__shared__ float f_Temp[SHARE_FLOAT_MAX];

	// 共有じゃないのにs_が付いているのは最初sharedにしていた名残。

#ifdef GPU_LOG
	const int Grain = nLog;
#else
	const int Grain = blockIdx.x;
#endif

	if (Grain < Grains)
	{
		if (threadIdx.x < ZCount)
		{
			int nOffset;

			int ksx, ksy, ksxh, ksyh;
			float tp = 0.0f, tn = 0.0f;

			// 乱数の格納にグローバルメモリを使っているので
			const int nRandIndex = (Grain * MAX_Z * RAND_ARRAY) + (RAND_ARRAY * threadIdx.x);
			int IX1 = 0, IX2 = 0, IX3 = 0;

			int iset = 0;
			float gset = 0.0f;

			float asy, tmpt;
			float htx, hty, htz, hax, hay, haz, huk, amk, hthi, hhw;
			float theta, phi, ssp, csp, sst, cst, bms, atmp;
			float tmptscling, alfbms;
			float fTemp;

			// s_取ると名前が被る。
			float s_alfg;

			float Yx, Yy, Yz;
			float EWTx, EWTy, EWTz;
			float SAVFx, SAVFy, SAVFz;
			float ACORx, ACORy, ACORz;

			float H0, TOL, SUM, W0;
			float TN, HU;
			int JSTART, NHNIL, NST, NQU;

			int NQ = 0, L = 0;
			float H;

			int IALTH;
			float RMAX, RC, EL0, CRATE, CONIT;

			int iht;
			float t, to;

			float yamp;

			int KFLAG;
			int IREDO, IRET, M, NCF, NEWQ;
			float DCON, DDN, DEL, DELP, DSM, DUP, EXSM, R, RH, RHDN, RHSM, RHUP, TOLD;

			int nCalculatingSLSODE;
			int nCalculatingSSTODE;

			int i, j;

			const int nIndex = Grain * ZCount + threadIdx.x;
			const int nYHIndex = threadIdx.x;

			// 出来る計算は先に行っておく
			const float div_1_3 = 1.0f / 3.0f;
			const float pi2 = 2.0f * PI;

			const float tmpc = d_tc[nIndex];
			const float htho = d_htho[threadIdx.x];
			float hfx = d_hfx[threadIdx.x];
			float hfy = d_hfy[threadIdx.x];
			float hfz = d_hfz[threadIdx.x];
			const float hko = d_hko[threadIdx.x];
			const float hks = d_hks[nIndex];
			const float ept = d_ept[threadIdx.x] / 3.0f;
			const float eptAex = (d_eptAex[threadIdx.x] - 1.0f) / 3.0f;
			const float heb = d_heb[threadIdx.x];
			const float alfg = d_alfg[threadIdx.x];

			Yx = d_ymx[nIndex];
			Yy = d_ymy[nIndex];
			Yz = d_ymz[nIndex];

			// 乱数の初期化
			nOffset = d_rs[nIndex];
			random(nOffset, &(d_rand[nRandIndex]), IX1, IX2, IX3);
			for (iht = 0; iht < nOffset; iht++)
				ran1(&(d_rand[nRandIndex]), IX1, IX2, IX3);
			gasdev(iset, gset, &(d_rand[nRandIndex]), IX1, IX2, IX3);

			// 座標
			if (0 == threadIdx.x)
			{
				// 丸めにより座標が変わってしまうため+1.0した後に-1してる。
				// 差し引き0だからと取ると座標が違ってくる。
				f_Temp[SHARE_GSX] = d_gsx[Grain] + fOffset;
				f_Temp[SHARE_GSY] = d_gsy[Grain];
				f_Temp[SHARE_GV] = d_gv[Grain];
			}
			__syncthreads();

			// X座標は0基準にされている。
			ksxh = (int)f_Temp[SHARE_GSX];

			ksx = ksxh - kxofs;

#if 1
// -----------------------------------------------------------------------
//  scfodeを直書き
// -----------------------------------------------------------------------
			if (0 == threadIdx.x)
			{
				s_ELCO[0] = 1.0f;
				s_ELCO[1] = 1.0f;
				s_TESCO[0] = 0.0f;
				s_TESCO[1] = 2.0f;
				s_TESCO[SIZE_TESCO_D1 * 1] = 1.0f;
				s_TESCO[SIZE_TESCO_D1 * 11 * 2] = 0.0f;

				float PC[12];
				PC[0] = 1.0f;

				float RQFAC = 1.0f;

				for (i = 1; i < 12; i++)
				{
					// -----------------------------------------------------------------------
					// The PC array will contain the coefficients of the polynomial
					//     p(x) = (x+1)*(x+2)*...*(x+nq-1).
					// Initially, p(x) = 1.
					// -----------------------------------------------------------------------
					float RQ1FAC = RQFAC;
					RQFAC /= (float)(i + 1);

					int NQM1 = i - 1;
					float FNQM1 = (float)i;

					int NQP1 = i + 1;

					// Form coefficients of p(x)*(x+nq-1). ----------------------------------

					PC[i] = 0.0f;
					for (j = 0; j <= NQM1; j++)
					{
						int nIndex = NQP1 - j - 1;
						PC[nIndex] = PC[nIndex - 1] + FNQM1 * PC[nIndex];
					}

					PC[0] *= FNQM1;

					// Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------

					float PINT = PC[0];
					float XPIN = PC[0] / 2.0f;
					float TSIGN = 1.0f;

					for(j = 1; j <= i; j++)
					{
						TSIGN = -TSIGN;
						PINT += (TSIGN * PC[j] / (float)(j + 1));
						XPIN += (TSIGN * PC[j] / (float)(j + 2));
					}

					// Store coefficients in ELCO and TESCO. --------------------------------

					s_ELCO[SIZE_ELCO_D1 * i] = PINT * RQ1FAC;
					s_ELCO[SIZE_ELCO_D1 * i + 1] = 1.0f;

					for (j = 1; j <= i; j++)
						s_ELCO[SIZE_ELCO_D1 * i + j + 1] = RQ1FAC * PC[j] / (float)(j + 1);

					float AGAMQ = RQFAC * XPIN;
					float RAGQ = 1.0f / AGAMQ;
					s_TESCO[SIZE_TESCO_D1 * i + 1] = RAGQ;

					if (NQP1 < 12)
						s_TESCO[SIZE_TESCO_D1 * NQP1] = RAGQ * RQFAC / (float)(NQP1 + 1);

					s_TESCO[SIZE_TESCO_D1 * NQM1 + 2] = RAGQ;
				}
			}
			__syncthreads();
#endif

			// 共有メモリの初期化
			// スレッド数が不確なため1つのスレッドで行う。
			if (0 == threadIdx.x)
			{
				for (i = 0; i < MAX_Z; i++)
				{
					s_Yx[i] = 0.0f;
					s_Yy[i] = 0.0f;
					s_Yz[i] = 0.0f;
					s_hea[i] = 0.0f;
					s_amk[i] = 0;
					s_Add[i] = 0.0f;
				}
			}
			__syncthreads();

			// 3Dヘッドを使う場合、初回だけ初期化する。
			if ((1 < nhfy) && (DUMMY_DOUBLE != nps))
				hfx = hfy = hfz = 0.0f;

			// 熱スポット移動
			for (iht = 0; iht < nht; iht++)
			{
				if (0 == threadIdx.x)
				{
					f_Temp[SHARE_HHW] = d_hhw[iht];
					f_Temp[SHARE_ASY] = d_ttw[iht] + f_Temp[SHARE_GSY];
				}
				__syncthreads();

				asy = f_Temp[SHARE_ASY];
				ksy = (int)asy;

				// caltemp
				// 初期温度
				tmpt = 275.0f;

				// 計算対象なら温度が入る。
				if ((asy >= 0.0f) && (ksy < ntpy) && (ksx >= 0) && (ksx < ntpx))
				{
					tmpt = d_tpf[(ZCount * ntpy * ksx) + (ZCount * ksy) + threadIdx.x];

					// 少数点下部分の補間。
					if (ksy > 0)
					{
						tp = tmpt;
						tn = d_tpf[(ZCount * ntpy * ksx) + (ZCount * (ksy - 1)) + threadIdx.x];
						tmpt = tn + ((tp - tn) * (asy - (float)ksy));
					}
				}

				// 3Dを使うならyは1より大きい。
				if ((1 < nhfy) && (DUMMY_DOUBLE != nps))
				{
					hfx = hfy = hfz = 0.0f;

					// 熱スポットより後方に配置される。
					//	y - nps - (MP_edge - tpf_edge)
					ksyh = ksy - (int)(nps - MP_edge - tpf_edge);

					if ((0 <= ksyh) && (ksyh < nhfy) && (0 <= ksxh) && (ksxh < nhfx))
					{
						hfx = d_hfx[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];
						hfy = d_hfy[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];
						hfz = d_hfz[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];

						// 少数点下部分の補間。
						// 熱の所と同じようになってもうまくいかなかったので自分がわかるやり方をした。
						// Y軸プラス方向があるか確認し補間の対象としている。
						ksyh++;
						if (nhfy > ksyh)
						{
							fTemp = asy - (float)ksyh + nps - MP_edge - tpf_edge;

							tn = hfx;
							tp = d_hfx[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];
							hfx = tn + ((tp - tn) * fTemp);

							tn = hfy;
							tp = d_hfy[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];
							hfy = tn + ((tp - tn) * fTemp);

							tn = hfz;
							tp = d_hfz[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];
							hfz = tn + ((tp - tn) * fTemp);
						}
					}
				}

				if (0 == threadIdx.x)
					f_Temp[SHARE_TMPT] = tmpt;

				__syncthreads();

				// 温度による計算の振り分け。

				htx = 0.0f;
				hty = 0.0f;
				htz = 0.0f;
				hax = 0.0f;
				hay = 0.0f;
				haz = 0.0f;
				huk = 0.0f;
				amk = 0.0f;
//				bms = 0.0f;
//				hthi = 0.0f;

//				s_hea[threadIdx.x + 1] = 0.0f;
//				s_alfg = alfg;

				s_amk[threadIdx.x] = 0;

				// 自信がキュリー温度を超えた
				if (tmpt > (tmpc - dtc))
				{
					theta = PI * ran1(&(d_rand[nRandIndex]), IX1, IX2, IX3);
					phi = pi2 * ran1(&(d_rand[nRandIndex]), IX1, IX2, IX3);
					ssp = sin(phi);
					csp = cos(phi);
					sst = sin(theta);
					cst = cos(theta);
					Yx = sst * csp;
					Yy = sst * ssp;
					Yz = cst;
				}
				// 一番上の層が指定温度より高い時
				else if (f_Temp[SHARE_TMPT] >= tambient)
				{
					tmptscling = alftscaling * tmpt;
					alfbms = pow((1.0f - tmptscling / tmpc), div_1_3);
					s_alfg = alfg / alfbms * (1.0f - tmptscling / (3.0f * tmpc));

					fTemp = 1.0f - (tmpt / tmpc);

					bms = pow(fTemp, div_1_3);
					hthi = htho * sqrt(s_alfg * tmpt / bms) * f_Temp[SHARE_GV];

					huk = hko * hks * pow(fTemp, ept);
					s_hea[threadIdx.x + 1] = heb * pow(fTemp, eptAex);

					atmp = 10.0f;
					while (abs(atmp) > 3.0f)
						atmp = gasdev(iset, gset, &(d_rand[nRandIndex]), IX1, IX2, IX3);

					theta = PI * ran1(&(d_rand[nRandIndex]), IX1, IX2, IX3);
					phi = pi2 * ran1(&(d_rand[nRandIndex]), IX1, IX2, IX3);

					ssp = sin(phi);
					csp = cos(phi);
					sst = sin(theta);
					cst = cos(theta);

					fTemp = hthi * atmp;
					htx = fTemp * sst * csp;
					hty = fTemp * sst * ssp;
					htz = fTemp * cst;

					hhw = f_Temp[SHARE_HHW];
					hax = hhw * hfx;
					hay = hhw * hfy;
					haz = hhw * hfz;

					amk = 1.0f;
					s_amk[threadIdx.x] = 1;
				}
				// 温度が計算範囲外
//				else
//					;

				// Warp内で値をまとめる
				for (i = (blockDim.x >> 1); i > 0; i >>= 1)
				{
					__syncthreads();

					if (i > threadIdx.x)
						s_amk[threadIdx.x] += s_amk[threadIdx.x + i];
				}
				__syncthreads();

				// 全層が計算結果を反映させないならこの位置はスキップする
#ifdef GPU_LOG
#else
				if (0 == s_amk[0])
				{
					if (1 == animation)
						animation_log(Grains, Grain, ZCount, s_Add, iht, nmwskip, npwskip, Yx, Yy, Yz, tmpt, hfz, d_amdx, d_amdy, d_amdz, d_ztmp, d_zhed);

					continue;
				}
#endif

				if (0 == threadIdx.x)
					n_Temp[SHARE_ISTATE] = 1;

				__syncthreads();

				// 時間毎の繰り返し
#ifdef GPU_LOG
				if (0 != s_amk[0])
#endif
				for (t = 0.0f; t < tf; t += dtf)
				{
					__syncthreads();

					to = t + dtf;

// -----------------------------------------------------------------------
//  slsodeを直書き
// -----------------------------------------------------------------------
					nCalculatingSLSODE = 1;

					// -----------------------------------------------------------------------
					// The following internal Common block contains variables which are
					// communicated between subroutines.  All real variables are listed
					// first, followed by all integers.  The block is declared in
					// Subroutines SLSODE, SINTDY, SSTODE, SPREPJ, and SSOLSY.
					// -----------------------------------------------------------------------
					// 以下の内部コモン・ブロックには、サブルーチン間で通信される変数が含まれている。サブルーチン間で通信されます。
					// すべての実数変数は 次にすべての整数が続く。
					// このブロックは サブルーチン SLSODE, SINTDY, SSTODE, SPREPJ, SSOLSY.
					// -----------------------------------------------------------------------
	 
					// -----------------------------------------------------------------------
					// Block A.
					// This code block is executed on every call.
					// It tests ISTATE and ITASK for legality and branches appropriately.
					// If ISTATE .GT. 1 but the flag INIT shows that initialization has
					// not yet been done, an error return occurs.
					// If ISTATE = 1 and TOUT = T, return immediately.
					// -----------------------------------------------------------------------
					// ブロックA。
					// このコードブロックは、呼び出しごとに実行される。
					// ISTATEとITASKの正当性をテストし、適切に分岐する。
					// もしISTATE .GT. 1であるが、フラグINITが初期化がまだ行われていないことを示す場合 1であるが、フラグINITが初期化がまだ行われていないことを示している場合、エラーリターンが発生する。
					// ISTATE = 1かつTOUT = Tの場合、直ちにリターンする。
					// -----------------------------------------------------------------------

					// ***FIRST EXECUTABLE STATEMENT  SLSODE

					if (1 == n_Temp[SHARE_ISTATE])
					{
						// -----------------------------------------------------------------------
						// Block B.
						// The next code block is executed for the initial call (ISTATE = 1),
						// or for a continuation call with parameter changes (ISTATE = 3).
						// It contains checking of all inputs and various initializations.
						// 
						// First check legality of the non-optional inputs NEQ, ITOL, IOPT,
						// MF, ML, and MU.
						// -----------------------------------------------------------------------
						// ブロックB。
						// 次のコードブロックは、最初の呼び出し(ISTATE = 1)に対して実行される、
						// またはパラメータの変更を伴う継続呼び出し(ISTATE = 3)に対して実行される。
						// これには、すべての入力のチェックとさまざまな初期化が含まれる。
						// 
						// 最初に、非オプション入力NEQ、ITOL、IOPTの合法性をチェックする、 MF、ML、およびMU。
						// -----------------------------------------------------------------------

						// Next process and check the optional inputs. --------------------------

						H0 = 0.0f;

						// -----------------------------------------------------------------------
						// Set work array pointers and check lengths LRW and LIW.
						// Pointers to segments of RWORK and IWORK are named by prefixing L to
						// the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
						// Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
						// -----------------------------------------------------------------------
						// ワーク配列のポインタを設定し、長さ LRW と LIW をチェックする。
						// RWORKとIWORKのセグメントへのポインタは、Lを先頭につけて命名される。を付けて命名する。 例えば、セグメントYHはRWORK(LYH)から始まる。
						// RWORKのセグメントは（順に）YH、WM、EWT、SAVF、ACORと表記される。
						// -----------------------------------------------------------------------

						// -----------------------------------------------------------------------
						// Block C.
						// The next block is for the initial call only (ISTATE = 1).
						// It contains all remaining initializations, the initial call to F,
						// and the calculation of the initial step size.
						// The error weights in EWT are inverted after being loaded.
						// -----------------------------------------------------------------------
						// Cブロック。
						// 次のブロックは、最初の呼び出しのみである(ISTATE = 1)。
						// このブロックには、残りのすべての初期化、F、 および初期ステップ・サイズの計算が含まれる。
						// EWTの誤差重みは、ロード後に反転される。
						// -----------------------------------------------------------------------

						TN = t;
						JSTART = 0;
						NHNIL = 0;
						NST = 0;
						HU = 0.0f;
						NQU = 0;

						// Initial call to F.  (LF0 points to YH(*,2).) -------------------------

						s_Yx[threadIdx.x + 1] = Yx;
						s_Yy[threadIdx.x + 1] = Yy;
						s_Yz[threadIdx.x + 1] = Yz;

						nOffset = ZCount + nYHIndex;
						fun(ZCount, d_YHx[nOffset], d_YHy[nOffset], d_YHz[nOffset], s_hea, s_Yx, s_Yy, s_Yz, huk, hax, hay, haz, s_alfg, htx, hty, htz, amk);

						// Load the initial value vector in YH. ---------------------------------

						d_YHx[nYHIndex] = Yx;
						d_YHy[nYHIndex] = Yy;
						d_YHz[nYHIndex] = Yz;

						// Load and invert the EWT array.  (H is temporarily set to 1.0.) -------

						NQ = 1;
						H = 1.0f;

						// sewset();
						EWTx = 1.0f / (RTO * abs(d_YHx[nYHIndex]) + ATO);
						EWTy = 1.0f / (RTO * abs(d_YHy[nYHIndex]) + ATO);
						EWTz = 1.0f / (RTO * abs(d_YHz[nYHIndex]) + ATO);

						// -----------------------------------------------------------------------
						// The coding below computes the step size, H0, to be attempted on the
						// first step, unless the user has supplied a value for this.
						// First check that TOUT - T differs significantly from zero.
						// A scalar tolerance quantity TOL is computed, as MAX(RTOL(I))
						// if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted
						// so as to be between 100*UROUND and 1.0E-3.
						// Then the computed value H0 is given by..
						//                                      NEQ
						//   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  )
						//                                       1
						// where   w0     = MAX ( ABS(T), ABS(TOUT) ),
						//         f(i)   = i-th component of initial value of f,
						//         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
						// The sign of H0 is inferred from the initial values of TOUT and T.
						// -----------------------------------------------------------------------
						// 以下のコーディングは、最初のステップで試行されるステップサイズH0を計算する。を計算する。
						// 最初に、TOUT - T がゼロと大きく異なることをチェックする。
						// スカラー許容量TOLが、MAX(RTOL(I)) が正であれば MAX(RTOL(I))として，そうでなければ MAX(ATOL(I)/ABS(Y(I)))として計算される。100*UROUNDと1.0E-3の間になるように調整される。
						// 計算値H0は次式で与えられる。
						// 										NEQ
						//   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2 )
						// 										 1
						// ここで、	w0		= MAX ( ABS(T), ABS(TOUT) )、
						// 			f(i)	= fの初期値のi番目の成分、
						// 			ywt(i)	= EWT(i)/TOL (y(i)の重み)。
						// H0の符号は、TOUTとTの初期値から推測される。
						// -----------------------------------------------------------------------

						if (H0 == 0.0f)
						{
							W0 = max(abs(t), abs(to));

							if (dtf < (2.0f * rumach * W0))
							{
// 計算失敗時の確認用
d_ymx[nIndex] = to;
d_ymy[nIndex] = W0;
d_ymz[nIndex] = 11111.0f;
								return;
//								break;
							}

							TOL = max(RTO, (100.0f * rumach));
							TOL = min(TOL, 0.001f);

							nOffset = ZCount + nYHIndex;
							SUM = svnorm(ZCount, s_Add, d_YHx[nOffset], d_YHy[nOffset], d_YHz[nOffset], EWTx, EWTy, EWTz);
							SUM = 1.0f / (TOL * W0 * W0) + TOL * pow(SUM, 2.0f);

							H0 = 1.0f / sqrt(SUM);
							H0 = min(H0, dtf);
							H0 = abs(H0) * ((0.0f <= dtf) ? 1.0f: -1.0f);
						}

						// Load H with H0 and scale YH(*,2) by H0. ------------------------------

						H = H0;

						nOffset = ZCount + nYHIndex;
						d_YHx[nOffset] *= H0;
						d_YHy[nOffset] *= H0;
						d_YHz[nOffset] *= H0;
					}
					else
					{
						// 2回目以降
						// -----------------------------------------------------------------------
						// Block D.
						// The next code block is for continuation calls only (ISTATE = 2 or 3)
						// and is to check stop conditions before taking a step.
						// -----------------------------------------------------------------------
						// Dブロック
						// 次のコードブロックは、継続コール（ISTATE = 2または3）専用である。で、ステップを実行する前に停止条件をチェックする。
						// -----------------------------------------------------------------------

						NST = 0;

						if (((TN - to) * H) >= 0.0f)
						{
							if (0 != sintdy(ZCount, L, NQ, to, TN, HU, rumach, H, Yx, Yy, Yz, d_YHx, d_YHy, d_YHz, nYHIndex))
							{
// 計算失敗時の確認用
d_ymx[nIndex] = H;
d_ymy[nIndex] = t;
d_ymz[nIndex] = 22222.0f;
								return;
//								break;
							}

							nCalculatingSLSODE = 0;
						}
					}

					while (nCalculatingSLSODE)
					{
						// -----------------------------------------------------------------------
						// Block E.
						// The next block is normally executed for all calls and contains
						// the call to the one-step core integrator SSTODE.
						//
						// This is a looping point for the integration steps.
						//
						// First check for too many steps being taken, update EWT (if not at
						// start of problem), check for too much accuracy being requested, and
						// check for H below the roundoff level in T.
						// -----------------------------------------------------------------------
						// ブロックE。
						// 次のブロックは、通常すべての呼び出しに対して実行され、次の内容を含む。ワンステップ・コア積分器SSTODEの呼び出しが含まれる。
						// これは統合ステップのループポイントである。
						// 最初にステップが多すぎないかチェックし、EWTを更新する。を更新し、要求される精度が高すぎないかチェックする。Tのラウンドオフ・レベル以下のHをチェックする。
						// -----------------------------------------------------------------------

// 250
						if (NST >= MXSTP0)
						{
// 計算失敗時の確認用
d_ymx[nIndex] = (float)NST;
d_ymy[nIndex] = t;
d_ymz[nIndex] = 33333.0f;
							return;
//							break;
						}

						//sewset();
						EWTx = 1.0f / (RTO * abs(d_YHx[nYHIndex]) + ATO);
						EWTy = 1.0f / (RTO * abs(d_YHy[nYHIndex]) + ATO);
						EWTz = 1.0f / (RTO * abs(d_YHz[nYHIndex]) + ATO);

// 270
						if (1.0f < (rumach * (fTemp = svnorm(ZCount, s_Add, d_YHx[nYHIndex], d_YHy[nYHIndex], d_YHz[nYHIndex], EWTx, EWTy, EWTz))))
						{

// 計算失敗時の確認用
d_ymx[nIndex] = fTemp;
d_ymy[nIndex] = t;
d_ymz[nIndex] = 44444.0f;
							return;
//							break;
						}

						if ((TN + H) == TN)
						{
							NHNIL++;
							if (NHNIL <= (MXHNL0 - 1))
							{
								//printf("SLSODE-  Warning..internal T (=R1) and H (=R2) are\n");
								//printf("      such that in the machine, T + H = T on the next step  \n");
								//printf("      (H = step size). Solver will continue anyway\n");

								if (NHNIL >= (MXHNL0 - 1))
								{
									//printf("SLSODE-  Above warning has been issued I1 times.  \n");
									//printf("      It will not be issued again for this problem\n");

// 計算失敗時の確認用
d_ymx[nIndex] = TN;
d_ymy[nIndex] = H;
d_ymz[nIndex] = 55555.0f;
return;							
								}
							}
						}

						// -----------------------------------------------------------------------
						//  CALL SSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,SPREPJ,SSOLSY)
						// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//  sstodeを直書き
// -----------------------------------------------------------------------
						nCalculatingSSTODE = 0;

						KFLAG = 0;
						TOLD = TN;
						NCF = 0;
						DELP = 0.0f;

						DEL = 0.0f;

						IRET = 0;

						if (JSTART == 0)
						{
							// -----------------------------------------------------------------------
							// On the first call, the order is set to 1, and other variables are
							// initialized.  RMAX is the maximully 1.E4 to compensate for the small
							// initial H, but then is normally equal to 10.  If a failure
							// occurs (in corrector convergence or error test), RMAX is set to 2
							// for the next increase.
							// -----------------------------------------------------------------------
							// 最初の呼び出しで、オーダーが1にセットされ、他の変数が初期化される。が初期化される。
							// RMAXは、1ステップでHを増加させることができる最大比率である。である。
							// 初期値は1.E4である。であるが、通常は10に等しい。
							// もし が発生した場合(コレクタの収束またはエラーテスト)、 RMAX は 2 に設定される。に設定される。
							// -----------------------------------------------------------------------

							NQ = 0;
							L = 1;
							IALTH = 2;
							RMAX = 10000.0f;
							RC = 0.0f;
							EL0 = 1.0f;
							CRATE = 0.7f;
							IRET = 3;

							// -----------------------------------------------------------------------
							// SCFODE is called to get all the integration coefficients for the
							// current METH.  Then the EL vector and related constants are reset
							// whenever the order NQ is changed, or at the start of the problem.
							// -----------------------------------------------------------------------
							// Hが変更される場合、H比RHがRMAX、HMIN、HMXIと照合され、YH配列が再スケーリングされる。RMAX，HMIN，HMXIと照合され，YH配列が再スケーリングされる。
							// IALTHは に設定される。に設定される。
							// -----------------------------------------------------------------------

					// 140

							nCalculatingSSTODE = 150;
						}

						// 700になると抜けられる
						while (700 != nCalculatingSSTODE)
						{
							__syncthreads();

							if (150 == nCalculatingSSTODE)
							{
					// 150
								if (threadIdx.x <= L)
									s_EL[threadIdx.x] = s_ELCO[SIZE_ELCO_D1 * NQ + threadIdx.x];
								__syncthreads();

								RC = RC * s_EL[0] / EL0;
								EL0 = s_EL[0];
								CONIT = 0.5f / (float)(NQ + 2 + 1);

								nCalculatingSSTODE = (2 == IRET) ? 170: 0;
							}

							if (170 == nCalculatingSSTODE)
							{
								// -----------------------------------------------------------------------
								// If H is being changed, the H ratio RH is checked against
								// RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
								// L = NQ + 1 to prevent a change of H for that many steps, unless
								// forced by a convergence or error test failure.
								// -----------------------------------------------------------------------
								// H が変更される場合、H 比 RH が RMAX、HMIN、HMXI と照合され、YH 配列が再スケーリングされる。
								// IALTH は L = NQ + 1 に設定され、収束またはエラーテストの失敗によって強制されない限り、そのステップ数分の H の変更を防止する。
								// -----------------------------------------------------------------------

					// 170
								RH = max(RH, 0.0f);

					// 175
								RH = min(RH, RMAX);
								//RH = RH / max(1.0f, (abs(H) * 0.0f * RH));	1で割るので無駄
								R = 1.0;

								for (i = 1; i <= L; i++)
								{
									R *= RH;

									nOffset = ZCount * i + nYHIndex;
									d_YHx[nOffset] *= R;
									d_YHy[nOffset] *= R;
									d_YHz[nOffset] *= R;
								}

								H *= RH;
								RC *= RH;
								IALTH = L + 1;

								if (IREDO == 0)
								{
					// 690
									RMAX = 10.0f;
									nCalculatingSSTODE = 700;
									continue;
								}

								nCalculatingSSTODE = 0;
							}

							// -----------------------------------------------------------------------
							// This section computes the predicted values by effectively
							// multiplying the YH array by the Pascal Triangle matrix.
							// RC is the ratio of new to old values of the coefficient  H*EL(1).
							// When RC differs from 1 by more than CCMAX, IPUP is set to MITER
							// to force PJAC to be called, if a Jacobian is involved.
							// In any case, PJAC is called at least every MSBP steps.
							// -----------------------------------------------------------------------
							// このセクションでは、YH 配列にパスカル・トライアングル行列を効果的に掛け合わせることにより、予測値を計算する。YH 配列に Pascal Triangle 行列を乗算する。
							// RC は係数 H*EL(1) の新しい値と古い値の比である。
							// RC が 1 と CCMAX 以上異なる場合、IPUP は MITER に設定され、PJAC が強制的に呼び出されます。に設定され、ヤコビアンが関与する場合はPJACが強制的に呼び出される。
							// いずれにせよ、PJAC は少なくとも MSBP ステップごとに呼び出される。
							// -----------------------------------------------------------------------
					// 200
							TN += H;

							for (i = 0; i <= NQ; i++)
							{
								for (j = (NQ - i); j <= NQ; j++)
								{
									nOffset = ZCount * (j + 1) + nYHIndex;
									d_YHx[nOffset - ZCount] += d_YHx[nOffset];
									d_YHy[nOffset - ZCount] += d_YHy[nOffset];
									d_YHz[nOffset - ZCount] += d_YHz[nOffset];
								}
							}

							// -----------------------------------------------------------------------
							// Up to MAXCOR corrector iterations are taken.  A convergence test is
							// made on the R.M.S. norm of each correction, weighted by the error
							// weight vector EWT.  The sum of the corrections is accumulated in the
							// vector ACOR(i).  The YH array is not altered in the corrector loop.
							// -----------------------------------------------------------------------
							// 最大 MAXCOR の補正の反復が行われる。
							// 収束テストは 各補正のR.M.S.ノルムについて収束テストが行われる。重みベクトル EWT によって重み付けされます。
							// 補正量の合計は ベクトルACOR(i)に累積される。
							// YH 配列は補正ループの中では変更されない。
							// -----------------------------------------------------------------------
					// 220
							M = 0;
							Yx = d_YHx[nYHIndex];
							Yy = d_YHy[nYHIndex];
							Yz = d_YHz[nYHIndex];

							s_Yx[threadIdx.x + 1] = Yx;
							s_Yy[threadIdx.x + 1] = Yy;
							s_Yz[threadIdx.x + 1] = Yz;
							fun(ZCount, SAVFx, SAVFy, SAVFz, s_hea, s_Yx, s_Yy, s_Yz, huk, hax, hay, haz, s_alfg, htx, hty, htz, amk);

							// -----------------------------------------------------------------------
							// If indicated, the matrix P = I - h*el(1)*J is reevaluated and
							// preprocessed before starting the corrector iteration.  IPUP is set
							// to 0 as an indicator that this has been done.
							// -----------------------------------------------------------------------
							// 指示された場合、P = I - h*el(1)*J行列が再評価され、補正器の反復を開始する前に前処理が行われる。 が再評価され, 補正器の反復を開始する前に前処理される.
							// IPUP は に 0 がセットされる。
							// -----------------------------------------------------------------------

							ACORx = 0.0;
							ACORy = 0.0;
							ACORz = 0.0;

							// 次に進めるなら0ではなくなる
							while (0 == nCalculatingSSTODE)
							{
								__syncthreads();

								//-----------------------------------------------------------------------
								// In the case of functional iteration, update Y directly from
								// the result of the last function evaluation.
								// -----------------------------------------------------------------------
								// 関数反復の場合、最後の関数評価結果Cから直接Yを更新する。を直接更新します。
								// -----------------------------------------------------------------------
					// 270
								nOffset = ZCount + nYHIndex;
								SAVFx = H * SAVFx - d_YHx[nOffset];
								SAVFy = H * SAVFy - d_YHy[nOffset];
								SAVFz = H * SAVFz - d_YHz[nOffset];

								Yx = SAVFx - ACORx;
								Yy = SAVFy - ACORy;
								Yz = SAVFz - ACORz;

								DEL = svnorm(ZCount, s_Add, Yx, Yy, Yz, EWTx, EWTy, EWTz);

								Yx = d_YHx[nYHIndex] + EL0 * SAVFx;
								Yy = d_YHy[nYHIndex] + EL0 * SAVFy;
								Yz = d_YHz[nYHIndex] + EL0 * SAVFz;

								ACORx = SAVFx;
								ACORy = SAVFy;
								ACORz = SAVFz;

								// -----------------------------------------------------------------------
								// Test for convergence.  If M.gt.0, an estimate of the convergence
								// rate constant is stored in CRATE, and this is used in the test.
								// -----------------------------------------------------------------------
								// 収束のテスト。
								// M.gt.0の場合、収束率定数の推定値がCRATEに格納され、テストに使用される。の推定値が CRATE に格納され、これがテストに使用される。
								// -----------------------------------------------------------------------
					// 400
								if (M != 0)
									CRATE = max((0.2f * CRATE), (DEL / DELP));

								DCON = DEL * min(1.0f, (1.5f * CRATE)) / (s_TESCO[SIZE_TESCO_D1 * NQ + 1] * CONIT);

								if (DCON > 1.0f)
								{
									M++;

									if (M >= MAXCOR)
									{
										nCalculatingSSTODE = 430;
										break;
									}

									if ((M > 1) && (DEL > (2.0f * DELP)))
									{
										nCalculatingSSTODE = 430;
										break;
									}

									DELP = DEL;

									s_Yx[threadIdx.x + 1] = Yx;
									s_Yy[threadIdx.x + 1] = Yy;
									s_Yz[threadIdx.x + 1] = Yz;
									fun(ZCount, SAVFx, SAVFy, SAVFz, s_hea, s_Yx, s_Yy, s_Yz, huk, hax, hay, haz, s_alfg, htx, hty, htz, amk);

									// -> 270
									continue;
								}

								// -----------------------------------------------------------------------
								// The corrector has converged.  JCUR is set to 0
								// to signal that the Jacobian involved may need updating later.
								// The local error test is made and control passes to statement 500
								// if it fails.
								// -----------------------------------------------------------------------
								// 補正器は収束した。
								// JCUR は 0 にセットされる。
								// ローカル・エラー・テストが行われ、失敗した場合はステートメント500 に渡される。
								// -----------------------------------------------------------------------
					// 450
								if (0 == M)
									DSM = DEL / s_TESCO[SIZE_TESCO_D1 * NQ + 1];
								else
									DSM = svnorm(ZCount, s_Add, ACORx, ACORy, ACORz, EWTx, EWTy, EWTz) / s_TESCO[SIZE_TESCO_D1 * NQ + 1];

								if (1.0f >= DSM)
								{
									// -----------------------------------------------------------------------
									// After a successful step, update the YH array.
									// Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
									// If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
									// use in a possible order increase on the next step.
									// If a change in H is considered, an increase or decrease in order
									// by one is considered also.  A change in H is made only if it is by a
									// factor of at least 1.1.  If not, IALTH is set to 3 to prevent
									// testing for that many steps.
									// -----------------------------------------------------------------------
									// ステップ成功後、YH配列を更新する。
									// IALTH = 1であればHの変更を検討する。
									// そうでなければIALTHを1減らす。
									// IALTHが1であり、NQ .lt. MAXORDの場合、ACORは が保存される。
									// Hの変更が考慮される場合、オーダーCの1増減も考慮される。の増減も考慮される。
									// Hの変更は 倍以上である場合にのみ行われる。
									// そうでない場合、IALTH は 3 に設定され、その数ステップの C テストが行われないようにする。に設定される。
									// -----------------------------------------------------------------------

									KFLAG = 0;
									IREDO = 0;
									NST++;
									HU = H;
									NQU = NQ;

									for (i = 0; i <= L; i++)
									{
										nOffset = ZCount * i + nYHIndex;
										d_YHx[nOffset] += (s_EL[i] * ACORx);
										d_YHy[nOffset] += (s_EL[i] * ACORy);
										d_YHz[nOffset] += (s_EL[i] * ACORz);
									}

									IALTH--;

									if (0 == IALTH)
									{
					// 520
										RHUP = 0.0f;
										if (L != (MAXORD - 1))
										{
											nOffset = ZCount * (MAXORD - 1) + nYHIndex;
											SAVFx = ACORx - d_YHx[nOffset];
											SAVFy = ACORy - d_YHy[nOffset];
											SAVFz = ACORz - d_YHz[nOffset];

											DUP = svnorm(ZCount, s_Add, SAVFx, SAVFy, SAVFz, EWTx, EWTy, EWTz) / s_TESCO[SIZE_TESCO_D1 * NQ + 2];
											RHUP = 1.0f / (1.4f * pow(DUP, (1.0f / (L + 1 + 1))) + 0.0000014f);
										}
									}
									else
									{
										if ((IALTH <= 1) && (L != (MAXORD - 1)))
										{
											nOffset = ZCount * (MAXORD - 1) + nYHIndex;
											d_YHx[nOffset] = ACORx;
											d_YHy[nOffset] = ACORy;
											d_YHz[nOffset] = ACORz;
										}

										nCalculatingSSTODE = 700;
										break;
									}
								}
								else
								{
									// -----------------------------------------------------------------------
									// The error test failed.  KFLAG keeps track of multiple failures.
									// Restore TN and the YH array to their previous values, and prepare
									// to try the step again.  Compute the optimum step size for this or
									// one lower order.  After 2 or more failures, H is forced to decrease
									// by a factor of 0.2 or less.
									// -----------------------------------------------------------------------
									// エラーテストが失敗した。 KFLAGは複数の失敗を記録する。
									// TNとYHアレイを以前の値に戻し、ステップを再試行する準備をする。ステップを再試行する準備をする。
									// 最適なステップ・サイズを計算する。つ下の次数について最適なステップサイズを計算する。
									// 2回以上の失敗の後、Hは強制的に を0.2倍以下にする。
									// -----------------------------------------------------------------------

					// 500
									KFLAG--;
									TN = TOLD;

									for (i = 0; i <= NQ; i++)
									{
										for (j = (NQ - i); j <= NQ; j++)
										{
											nOffset = ZCount * (j + 1) + nYHIndex;
											d_YHx[nOffset - ZCount] -= d_YHx[nOffset];
											d_YHy[nOffset - ZCount] -= d_YHy[nOffset];
											d_YHz[nOffset - ZCount] -= d_YHz[nOffset];
										}
									}

									RMAX = 2.0f;
									if (abs(H) <= 0.0f)
									{
										KFLAG = -11;
										return;
									}

									if (KFLAG <= -3)
									{
										// -----------------------------------------------------------------------
										// Control reaches this section if 3 or more failures have occured.
										// If 10 failures have occurred, exit with KFLAG = -1.
										// It is assumed that the derivatives that have accumulated in the
										// YH array have errors of the wrong order.  Hence the first
										// derivative is recomputed, and the order is set to 1.  Then
										// H is reduced by a factor of 10, and the step is retried,
										// until it succeeds or H reaches HMIN.
										// -----------------------------------------------------------------------
										// 失敗が3回以上発生した場合、制御はこのセクションに到達する。
										// 失敗が 10 回発生した場合は、KFLAG = -1 で終了する。
										// YH配列に蓄積された導関数は、間違った順序のエラーを持っていると仮定される。
										// したがって、最初の導関数が再計算され、順序が1に設定される。
										// そして Hを10倍し、ステップを再試行する、 成功するかHがHMINに達するまで。
										// -----------------------------------------------------------------------
					// 640
										if (KFLAG <= -10)
										{
											KFLAG = -10;
											return;
										}

										RH = 0.1f;
										RH = max(0.0f, RH);
										H *= RH;

										Yx = d_YHx[nYHIndex];
										Yy = d_YHy[nYHIndex];
										Yz = d_YHz[nYHIndex];

										s_Yx[threadIdx.x + 1] = Yx;
										s_Yy[threadIdx.x + 1] = Yy;
										s_Yz[threadIdx.x + 1] = Yz;
										fun(ZCount, SAVFx, SAVFy, SAVFz, s_hea, s_Yx, s_Yy, s_Yz, huk, hax, hay, haz, s_alfg, htx, hty, htz, amk);

					// 650
										nOffset = ZCount + nYHIndex;
										d_YHx[nOffset] = H * SAVFx;
										d_YHy[nOffset] = H * SAVFy;
										d_YHz[nOffset] = H * SAVFz;

										IALTH = 5;

										if (0 == NQ)
										{
											// Next 200
											nCalculatingSSTODE = 0;
											break;
										}

										NQ = 0;
										L = 1;
										IRET = 3;

										nCalculatingSSTODE = 150;
										break;
									}

									IREDO = 2;
									RHUP = 0.0f;

									// -> 540
								}

								// -----------------------------------------------------------------------
								// Regardless of the success or failure of the step, factors RHDN, RHSM, and RHUP are computed, by which H could be multiplied at order NQ - 1, order NQ, or order NQ + 1, respectively.
								// In the case of failure, RHUP = 0.0 to avoid an order increase.
								// The largest of these is determined and the new order chosen 	accordingly.
								// If the order is to be increased, we compute one additional scaled derivative.
								// -----------------------------------------------------------------------
								// ステップの成否にかかわらず、係数 RHDN、RHSM、RHUPが計算され、それぞれHが次数NQ - 1または次数NQ + 1で を乗じることができる。
								// 失敗の場合、次数の増加を避けるためにRHUP=0.0とする。
								// これらのうち最大のものが決定され、それに応じて新しいオーダーが選択される。を選択する。
								// 次数を増加させる場合は、次のように計算する。を計算する。
								// -----------------------------------------------------------------------
					// 540
								EXSM = 1.0f / (float)(L + 1);
								RHSM = 1.0f / (1.2f * pow(DSM, EXSM) + 0.0000012f);
								RHDN = 0.0f;

								if (0 != NQ)
								{
									nOffset = ZCount * L + nYHIndex;
									DDN = svnorm(ZCount, s_Add, d_YHx[nOffset], d_YHy[nOffset], d_YHz[nOffset], EWTx, EWTy, EWTz) / s_TESCO[SIZE_TESCO_D1 * NQ];

									RHDN = 1.0f / (1.3f * pow(DDN, (1.0f / (float)(NQ + 1))) + 0.0000013f);
								}

					// 560
								if ((RHSM < RHUP) && (RHUP > RHDN))
								{
					// 590
									NEWQ = L;
									RH = RHUP;

									if (RH < 1.1f)
									{
					// 610
										IALTH = 3;
										nCalculatingSSTODE = 700;
										break;
									}

									R = s_EL[L] / (float)(L + 1);

									nOffset = ZCount * (NEWQ + 1) + nYHIndex;
									d_YHx[nOffset] = ACORx * R;
									d_YHy[nOffset] = ACORy * R;
									d_YHz[nOffset] = ACORz * R;
								}
								else
								{
					// 570
									if (RHSM >= RHUP)
									{
										NEWQ = NQ;
										RH = RHSM;
									}
									else
									{
					// 580
										NEWQ = NQ - 1;
										RH = RHDN;

										if ((KFLAG < 0) && (RH > 1.0f))
											RH = 1.0;
									}

					// 620
									if ((KFLAG == 0) && (RH < 1.1f))
									{
										IALTH = 3;
										nCalculatingSSTODE = 700;
										break;
									}

									if (KFLAG <= -2)
										RH = min(RH, 0.2f);

									// -----------------------------------------------------------------------
									// If there is a change of order, reset NQ, l, and the coefficients.
									// In any case H is reset according to RH and the YH array is rescaled.
									// Then exit from 690 if the step was OK, or redo the step otherwise.
									// -----------------------------------------------------------------------
									// オーダーが変更された場合は、NQ、l、係数をリセットする。
									// いずれの場合も、H は RH に従ってリセットされ、YH 配列は再スケーリングされる。
									// ステップに問題がなければ690を終了し、そうでなければステップをやり直す。
									// -----------------------------------------------------------------------

									if (NEWQ == NQ)
									{
										nCalculatingSSTODE = 170;
										break;
									}
								}

					// 630
								NQ = NEWQ;
								L = NQ + 1;
								IRET = 2;

								nCalculatingSSTODE = 150;
							}

							__syncthreads();

							switch (nCalculatingSSTODE)
							{
							case 430:		break;
							default:		continue;
							}

							// -----------------------------------------------------------------------
							// The corrector iteration failed to converge.
							// If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
							// the next try.  Otherwise the YH array is retracted to its values
							// before prediction, and H is reduced, if possible.  If H cannot be
							// reduced or MXNCF failures have occurred, exit with KFLAG = -2.
							// -----------------------------------------------------------------------
							// 補正器の反復が収束しなかった。
							// もし MITER .ne. 0 でヤコビアンが古い場合、PJAC が呼ばれる。が呼び出される。
							// そうでなければ、YH 配列は予測前の値に戻される。に戻され、可能であれば H が減らされる。
							// Hが縮小できない場合 Hが縮小できないか、MXNCFに失敗した場合は、KFLAG = -2で終了する。
							// -----------------------------------------------------------------------

					// 430
							NCF++;
							RMAX = 2.0f;
							TN = TOLD;

							for (i = 0; i <= NQ; i++)
							{
								for (j = (NQ - i); j <= NQ; j++)
								{
									nOffset = ZCount * (j + 1) + nYHIndex;
									d_YHx[nOffset - ZCount] -= d_YHx[nOffset];
									d_YHy[nOffset - ZCount] -= d_YHy[nOffset];
									d_YHz[nOffset - ZCount] -= d_YHz[nOffset];
								}
							}

					// 445
							if (abs(H) < 0.0f)
							{
								KFLAG = -12;
								return;
							}

							if (NCF == MXNCF)
							{
								KFLAG = -13;
								return;
							}

							RH = 0.25f;
							IREDO = 1;

							nCalculatingSSTODE = 170;
						}

					// 700
						R = 1.0f / s_TESCO[SIZE_TESCO_D1 * NQU + 1];

						ACORx *= R;
						ACORy *= R;
						ACORz *= R;

					// 720
						JSTART = 1;
// -----------------------------------------------------------------------
//  sstodeの終わり
// -----------------------------------------------------------------------

						// sstodeでのエラーチェック
						if (1 != (1 - KFLAG))
						{
// 計算失敗時の確認用
d_ymx[nIndex] = (float)KFLAG;
d_ymy[nIndex] = H;
d_ymz[nIndex] = 66666.0f;
							return;
//							break;
						}

						// 収束の異常判定
						if (0.0f == H)
						{
// 計算失敗時の確認用
d_ymx[nIndex] = H;
d_ymy[nIndex] = t;
d_ymz[nIndex] = 77777.0f;
							return;
						}

						// -----------------------------------------------------------------------
						// Block F.
						// The following block handles the case of a successful return from the
						// core integrator (KFLAG = 0).  Test for stop conditions.
						// -----------------------------------------------------------------------
						// Fブロック
						// 次のブロックは、コアインテグレータからのリターン成功（KFLAG = 0）の場合を処理する。(KFLAG=0)。 停止条件をテストする。
						// -----------------------------------------------------------------------
// 300
						// ITASK = 1.  If TOUT has been reached, interpolate. -------------------
						if (((TN - to) * H) < 0.0)
						{
							continue;
						}

						if (0 != sintdy(ZCount, L, NQ, to, TN, HU, rumach, H, Yx, Yy, Yz, d_YHx, d_YHy, d_YHz, nYHIndex))
						{
// 計算失敗時の確認用
d_ymx[nIndex] = TN;
d_ymy[nIndex] = H;
d_ymz[nIndex] = 88888.0f;
							return;
//							break;
						}


						nCalculatingSLSODE = 0;
					}

					// -----------------------------------------------------------------------
					// Block G.
					// The following block handles all successful returns from SLSODE.
					// If ITASK .NE. 1, Y is loaded from YH and T is set accordingly.
					// ISTATE is set to 2, and the optional outputs are loaded into the
					// work arrays before returning.
					// -----------------------------------------------------------------------
					// Gブロック
					// 次のブロックは、SLSODEからのすべての成功リターンを処理する。
					// もしITASK .NE. 1 の場合、YH から Y がロードされ、それに応じて T が設定される。
					// ISTATE が 2 に設定され、オプション出力がワーク・アレイにロードされる。
					// -----------------------------------------------------------------------
// 420
					__syncthreads();
					if (0 == threadIdx.x)
						n_Temp[SHARE_ISTATE] = 2;

// -----------------------------------------------------------------------
//  slsodeの終わり
// -----------------------------------------------------------------------

					yamp = sqrt((Yx * Yx) + (Yy * Yy) + (Yz * Yz));

					// 1つでもおかしい層があれば全体をはじめからやり直します。 (ただし時間は巻き戻さない)
					if (0.02f < abs(yamp - 1.0f))
						n_Temp[SHARE_ISTATE] = 1;

					Yx /= yamp;
					Yy /= yamp;
					Yz /= yamp;
				}

				// アニメーション用ログ
				if (1 == animation)
					animation_log(Grains, Grain, ZCount, s_Add, iht, nmwskip, npwskip, Yx, Yy, Yz, tmpt, hfz, d_amdx, d_amdy, d_amdz, d_ztmp, d_zhed);

#ifdef GPU_LOG
				__syncthreads();
				if (f_Temp[SHARE_TMPT] > tambient)
				{
					int nTemp = iht * ZCount + threadIdx.x;
					d_asx[nTemp] = f_Temp[SHARE_GSX];
					d_asy[nTemp] = asy;
					d_bsy[nTemp] = f_Temp[SHARE_GSY];
					d_alft[nTemp] = s_alfg;
					d_amk[nTemp] = amk;
					d_hhwl[nTemp] = hhw;
//					d_hhwl[nTemp] = f_Temp[SHARE_HHW];
					d_tmpt[nTemp] = tmpt;
//					d_tmpt[nTemp] = f_Temp[SHARE_TMPT];
					d_htz[nTemp] = htz;
					d_hax[nTemp] = hfx;
					d_hay[nTemp] = hfy;
					d_haz[nTemp] = hfz;
//					d_haz[nTemp] = haz;
					d_huk[nTemp] = huk;
					d_hea[nTemp] = s_hea[threadIdx.x + 1];
					d_bms[nTemp] = bms;
					d_yy[nTemp] = s_Yy[threadIdx.x + 1];
					d_yz[nTemp] = s_Yz[threadIdx.x + 1];
					d_amy[nTemp] = 0.0f;
					d_amz[nTemp] = 0.0f;
					d_tmpc[nTemp] = tmpc;
					d_hthi[nTemp] = hthi;
				}
				__syncthreads();
#endif
			}

			d_ymx[nIndex] = Yx;
			d_ymy[nIndex] = Yy;
			d_ymz[nIndex] = Yz;
		}
	}
}


