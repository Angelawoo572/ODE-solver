#include "CUDAFunction.cuh"

// �����p
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
	// ���ʍ���

	// ����
	SHARE_ISTATE,				// 
	SHARE_INT_MAX,

	// ����
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

	// +1���Ă���̂�if�������炷���߂̍H�v�B
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

	// �{���Ȃ炱���ň�ԏ㉺�̑w��if�����������B
	// �z��̑O��ɗ]����0�̗̈��������if�������炵�Ă���B
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
// 0�`1�܂ł̗���
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
// ���-1����1�܂łɎ��܂闐��
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
// �����̏����ݒ�
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
// SINTDY �͏]���ϐ��x�N�g�� y �� K �Ԗڂ̓��֐��̕�Ԓl���v�Z���ADKY �Ɋi�[����B
// ���̃��[�`���̓p�b�P�[�W���� K = 0�AT = TOUT �ŌĂяo����܂����A���[�U�����݂̎����܂ł̔C�ӂ� K �ɑ΂��ČĂяo�����Ƃ��ł��܂��B
// (�g�p�������̏ڍׂȐ������Q��)�B
// 
// DKY�̌v�Z�l�́ANordsieck�̗���z��YH���g������Ԃɂ���ē�����B
// ���̔z��͎���NQCUR�ȉ��̃x�N�g���l�������Ɉ�ӂɑΉ����ADKY��T�ɂ����邱�̑�������K�Ԗڂ̓��֐��ɐݒ肳���B
// DKY�̎���
//             q
// DKY(i) = �a c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
//             j=K
// �����ŁAc(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR�ł���B
// 
// ��nq = NQCUR�Al = nq+1�AN = NEQ�Atn�Ah��COMMON�ŒʐM�����B
// ��L�̍��v�͋t���ɍs����B
// IFLAG�́AK�܂���T�̂����ꂩ���͈͊O�̏ꍇ�A����Ԃ��B
// 
__device__ int sintdy(const int &ZCount, const int &L, const int &NQ, const float &T, const float &TN, const float &HU, const float &UROUND, const float &H, float &Yx, float &Yy, float &Yz, const float* YHx, const float* YHy, const float* YHz, const int &nYHIndex)
{
	// ���̃R�[�h����AC��1�Œ�̂��ߌv�Z������폜���Ă���B

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
// �d�ݕt����敽�ϕ������x�N�g��
__device__ float svnorm(const int &ZCount, float* Add, const float &Vx, const float &Vy, const float &Vz, const float &Wx, const float &Wy, const float &Wz)
{
	Add[threadIdx.x] = pow((Vx * Wx), 2.0f) + pow((Vy * Wy), 2.0f) + pow((Vz * Wz), 2.0f);

	// Warp���Œl���܂Ƃ߂�
	int i;
	for (i = (blockDim.x >> 1); i > 0; i >>= 1)
	{
		__syncthreads();

		if (i > threadIdx.x)
			Add[threadIdx.x] += Add[threadIdx.x + i];
	}
	__syncthreads();

	// xyz���Ȃ̂�Z*3���Ă���B
	if (0 == threadIdx.x)
		Add[0] = sqrt(Add[0] / (float)(ZCount * 3));
	__syncthreads();

	return Add[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// �A�j���[�V�����p���O
__device__ void animation_log(const int &Grains, const int &Grain, const int &ZCount, float* Add, const int &iht, const int &nmwskip, const int &npwskip, const float &Yx, const float &Yy, const float &Yz, const float &tmpt, const float &hfz, float *amdx, float *amdy, float *amdz, float *ztmp, float *d_zhed)
{
	int i;

#ifdef ALL_ANIMATION
	// Warp���Œl���܂Ƃ߂�

	if (0 == threadIdx.x)
	{
		for (i = 0; i < MAX_Z; i++)
			Add[i] = 0.0f;
	}

	// ���� Z
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

		// �M
		ztmp[nIndex] = tmpt;

		// �w�b�h
		d_zhed[nIndex] = hfz;
	}

#else
	int nIndex;
	if (0 == (iht % nmwskip))
	{
		// Warp���Œl���܂Ƃ߂�
		nIndex = (iht / nmwskip) * Grains + Grain;

		if (0 == threadIdx.x)
		{
			for (i = 0; i < MAX_Z; i++)
				Add[i] = 0.0f;
		}

		// ���� Z
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

	// �M
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
// �v�Z�̃��C��
// 
// Fortran��goto���t���O�ɒu���ς��������Ă���
// �t���O�̔ԍ���goto�̔�ѐ�̔ԍ��Ƃ��Ă���B
// 
// �v�Z���s���͋����I��������
// ����Ɏ���z�ɏꏊ��\���ُ�l�����Ă���
// x�Ay�ɂ͂��̎��̏󋵂������l�����Ă���
// 
// �������̂���
// ���L�����������g�p���Ă� (��肷����ƌv�Z�����������Ȃ����̂ňꕔ�̓O���[�o��������)
// �֐��̓R�[������̂ł͂Ȃ��A�Ȃ�ׂ��W�J���Ă���
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
	// s_Yx�As_Yy�As_Yz�As_hea�͑}������threadIdx.x��+1���Ă���B
	// �����fun�ł̌v�Z����if�����Ȃ������߂ł���B
	// ���̎�������AZ�̍ő吔�́uMAX_Z - 2�v�ƂȂ�B

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

	// ���L����Ȃ��̂�s_���t���Ă���͍̂ŏ�shared�ɂ��Ă������c�B

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

			// �����̊i�[�ɃO���[�o�����������g���Ă���̂�
			const int nRandIndex = (Grain * MAX_Z * RAND_ARRAY) + (RAND_ARRAY * threadIdx.x);
			int IX1 = 0, IX2 = 0, IX3 = 0;

			int iset = 0;
			float gset = 0.0f;

			float asy, tmpt;
			float htx, hty, htz, hax, hay, haz, huk, amk, hthi, hhw;
			float theta, phi, ssp, csp, sst, cst, bms, atmp;
			float tmptscling, alfbms;
			float fTemp;

			// s_���Ɩ��O�����B
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

			// �o����v�Z�͐�ɍs���Ă���
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

			// �����̏�����
			nOffset = d_rs[nIndex];
			random(nOffset, &(d_rand[nRandIndex]), IX1, IX2, IX3);
			for (iht = 0; iht < nOffset; iht++)
				ran1(&(d_rand[nRandIndex]), IX1, IX2, IX3);
			gasdev(iset, gset, &(d_rand[nRandIndex]), IX1, IX2, IX3);

			// ���W
			if (0 == threadIdx.x)
			{
				// �ۂ߂ɂ����W���ς���Ă��܂�����+1.0�������-1���Ă�B
				// ��������0������Ǝ��ƍ��W������Ă���B
				f_Temp[SHARE_GSX] = d_gsx[Grain] + fOffset;
				f_Temp[SHARE_GSY] = d_gsy[Grain];
				f_Temp[SHARE_GV] = d_gv[Grain];
			}
			__syncthreads();

			// X���W��0��ɂ���Ă���B
			ksxh = (int)f_Temp[SHARE_GSX];

			ksx = ksxh - kxofs;

#if 1
// -----------------------------------------------------------------------
//  scfode�𒼏���
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

			// ���L�������̏�����
			// �X���b�h�����s�m�Ȃ���1�̃X���b�h�ōs���B
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

			// 3D�w�b�h���g���ꍇ�A���񂾂�����������B
			if ((1 < nhfy) && (DUMMY_DOUBLE != nps))
				hfx = hfy = hfz = 0.0f;

			// �M�X�|�b�g�ړ�
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
				// �������x
				tmpt = 275.0f;

				// �v�Z�ΏۂȂ牷�x������B
				if ((asy >= 0.0f) && (ksy < ntpy) && (ksx >= 0) && (ksx < ntpx))
				{
					tmpt = d_tpf[(ZCount * ntpy * ksx) + (ZCount * ksy) + threadIdx.x];

					// �����_�������̕�ԁB
					if (ksy > 0)
					{
						tp = tmpt;
						tn = d_tpf[(ZCount * ntpy * ksx) + (ZCount * (ksy - 1)) + threadIdx.x];
						tmpt = tn + ((tp - tn) * (asy - (float)ksy));
					}
				}

				// 3D���g���Ȃ�y��1���傫���B
				if ((1 < nhfy) && (DUMMY_DOUBLE != nps))
				{
					hfx = hfy = hfz = 0.0f;

					// �M�X�|�b�g������ɔz�u�����B
					//	y - nps - (MP_edge - tpf_edge)
					ksyh = ksy - (int)(nps - MP_edge - tpf_edge);

					if ((0 <= ksyh) && (ksyh < nhfy) && (0 <= ksxh) && (ksxh < nhfx))
					{
						hfx = d_hfx[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];
						hfy = d_hfy[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];
						hfz = d_hfz[(ZCount * nhfy * ksxh) + (ZCount * ksyh) + threadIdx.x];

						// �����_�������̕�ԁB
						// �M�̏��Ɠ����悤�ɂȂ��Ă����܂������Ȃ������̂Ŏ������킩������������B
						// Y���v���X���������邩�m�F����Ԃ̑ΏۂƂ��Ă���B
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

				// ���x�ɂ��v�Z�̐U�蕪���B

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

				// ���M���L�����[���x�𒴂���
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
				// ��ԏ�̑w���w�艷�x��荂����
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
				// ���x���v�Z�͈͊O
//				else
//					;

				// Warp���Œl���܂Ƃ߂�
				for (i = (blockDim.x >> 1); i > 0; i >>= 1)
				{
					__syncthreads();

					if (i > threadIdx.x)
						s_amk[threadIdx.x] += s_amk[threadIdx.x + i];
				}
				__syncthreads();

				// �S�w���v�Z���ʂ𔽉f�����Ȃ��Ȃ炱�̈ʒu�̓X�L�b�v����
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

				// ���Ԗ��̌J��Ԃ�
#ifdef GPU_LOG
				if (0 != s_amk[0])
#endif
				for (t = 0.0f; t < tf; t += dtf)
				{
					__syncthreads();

					to = t + dtf;

// -----------------------------------------------------------------------
//  slsode�𒼏���
// -----------------------------------------------------------------------
					nCalculatingSLSODE = 1;

					// -----------------------------------------------------------------------
					// The following internal Common block contains variables which are
					// communicated between subroutines.  All real variables are listed
					// first, followed by all integers.  The block is declared in
					// Subroutines SLSODE, SINTDY, SSTODE, SPREPJ, and SSOLSY.
					// -----------------------------------------------------------------------
					// �ȉ��̓����R�����E�u���b�N�ɂ́A�T�u���[�`���ԂŒʐM�����ϐ����܂܂�Ă���B�T�u���[�`���ԂŒʐM����܂��B
					// ���ׂĂ̎����ϐ��� ���ɂ��ׂĂ̐����������B
					// ���̃u���b�N�� �T�u���[�`�� SLSODE, SINTDY, SSTODE, SPREPJ, SSOLSY.
					// -----------------------------------------------------------------------
	 
					// -----------------------------------------------------------------------
					// Block A.
					// This code block is executed on every call.
					// It tests ISTATE and ITASK for legality and branches appropriately.
					// If ISTATE .GT. 1 but the flag INIT shows that initialization has
					// not yet been done, an error return occurs.
					// If ISTATE = 1 and TOUT = T, return immediately.
					// -----------------------------------------------------------------------
					// �u���b�NA�B
					// ���̃R�[�h�u���b�N�́A�Ăяo�����ƂɎ��s�����B
					// ISTATE��ITASK�̐��������e�X�g���A�K�؂ɕ��򂷂�B
					// ����ISTATE .GT. 1�ł��邪�A�t���OINIT�����������܂��s���Ă��Ȃ����Ƃ������ꍇ 1�ł��邪�A�t���OINIT�����������܂��s���Ă��Ȃ����Ƃ������Ă���ꍇ�A�G���[���^�[������������B
					// ISTATE = 1����TOUT = T�̏ꍇ�A�����Ƀ��^�[������B
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
						// �u���b�NB�B
						// ���̃R�[�h�u���b�N�́A�ŏ��̌Ăяo��(ISTATE = 1)�ɑ΂��Ď��s�����A
						// �܂��̓p�����[�^�̕ύX�𔺂��p���Ăяo��(ISTATE = 3)�ɑ΂��Ď��s�����B
						// ����ɂ́A���ׂĂ̓��͂̃`�F�b�N�Ƃ��܂��܂ȏ��������܂܂��B
						// 
						// �ŏ��ɁA��I�v�V��������NEQ�AITOL�AIOPT�̍��@�����`�F�b�N����A MF�AML�A�����MU�B
						// -----------------------------------------------------------------------

						// Next process and check the optional inputs. --------------------------

						H0 = 0.0f;

						// -----------------------------------------------------------------------
						// Set work array pointers and check lengths LRW and LIW.
						// Pointers to segments of RWORK and IWORK are named by prefixing L to
						// the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
						// Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
						// -----------------------------------------------------------------------
						// ���[�N�z��̃|�C���^��ݒ肵�A���� LRW �� LIW ���`�F�b�N����B
						// RWORK��IWORK�̃Z�O�����g�ւ̃|�C���^�́AL��擪�ɂ��Ė��������B��t���Ė�������B �Ⴆ�΁A�Z�O�����gYH��RWORK(LYH)����n�܂�B
						// RWORK�̃Z�O�����g�́i���ɁjYH�AWM�AEWT�ASAVF�AACOR�ƕ\�L�����B
						// -----------------------------------------------------------------------

						// -----------------------------------------------------------------------
						// Block C.
						// The next block is for the initial call only (ISTATE = 1).
						// It contains all remaining initializations, the initial call to F,
						// and the calculation of the initial step size.
						// The error weights in EWT are inverted after being loaded.
						// -----------------------------------------------------------------------
						// C�u���b�N�B
						// ���̃u���b�N�́A�ŏ��̌Ăяo���݂̂ł���(ISTATE = 1)�B
						// ���̃u���b�N�ɂ́A�c��̂��ׂĂ̏������AF�A ����я����X�e�b�v�E�T�C�Y�̌v�Z���܂܂��B
						// EWT�̌덷�d�݂́A���[�h��ɔ��]�����B
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
						// �ȉ��̃R�[�f�B���O�́A�ŏ��̃X�e�b�v�Ŏ��s�����X�e�b�v�T�C�YH0���v�Z����B���v�Z����B
						// �ŏ��ɁATOUT - T ���[���Ƒ傫���قȂ邱�Ƃ��`�F�b�N����B
						// �X�J���[���e��TOL���AMAX(RTOL(I)) �����ł���� MAX(RTOL(I))�Ƃ��āC�����łȂ���� MAX(ATOL(I)/ABS(Y(I)))�Ƃ��Čv�Z�����B100*UROUND��1.0E-3�̊ԂɂȂ�悤�ɒ��������B
						// �v�Z�lH0�͎����ŗ^������B
						// 										NEQ
						//   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2 )
						// 										 1
						// �����ŁA	w0		= MAX ( ABS(T), ABS(TOUT) )�A
						// 			f(i)	= f�̏����l��i�Ԗڂ̐����A
						// 			ywt(i)	= EWT(i)/TOL (y(i)�̏d��)�B
						// H0�̕����́ATOUT��T�̏����l���琄�������B
						// -----------------------------------------------------------------------

						if (H0 == 0.0f)
						{
							W0 = max(abs(t), abs(to));

							if (dtf < (2.0f * rumach * W0))
							{
// �v�Z���s���̊m�F�p
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
						// 2��ڈȍ~
						// -----------------------------------------------------------------------
						// Block D.
						// The next code block is for continuation calls only (ISTATE = 2 or 3)
						// and is to check stop conditions before taking a step.
						// -----------------------------------------------------------------------
						// D�u���b�N
						// ���̃R�[�h�u���b�N�́A�p���R�[���iISTATE = 2�܂���3�j��p�ł���B�ŁA�X�e�b�v�����s����O�ɒ�~�������`�F�b�N����B
						// -----------------------------------------------------------------------

						NST = 0;

						if (((TN - to) * H) >= 0.0f)
						{
							if (0 != sintdy(ZCount, L, NQ, to, TN, HU, rumach, H, Yx, Yy, Yz, d_YHx, d_YHy, d_YHz, nYHIndex))
							{
// �v�Z���s���̊m�F�p
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
						// �u���b�NE�B
						// ���̃u���b�N�́A�ʏ킷�ׂĂ̌Ăяo���ɑ΂��Ď��s����A���̓��e���܂ށB�����X�e�b�v�E�R�A�ϕ���SSTODE�̌Ăяo�����܂܂��B
						// ����͓����X�e�b�v�̃��[�v�|�C���g�ł���B
						// �ŏ��ɃX�e�b�v���������Ȃ����`�F�b�N���AEWT���X�V����B���X�V���A�v������鐸�x���������Ȃ����`�F�b�N����BT�̃��E���h�I�t�E���x���ȉ���H���`�F�b�N����B
						// -----------------------------------------------------------------------

// 250
						if (NST >= MXSTP0)
						{
// �v�Z���s���̊m�F�p
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

// �v�Z���s���̊m�F�p
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

// �v�Z���s���̊m�F�p
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
//  sstode�𒼏���
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
							// �ŏ��̌Ăяo���ŁA�I�[�_�[��1�ɃZ�b�g����A���̕ϐ��������������B�������������B
							// RMAX�́A1�X�e�b�v��H�𑝉������邱�Ƃ��ł���ő�䗦�ł���B�ł���B
							// �����l��1.E4�ł���B�ł��邪�A�ʏ��10�ɓ������B
							// ���� �����������ꍇ(�R���N�^�̎����܂��̓G���[�e�X�g)�A RMAX �� 2 �ɐݒ肳���B�ɐݒ肳���B
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
							// H���ύX�����ꍇ�AH��RH��RMAX�AHMIN�AHMXI�Əƍ�����AYH�z�񂪍ăX�P�[�����O�����BRMAX�CHMIN�CHMXI�Əƍ�����CYH�z�񂪍ăX�P�[�����O�����B
							// IALTH�� �ɐݒ肳���B�ɐݒ肳���B
							// -----------------------------------------------------------------------

					// 140

							nCalculatingSSTODE = 150;
						}

						// 700�ɂȂ�Ɣ�������
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
								// H ���ύX�����ꍇ�AH �� RH �� RMAX�AHMIN�AHMXI �Əƍ�����AYH �z�񂪍ăX�P�[�����O�����B
								// IALTH �� L = NQ + 1 �ɐݒ肳��A�����܂��̓G���[�e�X�g�̎��s�ɂ���ċ�������Ȃ�����A���̃X�e�b�v������ H �̕ύX��h�~����B
								// -----------------------------------------------------------------------

					// 170
								RH = max(RH, 0.0f);

					// 175
								RH = min(RH, RMAX);
								//RH = RH / max(1.0f, (abs(H) * 0.0f * RH));	1�Ŋ���̂Ŗ���
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
							// ���̃Z�N�V�����ł́AYH �z��Ƀp�X�J���E�g���C�A���O���s������ʓI�Ɋ|�����킹�邱�Ƃɂ��A�\���l���v�Z����BYH �z��� Pascal Triangle �s�����Z����B
							// RC �͌W�� H*EL(1) �̐V�����l�ƌÂ��l�̔�ł���B
							// RC �� 1 �� CCMAX �ȏ�قȂ�ꍇ�AIPUP �� MITER �ɐݒ肳��APJAC �������I�ɌĂяo����܂��B�ɐݒ肳��A���R�r�A�����֗^����ꍇ��PJAC�������I�ɌĂяo�����B
							// ������ɂ���APJAC �͏��Ȃ��Ƃ� MSBP �X�e�b�v���ƂɌĂяo�����B
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
							// �ő� MAXCOR �̕␳�̔������s����B
							// �����e�X�g�� �e�␳��R.M.S.�m�����ɂ��Ď����e�X�g���s����B�d�݃x�N�g�� EWT �ɂ���ďd�ݕt������܂��B
							// �␳�ʂ̍��v�� �x�N�g��ACOR(i)�ɗݐς����B
							// YH �z��͕␳���[�v�̒��ł͕ύX����Ȃ��B
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
							// �w�����ꂽ�ꍇ�AP = I - h*el(1)*J�s�񂪍ĕ]������A�␳��̔������J�n����O�ɑO�������s����B ���ĕ]������, �␳��̔������J�n����O�ɑO���������.
							// IPUP �� �� 0 ���Z�b�g�����B
							// -----------------------------------------------------------------------

							ACORx = 0.0;
							ACORy = 0.0;
							ACORz = 0.0;

							// ���ɐi�߂�Ȃ�0�ł͂Ȃ��Ȃ�
							while (0 == nCalculatingSSTODE)
							{
								__syncthreads();

								//-----------------------------------------------------------------------
								// In the case of functional iteration, update Y directly from
								// the result of the last function evaluation.
								// -----------------------------------------------------------------------
								// �֐������̏ꍇ�A�Ō�̊֐��]������C���璼��Y���X�V����B�𒼐ڍX�V���܂��B
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
								// �����̃e�X�g�B
								// M.gt.0�̏ꍇ�A�������萔�̐���l��CRATE�Ɋi�[����A�e�X�g�Ɏg�p�����B�̐���l�� CRATE �Ɋi�[����A���ꂪ�e�X�g�Ɏg�p�����B
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
								// �␳��͎��������B
								// JCUR �� 0 �ɃZ�b�g�����B
								// ���[�J���E�G���[�E�e�X�g���s���A���s�����ꍇ�̓X�e�[�g�����g500 �ɓn�����B
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
									// �X�e�b�v������AYH�z����X�V����B
									// IALTH = 1�ł����H�̕ύX����������B
									// �����łȂ����IALTH��1���炷�B
									// IALTH��1�ł���ANQ .lt. MAXORD�̏ꍇ�AACOR�� ���ۑ������B
									// H�̕ύX���l�������ꍇ�A�I�[�_�[C��1�������l�������B�̑������l�������B
									// H�̕ύX�� �{�ȏ�ł���ꍇ�ɂ̂ݍs����B
									// �����łȂ��ꍇ�AIALTH �� 3 �ɐݒ肳��A���̐��X�e�b�v�� C �e�X�g���s���Ȃ��悤�ɂ���B�ɐݒ肳���B
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
									// �G���[�e�X�g�����s�����B KFLAG�͕����̎��s���L�^����B
									// TN��YH�A���C���ȑO�̒l�ɖ߂��A�X�e�b�v���Ď��s���鏀��������B�X�e�b�v���Ď��s���鏀��������B
									// �œK�ȃX�e�b�v�E�T�C�Y���v�Z����B���̎����ɂ��čœK�ȃX�e�b�v�T�C�Y���v�Z����B
									// 2��ȏ�̎��s�̌�AH�͋����I�� ��0.2�{�ȉ��ɂ���B
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
										// ���s��3��ȏ㔭�������ꍇ�A����͂��̃Z�N�V�����ɓ��B����B
										// ���s�� 10 �񔭐������ꍇ�́AKFLAG = -1 �ŏI������B
										// YH�z��ɒ~�ς��ꂽ���֐��́A�Ԉ���������̃G���[�������Ă���Ɖ��肳���B
										// ���������āA�ŏ��̓��֐����Čv�Z����A������1�ɐݒ肳���B
										// ������ H��10�{���A�X�e�b�v���Ď��s����A �������邩H��HMIN�ɒB����܂ŁB
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
								// �X�e�b�v�̐��ۂɂ�����炸�A�W�� RHDN�ARHSM�ARHUP���v�Z����A���ꂼ��H������NQ - 1�܂��͎���NQ + 1�� ���悶�邱�Ƃ��ł���B
								// ���s�̏ꍇ�A�����̑���������邽�߂�RHUP=0.0�Ƃ���B
								// �����̂����ő�̂��̂����肳��A����ɉ����ĐV�����I�[�_�[���I�������B��I������B
								// �����𑝉�������ꍇ�́A���̂悤�Ɍv�Z����B���v�Z����B
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
									// �I�[�_�[���ύX���ꂽ�ꍇ�́ANQ�Al�A�W�������Z�b�g����B
									// ������̏ꍇ���AH �� RH �ɏ]���ă��Z�b�g����AYH �z��͍ăX�P�[�����O�����B
									// �X�e�b�v�ɖ�肪�Ȃ����690���I�����A�����łȂ���΃X�e�b�v����蒼���B
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
							// �␳��̔������������Ȃ������B
							// ���� MITER .ne. 0 �Ń��R�r�A�����Â��ꍇ�APJAC ���Ă΂��B���Ăяo�����B
							// �����łȂ���΁AYH �z��͗\���O�̒l�ɖ߂����B�ɖ߂���A�\�ł���� H �����炳���B
							// H���k���ł��Ȃ��ꍇ H���k���ł��Ȃ����AMXNCF�Ɏ��s�����ꍇ�́AKFLAG = -2�ŏI������B
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
//  sstode�̏I���
// -----------------------------------------------------------------------

						// sstode�ł̃G���[�`�F�b�N
						if (1 != (1 - KFLAG))
						{
// �v�Z���s���̊m�F�p
d_ymx[nIndex] = (float)KFLAG;
d_ymy[nIndex] = H;
d_ymz[nIndex] = 66666.0f;
							return;
//							break;
						}

						// �����ُ̈픻��
						if (0.0f == H)
						{
// �v�Z���s���̊m�F�p
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
						// F�u���b�N
						// ���̃u���b�N�́A�R�A�C���e�O���[�^����̃��^�[�������iKFLAG = 0�j�̏ꍇ����������B(KFLAG=0)�B ��~�������e�X�g����B
						// -----------------------------------------------------------------------
// 300
						// ITASK = 1.  If TOUT has been reached, interpolate. -------------------
						if (((TN - to) * H) < 0.0)
						{
							continue;
						}

						if (0 != sintdy(ZCount, L, NQ, to, TN, HU, rumach, H, Yx, Yy, Yz, d_YHx, d_YHy, d_YHz, nYHIndex))
						{
// �v�Z���s���̊m�F�p
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
					// G�u���b�N
					// ���̃u���b�N�́ASLSODE����̂��ׂĂ̐������^�[������������B
					// ����ITASK .NE. 1 �̏ꍇ�AYH ���� Y �����[�h����A����ɉ����� T ���ݒ肳���B
					// ISTATE �� 2 �ɐݒ肳��A�I�v�V�����o�͂����[�N�E�A���C�Ƀ��[�h�����B
					// -----------------------------------------------------------------------
// 420
					__syncthreads();
					if (0 == threadIdx.x)
						n_Temp[SHARE_ISTATE] = 2;

// -----------------------------------------------------------------------
//  slsode�̏I���
// -----------------------------------------------------------------------

					yamp = sqrt((Yx * Yx) + (Yy * Yy) + (Yz * Yz));

					// 1�ł����������w������ΑS�̂��͂��߂����蒼���܂��B (���������Ԃ͊����߂��Ȃ�)
					if (0.02f < abs(yamp - 1.0f))
						n_Temp[SHARE_ISTATE] = 1;

					Yx /= yamp;
					Yy /= yamp;
					Yz /= yamp;
				}

				// �A�j���[�V�����p���O
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


