#include "CHAMR.h"

#include "CUDAInterface.cuh"

#include <iostream>
#include <algorithm>
#include <filesystem>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

///////////////////////////////////////////////////////////////////////////////////////////////
// 
CHAMR::CHAMR()
{
	m_ZCount = 0;
	m_Grains = 0;
	m_nqs3 = 0;

	m_ntpx = 0;
	m_ntpy = 0;
	m_ntpz = 0;
	m_nhfx = 0;
	m_nhfy = 0;
	m_nhfz = 0;
	m_MP_edge = 0.0f;
	m_tpf_edge = 0.0f;
	m_gymax = 0.0;
	m_kxofs = 0;

	m_ntmpw = 0;
	m_nwr = 0;

	m_dtc = 0.0;
	m_tambient = 0.0;
	m_tpk = 0.0;
	m_tgk = 0.0;
	m_spd = 0.0;
	m_pmove = 0.0;
	m_pstart = 0.0;
	m_offset = 0.0;
	m_gsz = 0.0;
	m_asz = 0.0;
	m_dtsimi = 0.0;
	m_antf = 0.0;
	m_hapi = 0.0;
	m_trise = 0.0;
	m_nps = DUMMY_DOUBLE;
	m_sbit = 0.0;
	m_sdelay = 0.0;
	m_npwrite = 0;
	m_nmwrite = 0;
	m_ntwrite = 0;
	m_animation = 0;
	m_ichki = 0;
	m_sigma_gs = 0.0;
	m_gstc_d = 0.0;
	m_gstc_v = 0.0;

	memset(m_fpfx, 0, 256);
	memset(m_fpa_in, 0, 256);

	memset(m_fimed, 0, 256);
	memset(m_fihfd, 0, 256);
	memset(m_fitpf, 0, 256);
	memset(m_fitpf_re, 0, 256);
	memset(m_fiend, 0, 256);

	m_nLLG = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
CHAMR::~CHAMR()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
bool CHAMR::Run(const char* pszPath, const char* pszPattern, const char* pszParameter, const int nLogGrain)
{
	time_t timeStart = time(NULL);

	// �����̎�B
	char szBuffer[FILE_READ_BUFFER] = "";
	strcpy(szBuffer, pszPattern);

	// �v�Z�����������g�����Čv�Z�̏ꍇ�A�p�^�[������10�𒴂��鎞������B
	int nPatternIndex = 1;
	if (3 < strlen(szBuffer))
		nPatternIndex = (int)strlen(szBuffer) - 2;

	const int nrm = strtol(&szBuffer[nPatternIndex], NULL, 10);

	// for�񂵂ōs�I�����Ă���̂�0����for�����Ȃ��B
	szBuffer[nPatternIndex] = '\0';
	const int jgp = strtol(szBuffer, NULL, 10) + 1;

	// ���ʂ̃t�@�C�����ɂȂ�B
	strcpy(m_fpfx, pszPattern);

	// �����t�@�C���B
	sprintf(m_fpa_in, "%s/%s", pszPath, pszParameter);

	printf("Program Starts!\n");
	printf("___    \n");
	printf("Chosen Grain Pattern %d\n", jgp);
	printf("___    \n");

	printf("\n");
	printf("Reading the parameter file\n");
	printf("\n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// �����t�@�C���ǂݍ���

	FILE *pFilePAR = fopen(m_fpa_in, "rt");
	if (NULL == pFilePAR)
	{
		printf("Failed to file open. %s\n", m_fpa_in);
		return false;
	}

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f", &m_vhksigma, &m_alftscaling);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f %f %f", &m_dtc, &m_tambient, &m_tpk, &m_tgk);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f %f", &m_spd, &m_pmove, &m_pstart);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f", &m_offset);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f %d", &m_gsz, &m_asz, &m_ZCount);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f", &m_dtsimi, &m_antf);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f %f", &m_hapi, &m_trise, &m_nps);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f", &m_sbit, &m_sdelay);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%d %d %d %d", &m_npwrite, &m_nmwrite, &m_ntwrite, &m_animation);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%d", &m_ichki);

	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%f %f %f", &m_sigma_gs, &m_gstc_d, &m_gstc_v);

	//'.med'
	char szFile[256] = "";
	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%s", szFile);
	sprintf(m_fimed, "%s/%s%s", pszPath, szFile, EXTENSION_MED);

	//'.tpf'
	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%s", szFile);
	sprintf(m_fitpf, "%s/%s%s", pszPath, szFile, EXTENSION_TPF);
	sprintf(m_fitpf_re, "%s/%s_re%s", pszPath, szFile, EXTENSION_TPF);

	//'.hfd'
	fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
	printf("%s", szBuffer);
	sscanf(szBuffer, "%s", szFile);
	sprintf(m_fihfd, "%s/%s%s", pszPath, szFile, EXTENSION_HFD);

	// �g�p����p�^�[���t�@�C���܂ōs���ړ����܂��B
	for (int i = 0; i < jgp; i++)
	{
		fgets(szBuffer, FILE_READ_BUFFER, pFilePAR);
		printf("%s", szBuffer);
	}

	char pat[256] = "";
	sscanf(szBuffer, "%s", pat);

	char pat_ga[256] = "";
	sprintf(pat_ga, "%s/%s%s", pszPath, pat, EXTENSION_GA_DAT);	// '_ga.dat'
	printf("%s\n", pat_ga);

	char pat_gp[256] = "";
	sprintf(pat_gp, "%s/%s%s", pszPath, pat, EXTENSION_GP_DAT);	// '_gp.dat'
	printf("%s\n", pat_gp);

	fclose(pFilePAR);
	pFilePAR = nullptr;

	printf("\n");
	printf("Finished Reading Parameter File\n");
	printf("\n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// .med�t�@�C���ǂݍ���

	// Z�̌������ςɂł���p���med�t�@�C����ǂݍ��ށB
	// med�t�@�C����Z�̌����ɑΉ����Ă��Ȃ����͏]����parameter�t�@�C��������������B
	m_tempc = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_tempc)
		return false;

	m_sigmatc = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_sigmatc)
		return false;

	m_ahko = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_ahko)
		return false;

	m_sigmahk = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_sigmahk)
		return false;

	m_amsrt = std::make_unique<int[]>(m_ZCount);
	if (nullptr == m_amsrt)
		return false;

	m_aexo = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_aexo)
		return false;

	m_ept = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_ept)
		return false;

	m_eptAex = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_eptAex)
		return false;

	m_alf = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_alf)
		return false;

	m_alfg = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_alfg)
		return false;

	m_alft = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_alft)
		return false;

	m_amso = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_amso)
		return false;

	m_tcbulk = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_tcbulk)
		return false;

	m_sigmatc_gs = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_sigmatc_gs)
		return false;

	m_sigmatc_ots = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_sigmatc_ots)
		return false;

	m_Thickness = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_Thickness)
		return false;

	for (int z = 0; z < m_ZCount; z++)
		m_Thickness[z] = m_asz;

	FILE *pFileFIMED = fopen(m_fimed, "rt");

	for (int z = 0; z < m_ZCount; z++)
	{
		int iz = 0;
		float thickness = 0.0f;

		fgets(szBuffer, FILE_READ_BUFFER, pFileFIMED);
		sscanf(szBuffer, "%d %f %f %f %f %d %f %f %f %f %f", &iz, &m_tempc[z], &m_sigmatc[z], &m_ahko[z], &m_sigmahk[z], &m_amsrt[z], &m_aexo[z], &m_ept[z], &m_eptAex[z], &m_alf[z], &thickness);

		m_tcbulk[z] = m_tempc[z] / (1.0f - (m_gstc_d * pow(m_gsz, (-1.0f / m_gstc_v))));
		float partial_derivative = -(m_tempc[z] * m_gstc_d) / (m_gstc_v * pow(m_gsz, (1.0f + (1.0f / m_gstc_v))));
		float stddev_gs = m_gsz * m_sigma_gs;
		m_sigmatc_gs[z] = abs(partial_derivative) * stddev_gs / m_tempc[z];
		m_sigmatc_ots[z] = sqrt(pow(m_sigmatc[z], 2.0f) - pow(m_sigmatc_gs[z], 2.0f));

		m_alft[z] = m_alfg[z] = m_alf[z];
		m_amso[z] = (float)m_amsrt[z] / pow((1.0f - 300.0f / m_tempc[z]), (1.0f / 3.0f));

		if (0.0f != thickness)
			m_Thickness[z] = thickness;

		if (z != (iz - 1))
		{
			printf("Faulty reading in medium file\n");
			return false;
		}
	}

	fclose(pFileFIMED);
	pFileFIMED = nullptr;

	///////////////////////////////////////////////////////////////////////////////////////////////
	// .tpf�t�@�C���ǂݍ���

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Read-in thermal profile
	///////////////////////////////////////////////////////////////////////////////////////////////

	FILE *pFileTPF = fopen(m_fitpf, "rt");

	int ndumb = 0;
	fgets(szBuffer, FILE_READ_BUFFER, pFileTPF);
	sscanf(szBuffer, "%d %d %d %d %f", &m_ntpx, &m_ntpy, &m_ntpz, &ndumb, &m_tpf_edge);

	m_atpx = std::make_unique<float[]>(m_ntpx * m_ntpy * m_ntpz);
	if (nullptr == m_atpx)
		return false;

	m_atpy = std::make_unique<float[]>(m_ntpx * m_ntpy * m_ntpz);
	if (nullptr == m_atpy)
		return false;

	m_atpz = std::make_unique<float[]>(m_ntpx * m_ntpy * m_ntpz);
	if (nullptr == m_atpz)
		return false;

	m_atpt = std::make_unique<float[]>(m_ntpx * m_ntpy * m_ntpz);
	if (nullptr == m_atpt)
		return false;

	int indx = 0;
	float fMaxTpt = 0.0f;
	for (int z = 0; z < m_ntpz; z++)
	{
		for (int y = 0; y < m_ntpy; y++)
		{
			for (int x = 0; x < m_ntpx; x++)
			{
				indx = (m_ntpz * m_ntpy * x) + (m_ntpz * y) + z;
				fgets(szBuffer, FILE_READ_BUFFER, pFileTPF);
				sscanf(szBuffer, "%f %f %f %f", &m_atpx[indx], &m_atpy[indx], &m_atpz[indx], &m_atpt[indx]);

				// �ő剷�x�����߂�B
				if (fMaxTpt < m_atpt[indx])
					fMaxTpt = m_atpt[indx];
			}
		}
	}

	fclose(pFileTPF);
	pFileTPF = nullptr;

	// Tpf�Čv�Z
	// �p�����[�^�w�肪�������ꍇ�̂�
	if ((0.0f < m_tpk) && (0.0f < m_tgk))
	{
		printf("Recalculate tpf\n");
		printf("tpk=%.2f tgk=%.2fn\n", m_tpk, m_tgk);
		printf(" \n");

		const float fZeroTemp = 273.15f;
		const float fSubMaxTpt = fMaxTpt - fZeroTemp;

		// (T - 273.15[K]) / (Tmax - 273.15[K]) * (Tpk - Tgk) + Tgk

		pFileTPF = fopen(m_fitpf_re, "wt");
		if (nullptr != pFileTPF)
		{
			fprintf(pFileTPF, "%d %d %d %g\n", m_ntpx, m_ntpy, m_ntpz, m_tpf_edge);
			for (int z = 0; z < m_ntpz; z++)
			{
				for (int y = 0; y < m_ntpy; y++)
				{
					for (int x = 0; x < m_ntpx; x++)
					{
						indx = (m_ntpz * m_ntpy * x) + (m_ntpz * y) + z;
						m_atpt[indx] = (m_atpt[indx] - fZeroTemp) / fSubMaxTpt * (m_tpk - m_tgk) + m_tgk;

						fprintf(pFileTPF, "%g %g %g %g\n", m_atpx[indx], m_atpy[indx], m_atpz[indx], m_atpt[indx]);
					}
				}
			}
		}

		fclose(pFileTPF);
		pFileTPF = nullptr;
	}

	printf("unexpected ending\n");

	printf("Thermal Profile\n");
	printf("%d %d %d\n", m_ntpx, m_ntpy, m_ntpz);
	printf("Finished Reading Thermal File\n");
	printf(" \n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Z���̔M�쐬

	// Head-Medium Spacing
	float hms = (float)m_atpz[0];

	// Top Co Layer
	for (int z = 0; z < m_ZCount; z++)
		m_Thickness[z] = hms + (m_Thickness[z] * (float)z) + 0.0001f;

	m_tpf = std::make_unique<float[]>(m_ntpx * m_ntpy * m_ZCount);
	if (nullptr == m_tpf)
		return false;

	for (int z = 0; z < m_ZCount; z++)
	{
		int kzs = 0, kzt = 1;

		while (m_Thickness[z] > m_atpz[kzt])
		{
			if (kzt >= (m_ntpz - 1))
				break;

			kzs++;
			kzt++;
		}

		if (kzt < m_ntpz)
		{
			for (int i = 0; i < m_ntpx; i++)
			{
				const int nIndex = m_ntpy * m_ntpz * i;
				for (int j = 0; j < m_ntpy; j++)
				{
					const int nIndex1 = nIndex + (m_ntpz * j) + kzs;
					const int nIndex2 = nIndex + (m_ntpz * j) + kzt;
					float tmp1 = m_atpt[nIndex1];
					float tmp2 = m_atpt[nIndex2];
					float atz1 = m_atpz[nIndex1];
					float atz2 = m_atpz[nIndex2];
					m_tpf[(m_ZCount * m_ntpy * i) + (m_ZCount * j) + z] = tmp1 + (tmp2 - tmp1) * (m_Thickness[z] - atz1) / (atz2 - atz1);
				}
			}
		}
		else
		{
			for (int i = 0; i < m_ntpx; i++)
			{
				for (int j = 0; j < m_ntpy; j++)
				{
					m_tpf[(m_ZCount * m_ntpy * i) + (m_ZCount * j) + z] = m_atpt[(m_ntpz * m_ntpy * i) + (m_ntpz * j) + kzt];
				}
			}
		}
	}

	printf("Head Field Factor  %g\n", m_hapi);

	float gvm = m_gsz * m_gsz * m_asz * 1.0e-21f;
	float gvmn = m_gsz * m_gsz * m_asz * 1.0e-21f;

	float tpc_avg = 700.0f;
	float amsi = 1000.0f;
	float akb = 1.38e-16f;
	float akbtc = akb * tpc_avg;
	float hnm = akbtc / (gvmn * amsi);

	printf(" \n");
	printf(" Normalization field hnm %g\n", hnm);
	printf(" \n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// .hfd�t�@�C���ǂݍ���

	// Input head field
	FILE *pFileFIHFD = fopen(m_fihfd, "rt");

	int nzhfd = 0, hmshfd = 0, xphfd = 0;

	fgets(szBuffer, FILE_READ_BUFFER, pFileFIHFD);
	sscanf(szBuffer, "%d %d %d %f", &nzhfd, &hmshfd, &xphfd, &m_MP_edge);

	// 3D�̃w�b�h���g�����̐U�蕪���B
	// �t�H�[�}�b�g���Ⴂ���`�F�b�N���Ă���B
	float fDummy1 = 0.0, fDummy2 = DUMMY_DOUBLE;
	fgets(szBuffer, FILE_READ_BUFFER, pFileFIHFD);
	sscanf(szBuffer, "%f %f %f %f %f %f", &fDummy1, &fDummy1, &fDummy1, &fDummy1, &fDummy1, &fDummy2);

	// ���̓ǂݍ��݈ʒu�ɖ߂��B
	fseek(pFileFIHFD, SEEK_SET, 0);
	fgets(szBuffer, FILE_READ_BUFFER, pFileFIHFD);

	// �̂���̃w�b�h
	printf("check 3DHead\n");
	printf("nps=%f\n", m_nps);
	printf("tpf_edge=%f\n", m_tpf_edge);
	printf("MP_edge=%f\n", m_MP_edge);
	printf("r2c6=%f\n\n", fDummy2);
	if ((DUMMY_DOUBLE == m_nps) || (DUMMY_DOUBLE == fDummy2))
	{
		m_nhfx = 1;
		m_nhfy = 1;
		m_nhfz = m_ZCount;

		m_hfix = std::make_unique<float[]>(nzhfd);
		if (nullptr == m_hfix)
			return false;

		m_hfiy = std::make_unique<float[]>(nzhfd);
		if (nullptr == m_hfiy)
			return false;

		m_hfiz = std::make_unique<float[]>(nzhfd);
		if (nullptr == m_hfiz)
			return false;

		m_zhfd = std::make_unique<float[]>(nzhfd);
		if (nullptr == m_zhfd)
			return false;

		for (int i = 0; i < nzhfd; i++)
		{
			fgets(szBuffer, FILE_READ_BUFFER, pFileFIHFD);
			sscanf(szBuffer, "%f %f %f %f", &m_hfix[i], &m_hfiy[i], &m_hfiz[i], &m_zhfd[i]);
		}

		fclose(pFileFIHFD);
		pFileFIHFD = nullptr;

		if (m_Thickness[m_ZCount - 1] > m_zhfd[nzhfd - 1])
		{
			printf("Head Field Short of Depth!, Abort!\n");
			return false;
		}

		if (m_Thickness[0] < m_zhfd[0])
		{
			printf("Head Field HMS Too Large!, Abort!\n");
			return false;
		}

		m_hfx = std::make_unique<float[]>(m_ZCount);
		if (nullptr == m_hfx)
			return false;

		m_hfy = std::make_unique<float[]>(m_ZCount);
		if (nullptr == m_hfy)
			return false;

		m_hfz = std::make_unique<float[]>(m_ZCount);
		if (nullptr == m_hfz)
			return false;

		// Head Field Interpolation
		for (int z = 0; z < m_ZCount; z++)
		{
			int iz = 0;
			int izs = 0;
			float zpp = m_Thickness[z];
			float zhs = m_zhfd[iz];

			while (m_zhfd[iz] <= zpp)
			{
				zhs = m_zhfd[iz];
				izs = iz;
				iz++;
			}

			float zht = m_zhfd[iz];

			m_hfy[z] = m_hfiy[izs] + (m_hfiy[iz] - m_hfiy[izs]) * (1.0f - (zht - zpp) / (zht - zhs));
			m_hfz[z] = m_hfiz[izs] + (m_hfiz[iz] - m_hfiz[izs]) * (1.0f - (zht - zpp) / (zht - zhs));

			m_hfx[z] = 0.0;
			m_hfy[z] = m_hapi * m_hfy[z] / hnm;
			m_hfz[z] = m_hapi * m_hfz[z] / hnm;
		}
	}
	// 3D�w�b�h
	else
	{
		// 3D�w�b�h�̓w�b�_�[��x,y,z,MP_edge�ƂȂ��Ă���B
		// �̂̂�z,hms,?�Ȃ̂Œ��ӂ���B
		// ���ӂ��Ďg���Ζ��Ȃ��B
		// Y���̃w�b�h���E�͌㔼�����������Ȃ����߁A�Ō�̏�Ԃ����̂܂܉������Ă��܂��B
		// ���W�����̂܂܎g���ɂ͑O��ɉ�������K�v������܂��B
		const int nhfy = hmshfd;
		const int nhfy_offset = nhfy / 10;
		m_nhfx = nzhfd;
		m_nhfy = nhfy + (nhfy_offset * 2);
		m_nhfz = xphfd;

		printf("3DHead used. x=%d y=%d z=%d\n\n", m_nhfx, nhfy, m_nhfz);

		// ���������W�̊�łȂ��S�̂̒��S�����(0,0)�Ȃ��߁AY�ʒu�𒲐����܂��B
//		m_nps += (float)(m_nhfy / 2);
//		m_tpf_edge -= (float)(m_ntpy / 2);

		const int nSizeFileLine = m_nhfx * nhfy * m_nhfz;

		m_hfix = std::make_unique<float[]>(nSizeFileLine);
		if (nullptr == m_hfix)
			return false;

		m_hfiy = std::make_unique<float[]>(nSizeFileLine);
		if (nullptr == m_hfiy)
			return false;

		m_hfiz = std::make_unique<float[]>(nSizeFileLine);
		if (nullptr == m_hfiz)
			return false;

		m_hftx = std::make_unique<float[]>(nSizeFileLine);
		if (nullptr == m_hftx)
			return false;

		m_hfty = std::make_unique<float[]>(nSizeFileLine);
		if (nullptr == m_hfty)
			return false;

		m_hftz = std::make_unique<float[]>(nSizeFileLine);
		if (nullptr == m_hftz)
			return false;

		const int nSizeData = m_nhfx * m_nhfy * m_ZCount;
		m_hfx = std::make_unique<float[]>(nSizeData);
		if (nullptr == m_hfx)
			return false;

		m_hfy = std::make_unique<float[]>(nSizeData);
		if (nullptr == m_hfy)
			return false;

		m_hfz = std::make_unique<float[]>(nSizeData);
		if (nullptr == m_hfz)
			return false;

		// 3D�w�b�h��x,y,z,Hx,Hy,Hz�̏��ɓ����Ă���B
		indx = 0;
		for (int i = 0; i < nSizeFileLine; i++)
		{
			fgets(szBuffer, FILE_READ_BUFFER, pFileFIHFD);
			sscanf(szBuffer, "%f %f %f %f %f %f", &m_hfix[indx], &m_hfiy[indx], &m_hfiz[indx], &m_hftx[indx], &m_hfty[indx], &m_hftz[indx]);
			indx++;
		}

		// tpf�Ɠ����悤�ȍ��W�ɕ�Ԃ���B
		// �W���������킹�čs���B
		for (int z = 0; z < m_ZCount; z++)
		{
			int kzs = 0, kzt = 1;

			// �w�̌������ǂ̈ʒu���`�F�b�N�B
			while (m_Thickness[z] > m_hfiz[m_nhfx * nhfy * kzt])
			{
				if (kzt >= (m_nhfz - 1))
					break;

				kzs++;
				kzt++;
			}

			// Z�͈͓̔��Ȃ��ߕ�Ԃł���B
			if (kzt < m_nhfz)
			{
				for (int i = 0; i < m_nhfx; i++)
				{
					for (int j = 0; j < nhfy; j++)
					{
						// x -> y -> z �̏��ɃC���N�������g���ꂽ�f�[�^�̃C���f�b�N�X�쐬�B
						const int nIndex1 = i + (j * m_nhfx) + (kzs * m_nhfx * nhfy);
						const int nIndex2 = i + (j * m_nhfx) + (kzt * m_nhfx * nhfy);
//						const int nIndex1 = (i * m_nhfz * nhfy) + (j * m_nhfz) + kzs;
//						const int nIndex2 = (i * m_nhfz * nhfy) + (j * m_nhfz) + kzt;
						const int nIndex3 = (m_ZCount * m_nhfy * i) + (m_ZCount * (j + nhfy_offset)) + z;

						// ��Ԃ͐��`�B
						const float atz1 = m_hfiz[nIndex1];
						const float atz2 = m_hfiz[nIndex2];
						const float atzMove = m_Thickness[z] - atz1;
						const float atzDelta = atz2 - atz1;

						float tmp1 = m_hftx[nIndex1];
						float tmp2 = m_hftx[nIndex2];
						m_hfx[nIndex3] = m_hapi * (tmp1 + (tmp2 - tmp1) * atzMove / atzDelta) / hnm;

						tmp1 = m_hfty[nIndex1];
						tmp2 = m_hfty[nIndex2];
						m_hfy[nIndex3] = m_hapi * (tmp1 + (tmp2 - tmp1) * atzMove / atzDelta) / hnm;

						tmp1 = m_hftz[nIndex1];
						tmp2 = m_hftz[nIndex2];
						m_hfz[nIndex3] = m_hapi * (tmp1 + (tmp2 - tmp1) * atzMove / atzDelta) / hnm;
					}
				}
			}
			// �ȍ~��Z�������̂ŕ�Ԃł��Ȃ��B
			else
			{
				for (int i = 0; i < m_nhfx; i++)
				{
					for (int j = 0; j < nhfy; j++)
					{
						const int nIndex1 = i + (j * m_nhfx) + (kzs * m_nhfx * nhfy);
						const int nIndex3 = (m_ZCount * m_nhfy * i) + (m_ZCount * (j + nhfy_offset)) + z;
						m_hfx[nIndex3] = m_hapi * m_hftx[nIndex1] / hnm;
						m_hfy[nIndex3] = m_hapi * m_hfty[nIndex1] / hnm;
						m_hfz[nIndex3] = m_hapi * m_hftz[nIndex1] / hnm;
					}
				}
			}
		}

		// ����������㔼�̕������R�s�[�ō쐬���܂��B
		// �O���͉������܂���B
		for (int z = 0; z < m_ZCount; z++)
		{
			for (int i = 0; i < m_nhfx; i++)
			{
				for (int j = 0; j < nhfy_offset; j++)
				{
					const int nIndex3 = (m_ZCount * m_nhfy * i) + (m_ZCount * (nhfy + nhfy_offset + j)) + z;
					const int nIndex4 = (m_ZCount * m_nhfy * i) + (m_ZCount * (nhfy + nhfy_offset - 1)) + z;
					m_hfx[nIndex3] = m_hfx[nIndex4];
					m_hfy[nIndex3] = m_hfy[nIndex4];
					m_hfz[nIndex3] = m_hfz[nIndex4];
				}
			}
		}
	}

	// alf      Thermal bath coupling lambda in Damping vs Tem
	// alfg     Gilbert Damping Constant for the LLG equation
	// alft     Damping Constant Used for Calculating Thermal Field

	m_heb = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_heb)
		return false;

	m_hko = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_hko)
		return false;

	for (int i = 0; i < m_ZCount; i++)
	{
		m_heb[i] = 2.0f * m_aexo[i] / (hnm * m_amso[i] * m_asz * m_asz * 1.0e-14f);
		m_hko[i] = (m_ahko[i] / hnm) / pow((1.0f - 300.0f / m_tempc[i]), (m_ept[i] / 3.0f));
	}

	printf("  \n");
	printf("Exchange Field %g\n", (m_heb[1] * hnm));
	printf("  \n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Input Medium Grains
	///////////////////////////////////////////////////////////////////////////////////////////////


	m_nqs3 = 3 * m_ZCount;
	printf("nqt %d\n", m_nqs3);

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Setup output files
	///////////////////////////////////////////////////////////////////////////////////////////////

	const char* fo_avg = "avg_";		//fpfx
	const char* fo_dym = "dym_";		//fpfx

	printf("start simulation\n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Starting Simulation
	///////////////////////////////////////////////////////////////////////////////////////////////

	float tbit = m_sbit / m_spd;
	float tdelay = m_sdelay / m_spd;

	printf(" \n");
	printf(" dtsim = %g\n", m_dtsimi);
	printf(" ");

	printf(" \n");

	m_htho = std::make_unique<float[]>(m_ZCount);
	if (nullptr == m_htho)
		return false;

	float aaa= 2.0f * akb;
	for (int i = 0; i < m_ZCount; i++)
	{
		float bbb = m_amso[i] * gvm * GAMA * m_dtsimi;
		m_htho[i] = sqrt(aaa / bbb) / hnm;
		printf("Thermal Field %d %g\n", i, (m_htho[i] * hnm));
	}

	printf(" \n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Start read grain pattern files
	///////////////////////////////////////////////////////////////////////////////////////////////

	printf("Grain Pattern File\n");
	printf("%s\n", pat);
	printf(" \n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// open grain center position file
	// gps(i) is the area of ith grain in nm^2

	FILE *pFilePAT = fopen(pat_ga, "rt");
	indx = 0;
	float aindx = 0.0, dummy;
	while (NULL != fgets(szBuffer, FILE_READ_BUFFER, pFilePAT))
	{
		sscanf(szBuffer, "%f %f %f %f", &aindx, &dummy, &dummy, &dummy);
		if (0.0 != aindx)
			indx++;
	}

	m_Grains = indx;

	m_gpx = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gpx)
		return false;

	m_gpy = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gpy)
		return false;

	m_gps = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gps)
		return false;

	m_gvmi = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gvmi)
		return false;

	m_gdia = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gdia)
		return false;

	m_gv = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gv)
		return false;

	indx = 0;
	fseek(pFilePAT, 0, SEEK_SET);
	while (NULL != fgets(szBuffer, FILE_READ_BUFFER, pFilePAT))
	{
		sscanf(szBuffer, "%f %f %f %f", &aindx, &m_gpx[indx], &m_gpy[indx], &m_gps[indx]);
		m_gvmi[indx] = m_gps[indx] * m_asz * 1.0e-21f;
		m_gdia[indx] = sqrt(4.0f * m_gps[indx] / PI);
		m_gv[indx] = sqrt(gvm / m_gvmi[indx]);
		indx++;
	}

	fclose(pFilePAT);
	pFilePAT = nullptr;

	printf("Input grain area file finished \n");

	m_gymax = -1.0e4f;
	float gymin = 1.0e4f;
	float gxmax = -1.0e4f;
	float gxmin = 1.0e4f;

	for (int i = 0; i < m_Grains; i++)
	{
		if (m_gpy[i] > m_gymax)		m_gymax = m_gpy[i];
		if (m_gpy[i] < gymin)		gymin = m_gpy[i];
		if (m_gpx[i] > gxmax)		gxmax = m_gpx[i];
		if (m_gpx[i] < gxmin)		gxmin = m_gpx[i];
	}

	float gxmid = (gxmax - gxmin) / 2.0f;
	float gymid = (m_gymax - gymin) / 2.0f;

	printf(" \n");
	printf("X-coord of the Grains: %g %g\n", gxmin, gxmax);
	printf("Y-coord of the Grains: %g %g\n", gymin, m_gymax);

	printf("Number of grains = %d\n", m_Grains);

	m_gsx = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gsx)
		return false;

	m_gsy = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_gsy)
		return false;

	for (int i = 0; i < m_Grains; i++)
	{
		m_gsx[i] = m_gpx[i] - gxmin + m_offset;
		m_gsy[i] = m_gpy[i] - m_gymax;
	}

	float gxofs = -1.0f * gxmin - (float)m_ntpx / 2.0f;
	m_kxofs = (int)(gxofs + 0.01f);

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Setup TC distribution
	///////////////////////////////////////////////////////////////////////////////////////////////

	int idum = -13;
	int jdum = -17;
	int nrnd = (nrm + 3) * 38;

	for (int i = 0; i < nrnd; i++)
		ran1(idum);

	for (int i = 0; i < nrnd; i++)
		gasdev(jdum);

	m_tc = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_tc)
		return false;

	m_hks = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_hks)
		return false;

	FILE *pFileGP = fopen("grainpars.txt", "wt");
	fprintf(pFileGP, "ik,jz,tcbulk,gdia,gtc,rtc,tc,dtcK_ots,sigmatc,sigmatc_gs,sigmatc_ots,hks,sigmahk,vhkp\n");

	for (int i = 0; i < m_Grains; i++)
	{
		float ahkp = 10.0f;
		while (abs(ahkp) > 5.0f)
			ahkp = gasdev(jdum);

		float vhkp;
		for (int z = 0; z < m_ZCount; z++)
		{
			vhkp = 10.0f;
			if (abs(vhkp) > 5.0f)
				vhkp = gasdev(jdum);

			m_hks[(m_ZCount * i) + z] = (1.0f + m_sigmahk[z] * ahkp) * (1.0f + m_vhksigma * vhkp);
		}

		float atmp = 10.0f;
		while (abs(atmp) > 5.0f)
			atmp = gasdev(jdum);

		for (int z = 0; z < m_ZCount; z++)
		{
			float gtc = m_tcbulk[z] * (1.0f - m_gstc_d * pow(m_gdia[i], (-1.0f / m_gstc_v)));
			float rtc = m_sigmatc_ots[z] * atmp;
			m_tc[(m_ZCount * i) + z] = gtc * (1.0f + rtc);

			fprintf(pFileGP, "%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", (i + 1), (z + 1), m_tcbulk[z], m_gdia[i], gtc, rtc, m_tc[(m_ZCount * i) + z], 0.0f, m_sigmatc[z], m_sigmatc_gs[z], m_sigmatc_ots[z], m_hks[(m_ZCount * i) + z], m_sigmahk[z], vhkp);
		}
	}

	fclose(pFileGP);

	printf(" \n");
	printf(" nqt %d\n", m_nqs3);
	printf(" \n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// The following input file describes edge positions of all grains
	// �ȉ��̓��̓t�@�C���ɂ́A���ׂĂ̌������̃G�b�W�ʒu���L�q����Ă���B
	///////////////////////////////////////////////////////////////////////////////////////////////

	pFilePAT = fopen(pat_gp, "rt");

	indx = 0;
	float gnv = 0.0, dMax = 0.0;

	while (NULL != fgets(szBuffer, FILE_READ_BUFFER, pFilePAT))
	{
		// �ԍ��ƒ��_��
		sscanf(szBuffer, "%f %f", &aindx, &gnv);
		if (dMax < gnv)
			dMax = gnv;

		// ���_�����[�h
		for (int i = 0; i < (int)gnv; i++)
			fgets(szBuffer, FILE_READ_BUFFER, pFilePAT);
	}

	if ((int)(aindx + 0.01) != m_Grains)
	{
		printf("Number of grains does not match %d %g\n", m_Grains, (aindx + 0.01));
		fclose(pFilePAT);
		return false;
	}
	m_nnv = std::make_unique<int[]>(m_Grains);
	if (nullptr == m_nnv)
		return false;

	m_Vertexs = (int)(dMax + 0.01);

	m_vx = std::make_unique<float[]>(m_Grains * m_Vertexs);
	if (nullptr == m_vx)
		return false;

	m_vy = std::make_unique<float[]>(m_Grains * m_Vertexs);
	if (nullptr == m_vy)
		return false;

	indx = 0;
	fseek(pFilePAT, 0, SEEK_SET);
	while (NULL != fgets(szBuffer, FILE_READ_BUFFER, pFilePAT))
	{
		// �ԍ��ƒ��_��
		sscanf(szBuffer, "%f %f", &aindx, &gnv);

		int nnv = (int)(gnv + 0.01f);
		m_nnv[(int)(aindx + 0.01f) - 1] = nnv;

		// ���_
		for (int inv = 0; inv < nnv; inv++)
		{
			float tvx = 0.0f, tvy = 0.0f;
			fgets(szBuffer, FILE_READ_BUFFER, pFilePAT);
			sscanf(szBuffer, "%f %f", &tvx, &tvy);

			m_vx[indx] = tvx;
			m_vy[indx] = tvy;
			indx++;
		}
	}
	fclose(pFilePAT);

	printf(" \n");
	printf("Input grain edges file finished\n");
	printf(" \n");

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Set up movie output file
	///////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Initial magnetization patterns
	///////////////////////////////////////////////////////////////////////////////////////////////

	m_ymx = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_ymx)
		return false;

	m_ymy = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_ymy)
		return false;

	m_ymz = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_ymz)
		return false;

	indx = 0;
	switch (m_ichki)
	{
	// �����_��
	case 0:
		for (int i = 0; i < m_Grains; i++)
		{
			float theta = 0.0;
			float ra1 = ran1(idum);
			float ra2 = ran1(idum);
			float ra3 = ran1(idum);
			if (ra1 > 0.5)
				theta = (3.0f * ra2) * PI / 180.0f;
			else
				theta = (180.0f - 3.0f * ra2) * PI / 180.0f;

			float phi = ra3 * 2.0f * PI / 180.0f;

			float samx = sin(theta) * cos(phi);
			float samy = sin(theta) * sin(phi);
			float samz = cos(theta);

			for (int z = 0; z < m_ZCount; z++)
			{
				m_ymx[indx] = samx;
				m_ymy[indx] = samy;
				m_ymz[indx] = samz;
				indx++;
			}
		}
		break;

	// Positive
	case 1:
		for (int i = 0; i < m_Grains; i++)
		{
			for (int z = 0; z < m_ZCount; z++)
			{
				m_ymx[indx] = 0.0;
				m_ymy[indx] = 0.0;
				m_ymz[indx] = 1.0;
				indx++;
			}
		}
		break;

	// Negative
	default:
		for (int i = 0; i < m_Grains; i++)
		{
			for (int z = 0; z < m_ZCount; z++)
			{
				m_ymx[indx] = 0.0;
				m_ymy[indx] = 0.0;
				m_ymz[indx] = -1.0;
				indx++;
			}
		}
		break;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Set up simulatio time steps       
	//  Calculate total simulation time  
	///////////////////////////////////////////////////////////////////////////////////////////////

	float ypsiz = m_gymax - gymin + (float)m_ntpy;
	float totsim = m_pmove * ypsiz / m_spd;
	m_nht = (int)(totsim / m_dtsimi) + 1;
	printf("total sim, m_nht %g %d\n", totsim, m_nht);

	int npwskip = (int)((float)m_nht / (float)m_npwrite);
	int nmwskip = (int)((float)m_nht / (float)m_nmwrite);
	int ntwskip = (int)((float)m_nht / (float)m_ntwrite);

	if (1 == m_animation)
	{
		printf("\n");
		printf("animation log\n");
		printf("npwskip=%d -> %d\n", npwskip, (m_nht / npwskip));
		printf("nmwskip=%d -> %d\n", nmwskip, (m_nht / nmwskip));
		printf("ntwskip=%d -> %d\n", ntwskip, (m_nht / ntwskip));
		printf("\n");
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Set up head field
	///////////////////////////////////////////////////////////////////////////////////////////////

	m_ttw = std::make_unique<float[]>(m_nht);
	if (nullptr == m_ttw)
		return false;

	m_hhw = std::make_unique<float[]>(m_nht);
	if (nullptr == m_hhw)
		return false;

	// �������݃p�^�[��
	std::vector<size_t> vecPatternInitialLength;
	std::vector<std::vector<int>> vecPattern;
	vecPatternInitialLength.clear();
	vecPattern.clear();

	sprintf(szFile, "%s/%s", pszPath, FILE_PATTERN);
	FILE *pFilePattern = fopen(szFile, "rt");
	if (nullptr != pFilePattern)
	{
		// �n�߂�
		char szLine[1024];

		while (true)
		{
			szLine[0] = '\0';
			if (NULL == fgets(szLine, 1024, pFilePattern))
				break;

			std::vector<int> vec;

			int i = 0;
			while (('0' == szLine[i]) || ('1' == szLine[i]) || ('2' == szLine[i]))
			{
				vec.push_back(('0' == szLine[i]) ? 0: 1);
				i++;
			}

			if (0 == vec.size())
				break;

			vecPatternInitialLength.push_back(vec.size());

			// �{�`����
			szLine[0] = '\0';
			fgets(szLine, 1024, pFilePattern);

			i = 0;
			while (('0' == szLine[i]) || ('1' == szLine[i]) || ('2' == szLine[i]))
			{
				vec.push_back(('0' == szLine[i]) ? 0: 1);
				i++;
			}

			// �n�߂��Ȃ��{�`�����݂̂̏ꍇ
			if (vec.size() == vecPatternInitialLength[vecPatternInitialLength.size() - 1])
				vecPatternInitialLength[vecPatternInitialLength.size() - 1] = 0;

			vecPattern.push_back(vec);

			size_t size = vecPattern.size();
			printf("pattern %zd len:%zd initial: ", size, vecPattern[size - 1].size());
			size_t j = 0;
			for (; j < vecPatternInitialLength[size - 1]; j++)
				printf("%d", vecPattern[size - 1][j]);
			printf(" bit:");
			for (; j < vecPattern[size - 1].size(); j++)
				printf("%d", vecPattern[size - 1][j]);
			printf("\n");
		}

		fclose(pFilePattern);
	}

	float tsm = tdelay;
	float tsm_counter = tdelay;
	float sgni = -1.0f;
	float sgni_before = sgni;
	float hhw = 0.0;
	int nBitStep = -1;

	for (int i = 0; i < m_nht; i++)
	{
		m_ttw[i] = (float)i * m_dtsimi * m_spd + m_pstart;
		tsm += m_dtsimi;
		tsm_counter += m_dtsimi;

		if (tsm > tbit)
		{
			// �p�^�[���w�肪�����
			if (0 < vecPattern.size())
			{
				nBitStep++;

				// �p�^�[�����̈���Z����ΌJ��Ԃ��B
				if ((int)vecPattern[0].size() <= nBitStep)
					nBitStep = (int)vecPatternInitialLength[0];

				// ���]���������Ȃ���Ώ������΂��̂ł͂Ȃ�
				// ��ɔ��]�����Ă����Č��ɖ߂�����B
				if (0 == vecPattern[0][nBitStep])
					sgni = sgni * -1.0f;
			}

			tsm = 0.0f;
			sgni = sgni * -1.0f;

			// �r�b�g���]�̎��A�O��̒l����n�߂�
			// ���̍ہA�����ʒu�̋t���������Ă��܂������悤�ɂ��Ă���B
			if (sgni_before != sgni)
			{
				tsm_counter = 0.0f;
				hhw = m_hhw[i - 1] + sgni_before;
			}
			sgni_before = sgni;
		}

		m_hhw[i] = sgni * (2.0f * exp(-1.0f * tsm_counter / m_trise) - 1.0f) + hhw;

		if (1.0f < m_hhw[i])
			m_hhw[i] = 1.0f;
		else if (-1.0f > m_hhw[i])
			m_hhw[i] = -1.0f;
	}

	FILE *pf = fopen("hhw.csv", "wt");
	for (int i = 0; i < m_nht; i++)
		fprintf(pf, "%d,%f\n", (i + 1), m_hhw[i]);
	fclose(pf);

	// ���n�Ɏw��nm�Ńp�^�[�������
	if (1 < m_ichki)
	{
		for (int i = 0; i < m_Grains; i++)
		{
			for (int z = 0; z < m_ZCount; z++)
			{
				int nIndex = m_ZCount * i + z;
				if (0.0 <= m_gpy[i])
				{
					if (0 < (((int)m_gpy[i] / m_ichki) % 2))
					{
						m_ymx[nIndex] = 0.0f;
						m_ymy[nIndex] = 0.0f;
						m_ymz[nIndex] = 1.0f;
					}
					else
					{
						m_ymx[nIndex] = 0.0f;
						m_ymy[nIndex] = 0.0f;
						m_ymz[nIndex] = -1.0f;
					}
				}
				else
				{
					if (0 < (((int)abs(m_gpy[i]) / m_ichki) % 2))
					{
						m_ymx[nIndex] = 0.0f;
						m_ymy[nIndex] = 0.0f;
						m_ymz[nIndex] = -1.0f;
					}
					else
					{
						m_ymx[nIndex] = 0.0f;
						m_ymy[nIndex] = 0.0f;
						m_ymz[nIndex] = 1.0f;
					}
				}
			}
		}
	}

	// �t�@�C�����珉��������ǂݍ��݂܂��B
	sprintf(szFile, "%s/%s%s", pszPath, pat, EXTENSION_VEC_MZ);
	FILE *pFileMz = fopen(szFile, "rt");
	if (nullptr != pFileMz)
	{
		// �w���܂Ƃ߂��ŏI��Mz
		int nIndex = 0;
		while (NULL != fgets(szBuffer, FILE_READ_BUFFER, pFileMz))
		{
			if (0 == strlen(szBuffer))
				break;

			// �z��̓�����Ɠ����悤��Grain��Z���ɕ���ł���B
			m_ymz[nIndex] = (float)strtod(szBuffer, NULL);
			nIndex++;
		}

		fclose(pFileMz);
	}
	pFileMz = nullptr;

	// ���������̏�Ԃ��t�@�C���ɏo�͂��܂��B
/*	sprintf(szFile, "%s/%s/%s%s", pszPath, FOLDER_SIMON, m_fpfx, EXTENSION_KIN);
	FILE *pFileKin = fopen(szFile, "wt");
	if (nullptr != pFileKin)
	{
		for (int z = 0; z < m_ZCount; z++)
		{
			for (int i = 0; i < m_Grains; i++)
			{
				int nIndex = m_ZCount * i + z;
				fprintf(pFileKin, "%g %g %d %g %g %g\n", m_gpx[i], m_gpy[i], z, m_ymx[nIndex], m_ymy[nIndex], m_ymz[nIndex]);
			}
		}

		fclose(pFileKin);
	}
	pFileKin = nullptr;
*/
	// �Đ��p��Ms���t�@�C���o�͂��܂��B
	sprintf(szFile, "%s/%s/%s%s", pszPath, FOLDER_SIMON, m_fpfx, EXTENSION_MS);
	FILE *pFileMs = fopen(szFile, "wt");
	if (nullptr != pFileMs)
	{
		for (int z = 0; z < m_ZCount; z++)
		{
			for (int i = 0; i < m_Grains; i++)
				fprintf(pFileMs, "%d\n", m_amsrt[z]);
		}

		fclose(pFileMs);
	}
	pFileMs = nullptr;

	// �W�v�ɕK�v�ȃt�@�C���R�s�[

	try
	{
		sprintf(szFile, "%s/%s/%s%s", pszPath, FOLDER_SIMON, pat, EXTENSION_GA_DAT);
		std::filesystem::copy_file(pat_ga, szFile, std::filesystem::copy_options::overwrite_existing);
	}
	catch (std::filesystem::filesystem_error fe)
	{
		printf("**********************************************************************\n");
		printf("Failed to file copy. File : %s -> %s", pat_ga, szFile);
		printf("**********************************************************************\n");
	}

	try
	{
		sprintf(szFile, "%s/%s/%s%s", pszPath, FOLDER_SIMON, pat, EXTENSION_GP_DAT);
		std::filesystem::copy_file(pat_gp, szFile, std::filesystem::copy_options::overwrite_existing);
	}
	catch (std::filesystem::filesystem_error fe)
	{
		printf("**********************************************************************\n");
		printf("Failed to file copy. File : %s -> %s", pat_gp, szFile);
		printf("**********************************************************************\n");
	}

	try
	{
		sprintf(szFile, "%s/%s/%s", pszPath, FOLDER_SIMON, pszParameter);
		std::filesystem::copy_file(m_fpa_in, szFile, std::filesystem::copy_options::overwrite_existing);
	}
	catch (std::filesystem::filesystem_error fe)
	{
		printf("**********************************************************************\n");
		printf("Failed to file copy. File : %s -> %s", m_fpa_in, szFile);
		printf("**********************************************************************\n");
	}

	try
	{
		std::filesystem::path path(m_fimed);
		sprintf(szFile, "%s/%s/%s", pszPath, FOLDER_SIMON, path.filename().string().c_str());
		std::filesystem::copy_file(m_fimed, szFile, std::filesystem::copy_options::overwrite_existing);
	}
	catch (std::filesystem::filesystem_error fe)
	{
		printf("**********************************************************************\n");
		printf("Failed to file copy. File : %s -> %s", m_fpa_in, szFile);
		printf("**********************************************************************\n");
	}

	// �}���`�g���b�N�p
	VectorTrack vecTrack;
	vecTrack.clear();

	sprintf(szFile, "%s/%s", pszPath, FILE_OFFSET);
	FILE *pFileTrack = fopen(szFile, "rt");
	if (nullptr != pFileTrack)
	{
		fgets(szBuffer, FILE_READ_BUFFER, pFileTrack);

		while (NULL != fgets(szBuffer, FILE_READ_BUFFER, pFileTrack))
		{
			printf("Track %s", szBuffer);

			TRACK Track;

			// Name
			char *p = strtok(szBuffer, ",");
			if (NULL == p)
				break;
			strcpy(Track.szName, p);

			// Offset
			p = strtok(NULL, ",");
			if (NULL == p)
				break;
			Track.fOffset = (float)strtod(p, NULL);

			// Bitlength
			p = strtok(NULL, ",");
			if (NULL == p)
				break;
			Track.fBitlength = (float)strtod(p, NULL);

			vecTrack.push_back(Track);
		}

		fclose(pFileTrack);

		if (0 < vecTrack.size())
		{
			// ���n���c���Ă���
			sprintf(szFile, "%s/%s%s%s", pszPath, FILE_VEC, m_fpfx, m_fpfx);
			FILE *pFileVEC = fopen(szFile, "wt");
			if (nullptr != pFileVEC)
			{
				fprintf(pFileVEC, "%d %d %d\n", m_Grains, m_nht, m_nht);

				indx = 0;
				for (int i = 0; i < m_Grains; i++)
				{
					fprintf(pFileVEC, "%d %d %d\n", (i + 1), m_nnv[i], m_Grains);

					for (int j = 0; j < m_nnv[i]; j++)
					{
						fprintf(pFileVEC, "%g %g %g\n", m_vx[indx], m_vy[indx], 0.0);
						indx++;
					}
				}

				for (int i = 0; i < m_Grains; i++)
				{
					int nIndex = m_ZCount * i;
					fprintf(pFileVEC, "%g %g %g\n", m_ymx[nIndex], m_ymy[nIndex], m_ymz[nIndex]);
				}

				fclose(pFileVEC);
			}
			pFileVEC = nullptr;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	// Start Simulation
	//  Scan through each grains
	//  Each grain gets moved at speed: spd
	///////////////////////////////////////////////////////////////////////////////////////////////

	m_Yx = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_Yx)
		return false;

	m_Yy = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_Yy)
		return false;

	m_Yz = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_Yz)
		return false;

	m_YHx = std::make_unique<float[]>(m_Grains * m_ZCount * MAXORD);
	if (nullptr == m_YHx)
		return false;

	m_YHy = std::make_unique<float[]>(m_Grains * m_ZCount * MAXORD);
	if (nullptr == m_YHy)
		return false;

	m_YHz = std::make_unique<float[]>(m_Grains * m_ZCount * MAXORD);
	if (nullptr == m_YHz)
		return false;

	m_EWTx = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_EWTx)
		return false;

	m_EWTy = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_EWTy)
		return false;

	m_EWTz = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_EWTz)
		return false;

	m_SAVFx = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_SAVFx)
		return false;

	m_SAVFy = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_SAVFy)
		return false;

	m_SAVFz = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_SAVFz)
		return false;

	m_ACORx = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_ACORx)
		return false;

	m_ACORy = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_ACORy)
		return false;

	m_ACORz = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_ACORz)
		return false;

	m_tmpt = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_tmpt)
		return false;

	m_bms = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_bms)
		return false;

	m_hthi = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_hthi)
		return false;

	m_htx = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_htx)
		return false;

	m_hty = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_hty)
		return false;

	m_htz = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_htz)
		return false;

	m_hax = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_hax)
		return false;

	m_hay = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_hay)
		return false;

	m_haz = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_haz)
		return false;

	m_huk = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_huk)
		return false;

	m_hea = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_hea)
		return false;

	m_amk = std::make_unique<float[]>(m_Grains * m_ZCount);
	if (nullptr == m_amk)
		return false;

	m_amgx = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_amgx)
		return false;

	m_amgy = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_amgy)
		return false;

	m_amgz = std::make_unique<float[]>(m_Grains);
	if (nullptr == m_amgz)
		return false;

	m_amdz = std::make_unique<float[]>(m_Grains * m_nht);
	if (nullptr == m_amdz)
		return false;

	m_ztmp = std::make_unique<float[]>(m_Grains * m_nht);
	if (nullptr == m_ztmp)
		return false;

	m_zhed = std::make_unique<float[]>(m_Grains * m_nht);
	if (nullptr == m_zhed)
		return false;

	m_rs = std::make_unique<int[]>(m_Grains * m_ZCount);
	if (nullptr == m_rs)
		return false;

	m_ELCO = std::make_unique<float[]>(SIZE_ELCO_D1 * SIZE_ELCO_D2);
	if (nullptr == m_ELCO)
		return false;

	m_TESCO = std::make_unique<float[]>(SIZE_TESCO_D1 * SIZE_TESCO_D2);
	if (nullptr == m_TESCO)
		return false;

#if 0
	// �������Ōv�Z����ƂȂ������܂������Ȃ�
	// �z��̃C���f�b�N�X�����܂������Ă��Ȃ��B
	m_ELCO[0] = 1.0f;
	m_ELCO[1] = 1.0f;
	m_TESCO[0] = 0.0f;
	m_TESCO[1] = 2.0f;
	m_TESCO[SIZE_TESCO_D1 * 1] = 1.0f;
	m_TESCO[SIZE_TESCO_D1 * 11 + 2] = 0.0f;

	float PC[12];
	PC[0] = 1.0f;

	float RQFAC = 1.0f;

	for (int i = 1; i < 12; i++)
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
		for (int j = 0; j <= NQM1; j++)
		{
			int nIndex = NQP1 - j - 1;
			PC[nIndex] = PC[nIndex - 1] + FNQM1 * PC[nIndex];
		}

		PC[0] *= FNQM1;

		// Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------

		float PINT = PC[0];
		float XPIN = PC[0] / 2.0f;
		float TSIGN = 1.0f;

		for(int j = 1; j <= i; j++)
		{
			TSIGN = -TSIGN;
			PINT += (TSIGN * PC[j] / (float)(j + 1));
			XPIN += (TSIGN * PC[j] / (float)(j + 2));
		}

		// Store coefficients in ELCO and TESCO. --------------------------------

		m_ELCO[SIZE_ELCO_D1 * i] = PINT * RQ1FAC;
		m_ELCO[SIZE_ELCO_D1 + i + 1] = 1.0f;

		for (int j = 1; j <= i; j++)
			m_ELCO[SIZE_ELCO_D1 * i + j + 1] = RQ1FAC * PC[j] / (float)(j + 1);

		float AGAMQ = RQFAC * XPIN;
		float RAGQ = 1.0f / AGAMQ;
		m_TESCO[SIZE_TESCO_D1 + i + 1] = RAGQ;

		if (NQP1 < 12)
			m_TESCO[SIZE_TESCO_D1 * NQP1] = RAGQ * RQFAC / (float)(NQP1 + 1);

		m_TESCO[SIZE_TESCO_D1 * NQM1 + 2] = RAGQ;
	}
#endif

	int nResult = 0;

#ifdef CPU

	FILE *pFile = fopen("xxx.csv", "wt");

	fprintf(pFile, "ik,gpx(ik),gpy(ik),gsx(ik),gsy(ik),t,tf,iht,nht,ttw(iht),hhw(iht),kxofs,asx,asy,jz,tmpt(jz),bms(jz),hthi(jz),htz(jz),haz(jz),huk(jz),hea(jz),gms(jz),amk(jz),amz,tmpt(0),y(0)\n");

	m_ntmpw = 0;
	m_nwr = 0;

	for (int i = 0; i < m_Grains; i++)
//	for (int i = 0; i < 1; i++)
	{
		for (int z = 0; z < m_ZCount; z++)
		{
			const int nIndex = (m_ZCount * i) + z;
			m_Yx[z] = m_ymx[nIndex];
			m_Yy[z] = m_ymy[nIndex];
			m_Yz[z] = m_ymz[nIndex];
		}

		// get grain position
		float asx = m_gsx[i];
		float bsy = m_gsy[i];

		int itp = 0;
		int iwr = 0;

		for (int iht = 0; iht < m_nht; iht++)
//		for (int iht = 0; iht < 1; iht++)
		{
			float tsim = m_ttw[iht];
			float asy = bsy + tsim * m_spd + m_pstart;
			caltemp(asx, asy);

			if ((iht % npwskip) == 0)
			{
				itp++;
				m_ztmp[itp] = m_tmpt[0];
				m_zhed[itp] = m_haz[0];
			}

			float tmpt_1 = m_tmpt[0];

			float tmpc = 0.0f;
			for (int z = 0; z < m_ZCount; z++)
			{
				m_htx[z] = 0.0f;
				m_hty[z] = 0.0f;
				m_htz[z] = 0.0f;
				m_hax[z] = 0.0f;
				m_hay[z] = 0.0f;
				m_haz[z] = 0.0f;
				m_huk[z] = 0.0f;
				m_hea[z] = 0.0f;
				m_amk[z] = 0.0f;

				tmpc = m_tc[(m_ZCount * i) + z];

				if (m_tmpt[z] > (tmpc - m_dtc))
				{
					float theta = PI * ran1(idum);
					float phi = 2.0f * PI * ran1(idum);
					float ssp = sin(phi);
					float csp = cos(phi);
					float sst = sin(theta);
					float cst = cos(theta);
					m_Yx[z] = sst * csp;
					m_Yy[z] = sst * ssp;
					m_Yz[z] = cst;
				}
				else if (tmpt_1 >= m_tambient)
				{
					float tmptscling = m_alftscaling * m_tmpt[z];
					float alfbms = pow((1.0f - tmptscling / tmpc), (1.0f / 3.0f));
					float s_alfg = m_alfg[z] / alfbms * (1.0f - tmptscling / (3.0f * tmpc));

					float fTemp = 1.0f - (m_tmpt[z] / tmpc);

					float bms = pow(fTemp, (1.0f / 3.0f));
					m_hthi[z] = m_htho[z] * sqrt(s_alfg * m_tmpt[z] / bms) * m_gv[i];

					m_huk[z] = m_hko[z] * m_hks[(m_ZCount * i) + z] * pow(fTemp, m_ept[z]);
					m_hea[z] = m_heb[z] * pow(fTemp, m_eptAex[z]);

					float atmp = 10.0f;
					while (abs(atmp) > 3.0f)
						atmp = gasdev(jdum);

					float theta = PI * ran1(idum);
					float phi = 2.0f * PI * ran1(idum);

					float ssp = sin(phi);
					float csp = cos(phi);
					float sst = sin(theta);
					float cst = cos(theta);
					m_htx[z] = m_hthi[z] * atmp * sst * csp;
					m_hty[z] = m_hthi[z] * atmp * sst * ssp;
					m_htz[z] = m_hthi[z] * atmp * cst;

					m_hax[z] = m_hhw[iht] * m_hax[z];
					m_hay[z] = m_hhw[iht] * m_hay[z];
					m_haz[z] = m_hhw[iht] * m_haz[z];

					m_amk[z] = 1.0f;
				}
			}

			int ist = 1;
			float tf = m_dtsimi * GAMA * hnm;
			float dtf = tf / m_antf;

			for (float t = 0.0f; t < tf; t += dtf)
//			for (float t = 0.0f; t < tf; t = tf)
			{
				float to = t + dtf;
				if (false == slsode(t, to, ist))
				{
					printf("\n\nhere!\n\n");
					return false;
				}

				if (ist < 0)
				{
					printf("I was hoping not to end up here!\n");
					return false;
				}

				for (int z = 0; z < m_ZCount; z++)
				{
					float yamp = sqrt((m_Yx[z] * m_Yx[z]) + (m_Yy[z] * m_Yy[z]) + (m_Yz[z] * m_Yz[z]));
					float ydif = abs(yamp - 1.0f);
					if (ydif > 0.02f)
					{
						printf("\t%d\t%f\t%d\tYamp Anormny! %f\t%f\n", iht, t, z, yamp, ydif);
						ist = 1;
					}
					m_Yx[z] /= yamp;
					m_Yy[z] /= yamp;
					m_Yz[z] /= yamp;
				}
			}

			if (NULL != pFile)
			{
				int ik = 2;
				int jz = 26;

				if (0.0f != m_amk[0])
				{
					float amx = 0.0, amy = 0.0, amz = 0.0;
					myavg(amx, amy, amz);

					fprintf(pFile, "%d,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%d,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", ik, m_gpx[ik], m_gpy[ik], m_gsx[ik], m_gsy[ik], tf, tf, iht, m_nht, m_ttw[iht], m_hhw[iht], m_kxofs, asx, asy, jz, m_tmpt[jz], m_bms[jz], m_hthi[jz], m_htz[jz], m_haz[jz], m_huk[jz], m_hea[jz], 0.0, m_amk[jz], amz, m_tmpt[0], m_Yz[0]);
				}
			}

			// finish one m_dtsimi step

			if ((iht % nmwskip) == 0)
			{
				iwr++;
				float amx = 0.0f, amy = 0.0f, amz = 0.0f;
				myavg(amx, amy, amz);
				m_amdx[iwr] = amx;
				m_amdy[iwr] = amy;
				m_amdz[iwr] = amz;
			}

			if ((iht % ntwskip) == 0)
				printf("%f %f %d %d\n", tsim, totsim, i, m_Grains);
		}

		printf("Motion completed for Grain %d\n", i);

		// Finish all steps for a single grain

		float amx = 0.0f, amy = 0.0f, amz = 0.0f;
		myavg(amx, amy, amz);
		m_amgx[i] = amx;
		m_amgy[i] = amy;
		m_amgz[i] = amz;

		if (0 == m_ntmpw)
			m_ntmpw = itp;
		if (0 == m_nwr)
			m_nwr = iwr;

		for (int z = 0; z < m_ZCount; z++)
		{
			const int nIndex = (m_ZCount * i) + z;
			m_ymx[nIndex] = m_Yx[z];
			m_ymy[nIndex] = m_Yy[z];
			m_ymz[nIndex] = m_Yz[z];
		}

	}

	fclose(pFile);


#else

	printf("GPU Mode\n");

	CCUDAInterface* pInterface = new CCUDAInterface;
	if (nullptr != pInterface)
	{
		// GPU���g�p�ł��邩
		printf("GPU Check\n");
		if (cudaSuccess == pInterface->CheckGPU())
		{
			// �������m��
			printf("GPU Memory\n");
			size_t nMemorySize = 0;
			if (cudaSuccess != (nResult = pInterface->Malloc(m_ntpx, m_ntpy, m_ntpz, m_ZCount, m_nhfx, m_nhfy, m_nhfz, m_Grains, m_nht, nmwskip, npwskip, nMemorySize)))
			{
				// �G���[��
				printf("Failed to get GPU memory. %d\n", nResult);
				printf("Get GPU Memory %zd MB (%zd byte)\n", (nMemorySize / 1024 / 1024), nMemorySize);
			}
			else
			{
				printf("Get GPU Memory %zd MB (%zd byte)\n", (nMemorySize / 1024 / 1024), nMemorySize);
				printf("GPU Copy\n");

				// �e�X���b�h�����������闐���̃V�[�h
				srand(nrnd);
				for (int i = 0; i < (m_Grains * m_ZCount); i++)
					m_rs[i] = (rand() % 0xFF) + (rand() % 0xFF);
//					m_rs[i] = 0;

				int nSuccess = 0;
				if (cudaSuccess == pInterface->Upload2(m_Thickness.get(), m_tpf.get()))									nSuccess++;
				if (cudaSuccess == pInterface->Upload4(m_hfx.get(), m_hfy.get(), m_hfz.get()))							nSuccess++;
				if (cudaSuccess == pInterface->Upload5(m_ept.get(), m_eptAex.get(), m_alfg.get()))						nSuccess++;
				if (cudaSuccess == pInterface->Upload6(m_heb.get(), m_hko.get()))										nSuccess++;
				if (cudaSuccess == pInterface->Upload7(m_gv.get(), m_htho.get()))										nSuccess++;
				if (cudaSuccess == pInterface->Upload9(m_gsx.get(), m_gsy.get()))										nSuccess++;
				if (cudaSuccess == pInterface->Upload10(m_tc.get(), m_hks.get()))										nSuccess++;
				if (cudaSuccess == pInterface->Upload11(m_ymx.get(), m_ymy.get(), m_ymz.get()))							nSuccess++;
				if (cudaSuccess == pInterface->Upload12(m_ttw.get(), m_hhw.get()))										nSuccess++;
				if (cudaSuccess == pInterface->Upload13(m_rs.get()))													nSuccess++;
				if (cudaSuccess == pInterface->Upload14(m_ELCO.get(), m_TESCO.get()))									nSuccess++;

				if (11 != nSuccess)
				{
					printf("**********************************************************************\n");
					printf("Failed to copy. %d", nSuccess);
					printf("**********************************************************************\n");
				}
				else
				{
					const int nTrackLoop = (0 == vecTrack.size()) ? 1: (int)vecTrack.size();
					for (int l = 0; l < nTrackLoop; l++)
					{
						char szAddName[64] = "\0";
						float fOffset = 0.0f;

						// Bitlength�ɂ�菑�����݃p�^�[���̍Đݒ�
						if (0 < vecTrack.size())
						{
							const int nPatternIndex = (l < (int)vecPattern.size()) ? l: ((int)vecPattern.size() - 1);
							tbit = vecTrack[l].fBitlength / m_spd;
							tsm = tdelay;
							tsm_counter = tdelay;
							sgni = -1.0f;
							sgni_before = sgni;
							hhw = 0.0;
							nBitStep = -1;

							for (int i = 0; i < m_nht; i++)
							{
								tsm += m_dtsimi;
								tsm_counter += m_dtsimi;

								if (tsm > tbit)
								{
									if (0 < vecPattern.size())
									{
										nBitStep++;

										// �p�^�[�����̈���Z����ΌJ��Ԃ��B
										if ((int)vecPattern[nPatternIndex].size() <= nBitStep)
											nBitStep = (int)vecPatternInitialLength[nPatternIndex];

										if (0 == vecPattern[nPatternIndex][nBitStep])
											sgni = sgni * -1.0f;
									}

									tsm = 0.0;
									sgni = sgni * -1.0f;

									// �r�b�g���]�̎��A�O��̒l����n�߂�
									// ���̍ہA�����ʒu�̋t���������Ă��܂������悤�ɂ��Ă���B
									if (sgni_before != sgni)
									{
										tsm_counter = 0.0f;
										hhw = m_hhw[i - 1] + sgni_before;
									}
									sgni_before = sgni;
								}

								m_hhw[i] = sgni * (2.0f * exp(-1.0f * tsm_counter / m_trise) - 1.0f) + hhw;

								if (1.0f < m_hhw[i])
									m_hhw[i] = 1.0f;
								else if (-1.0f > m_hhw[i])
									m_hhw[i] = -1.0f;
							}

							pInterface->Upload12(m_ttw.get(), m_hhw.get());

							szAddName[0] = '_';
							strcpy(&szAddName[1], vecTrack[l].szName);

							fOffset = vecTrack[l].fOffset;
						}

						printf("GPU Start %d\n", cudaGetLastError());

						float tf = m_dtsimi * GAMA * hnm;
						float dtf = tf / m_antf;
						nResult = pInterface->Run((nLogGrain - 1), m_kxofs, m_dtc, m_tambient, tf, m_alftscaling, dtf, rumach(), fOffset, m_nps, m_MP_edge, m_tpf_edge, m_animation);

						printf("GPU %d\n", nResult);

						pInterface->Download1(m_ymx.get(), m_ymy.get(), m_ymz.get());

						if (1 == m_animation)
							pInterface->Download3(m_amdz.get(), m_ztmp.get(), m_zhed.get());

						for (int i = 0; i < m_Grains; i++)
						{
							for (int z = 0; z < m_ZCount; z++)
							{
								m_Yx[z] = m_ymx[m_ZCount * i + z];
								m_Yy[z] = m_ymy[m_ZCount * i + z];
								m_Yz[z] = m_ymz[m_ZCount * i + z];
							}

							myavg(m_amgx[i], m_amgy[i], m_amgz[i]);
						}

// 9009
						printf("All grains finished !!! %d\n", m_Grains);

						// �����ƌv�Z�o�������̂݃��O���o���܂��B
						if (cudaSuccess != nResult)
						{
							printf("**********************************************************************\n");
							printf("Failed to calculate. %d", nResult);
							printf("**********************************************************************\n");
						}
						else
						{
							// ����Track�������݂̏ꍇ
							if (0 < vecTrack.size())
							{
								// �W�~�[�搶�̏o��
								sprintf(szFile, "%s/%s%s%s", pszPath, FILE_VEC, m_fpfx, szAddName);
								FILE *pFileVEC = fopen(szFile, "wt");
								if (nullptr != pFileVEC)
								{
									fprintf(pFileVEC, "%d %d %d\n", m_Grains, m_nht, m_nht);

									indx = 0;
									for (int i = 0; i < m_Grains; i++)
									{
										fprintf(pFileVEC, "%d %d %d\n", (i + 1), m_nnv[i], m_Grains);

										for (int j = 0; j < m_nnv[i]; j++)
										{
											fprintf(pFileVEC, "%g %g %g\n", m_vx[indx], m_vy[indx], 0.0);
											indx++;
										}
									}

									for (int i = 0; i < m_Grains; i++)
									{
										fprintf(pFileVEC, "%g %g %g\n", m_amgx[i], m_amgy[i], m_amgz[i]);
									}

									fclose(pFileVEC);
								}
								pFileVEC = nullptr;

								// �A�j���[�V�����p�t�@�C��
								if (1 == m_animation)
								{
									// ����
#ifdef ALL_ANIMATION
									const int nwr = m_nht / m_nmwrite + 1;
									const int nSkipM = m_nmwrite;
#else
									const int nwr = m_nht / nmwskip + 1;
									const int nSkipM = 1;
#endif

									sprintf(szFile, "%s/%s%s%s", pszPath, FILE_DYM, m_fpfx, szAddName);
									FILE *pFileDYM = fopen(szFile, "wt");
									if (nullptr != pFileDYM)
									{
										fprintf(pFileDYM, "%d %d %d\n", m_Grains, nwr, m_nht);

										indx = 0;
										for (int i = 0; i < m_Grains; i++)
										{
											fprintf(pFileDYM, "%d %d %d\n", (i + 1), m_nnv[i], m_Grains);

											for (int j = 0; j < m_nnv[i]; j++)
											{
												fprintf(pFileDYM, "%g %g %g\n", m_vx[indx], m_vy[indx], 0.0);
												indx++;
											}
										}

										for (int j = 0; j < nwr; j++)
										{
											fprintf(pFileDYM, "%d %d %d\n", (j + 1), nwr, m_Grains);

											const int nIndex = m_Grains * j * nSkipM;
											for (int i = 0; i < m_Grains; i++)
												fprintf(pFileDYM, "0 0 %g\n", m_amdz[nIndex + i]);
										}

										fclose(pFileDYM);
									}
									pFileDYM = nullptr;

									// �M
									sprintf(szFile, "%s/%s%s%s", pszPath, FILE_TMP, m_fpfx, szAddName);
									FILE *pFileTMP = fopen(szFile, "wt");
									if (nullptr != pFileTMP)
									{
#ifdef ALL_ANIMATION
									const int ntmpw = m_nht / m_npwrite + 1;
									const int nSkipT = m_npwrite;
#else
									const int ntmpw = m_nht / npwskip + 1;
									const int nSkipT = 1;
#endif

										fprintf(pFileTMP, "%d %d %d\n", m_Grains, ntmpw, m_nht);

										indx = 0;
										for (int i = 0; i < m_Grains; i++)
										{
											fprintf(pFileTMP, "%d %d %d\n", (i + 1), m_nnv[i], m_Grains);

											for (int j = 0; j < m_nnv[i]; j++)
											{
												fprintf(pFileTMP, "%g %g %g\n", m_vx[indx], m_vy[indx], 0.0);
												indx++;
											}
										}

										for (int j = 0; j < ntmpw; j++)
										{
											fprintf(pFileTMP, "%d %d %d\n", (j + 1), ntmpw, m_Grains);

											const int nIndex = m_Grains * j * nSkipT;
											for (int i = 0; i < m_Grains; i++)
												fprintf(pFileTMP, "%g 0 %g\n", m_zhed[nIndex + i], m_ztmp[nIndex + i]);
										}

										fclose(pFileTMP);
									}
									pFileTMP = nullptr;
								}

								// �T�C�����搶�p�̍ŏI���ʃt�@�C��
								sprintf(m_fiend, "%s/%s/%s%s%s", pszPath, FOLDER_SIMON, m_fpfx, szAddName, EXTENSION_END);
								FILE *pFileEnd = fopen(m_fiend, "wt");
								if (nullptr != pFileEnd)
								{
									for (int z = 0; z < m_ZCount; z++)
									{
										for (int i = 0; i < m_Grains; i++)
										{
											int nIndex = m_ZCount * i + z;
											fprintf(pFileEnd, "%g %g %d %g %g %g\n", m_gpx[i], m_gpy[i], z, m_ymx[nIndex], m_ymy[nIndex], m_ymz[nIndex]);
										}
									}

									fclose(pFileEnd);
								}
								pFileEnd = nullptr;
							}
						}

#ifdef GPU_LOG
						std::unique_ptr<float[]> asx = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> asy = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> bsy = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> alft = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> amk = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> hhwl = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> tmpt = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> htz = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> hax = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> hay = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> haz = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> huk = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> hea = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> bms = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> yy = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> yz = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> amy = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> amz = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> tmpc = std::make_unique<float[]>(m_ZCount * m_nht);
						std::unique_ptr<float[]> hthi = std::make_unique<float[]>(m_ZCount * m_nht);

						pInterface->Download2(asx.get(), asy.get(), bsy.get(), alft.get(), amk.get(), hhwl.get(), tmpt.get(), htz.get(), hax.get(), hay.get(), haz.get(), huk.get(), hea.get(), bms.get(), yy.get(), yz.get(), amy.get(), amz.get(), tmpc.get(), hthi.get());

						sprintf(szFile, "%s/Spdym_G_ik%d.csv", pszPath, nLogGrain);
						FILE* pFile = fopen(szFile, "wt");
						if (nullptr != pFile)
						{
							fprintf(pFile, "asx_G,asy_G,bsy_G,iht_G,ik_G,jz_G,alfg_G,alft_G,amk(jz)_G,hhw(jz)_G,tmpt(jz)_G,htz(jz)_G,hax(jz)_G,hay(jz)_G,haz(jz)_G,huk(jz)_G,hea(jz)_G,heb(jz)_G,gms(jz)_G,y(jmy)_G,y(jmz)_G,amy_G,amz_G,tmpc_G,hks(ikjz)_G,hthi(jz)_G\n");

							for (int i = 0; i < m_nht; i++)
							{
								for (int j = 0; j < m_ZCount; j++)
								{
									const int nIndex = i * m_ZCount + j;
									if (0.0f != tmpc[nIndex])
									{
										fprintf(pFile, "%f,%f,%f,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
												asx[nIndex], asy[nIndex], bsy[nIndex], (i + 1), nLogGrain, (j + 1),
												alft[nIndex], alft[nIndex], 
												amk[nIndex], hhwl[nIndex], tmpt[nIndex], htz[nIndex], hax[nIndex], hay[nIndex], haz[nIndex],
												huk[nIndex], hea[nIndex], m_heb[j], bms[nIndex],
												yy[nIndex], yz[nIndex], amy[nIndex], amz[nIndex],
												tmpc[nIndex], m_hks[(nLogGrain - 1) * m_ZCount + j], hthi[nIndex]);
									}
								}
							}

							fclose(pFile);
						}
#endif
					}
				}
			}

			pInterface->Free();
			pInterface->ReleaseGPU();
		}

		delete pInterface;
	}
	pInterface = nullptr;

#endif

	printf("\nOutput log file.\n");

	// �����ƌv�Z�o�������̂݃��O���o���܂��B
	if (cudaSuccess != nResult)
	{
		printf("**********************************************************************\n");
		printf("Failed to calculate. %d", nResult);
		printf("**********************************************************************\n");
	}
	else
	{
		// �W�~�[�搶�̏o��
		sprintf(szFile, "%s/%s%s", pszPath, FILE_VEC, m_fpfx);
		FILE *pFileVEC = fopen(szFile, "wt");
		if (nullptr != pFileVEC)
		{
			fprintf(pFileVEC, "%d %d %d\n", m_Grains, m_nht, m_nht);

			indx = 0;
			for (int i = 0; i < m_Grains; i++)
			{
				fprintf(pFileVEC, "%d %d %d\n", (i + 1), m_nnv[i], m_Grains);

				for (int j = 0; j < m_nnv[i]; j++)
				{
					fprintf(pFileVEC, "%g %g %g\n", m_vx[indx], m_vy[indx], 0.0);
					indx++;
				}
			}

			for (int i = 0; i < m_Grains; i++)
			{
				fprintf(pFileVEC, "%g %g %g\n", m_amgx[i], m_amgy[i], m_amgz[i]);
			}

			fclose(pFileVEC);
		}
		pFileVEC = nullptr;

		// �A�j���[�V�����p�t�@�C��
		if (1 == m_animation)
		{
			// ����
#ifdef ALL_ANIMATION
			const int nwr = m_nht / m_nmwrite + 1;
			const int nSkipM = m_nmwrite;
#else
			const int nwr = m_nht / nmwskip + 1;
			const int nSkipM = 1;
#endif

			sprintf(szFile, "%s/%s%s", pszPath, FILE_DYM, m_fpfx);
			FILE *pFileDYM = fopen(szFile, "wt");
			if (nullptr != pFileDYM)
			{
				fprintf(pFileDYM, "%d %d %d\n", m_Grains, nwr, m_nht);

				indx = 0;
				for (int i = 0; i < m_Grains; i++)
				{
					fprintf(pFileDYM, "%d %d %d\n", (i + 1), m_nnv[i], m_Grains);

					for (int j = 0; j < m_nnv[i]; j++)
					{
						fprintf(pFileDYM, "%g %g %g\n", m_vx[indx], m_vy[indx], 0.0);
						indx++;
					}
				}

				for (int j = 0; j < nwr; j++)
				{
					fprintf(pFileDYM, "%d %d %d\n", (j + 1), nwr, m_Grains);

					const int nIndex = m_Grains * j * nSkipM;
					for (int i = 0; i < m_Grains; i++)
						fprintf(pFileDYM, "0 0 %g\n", m_amdz[nIndex + i]);
				}

				fclose(pFileDYM);
			}
			pFileDYM = nullptr;

			// �M
#ifdef ALL_ANIMATION
			const int ntmpw = m_nht / m_npwrite + 1;
			const int nSkipT = m_npwrite;
#else
			const int ntmpw = m_nht / npwskip + 1;
			const int nSkipT = 1;
#endif

			sprintf(szFile, "%s/%s%s", pszPath, FILE_TMP, m_fpfx);
			FILE *pFileTMP = fopen(szFile, "wt");
			if (nullptr != pFileTMP)
			{
				fprintf(pFileTMP, "%d %d %d\n", m_Grains, ntmpw, m_nht);

				indx = 0;
				for (int i = 0; i < m_Grains; i++)
				{
					fprintf(pFileTMP, "%d %d %d\n", (i + 1), m_nnv[i], m_Grains);

					for (int j = 0; j < m_nnv[i]; j++)
					{
						fprintf(pFileTMP, "%g %g %g\n", m_vx[indx], m_vy[indx], 0.0);
						indx++;
					}
				}

				for (int j = 0; j < ntmpw; j++)
				{
					fprintf(pFileTMP, "%d %d %d\n", (j + 1), ntmpw, m_Grains);

					const int nIndex = m_Grains * j * nSkipT;
					for (int i = 0; i < m_Grains; i++)
						fprintf(pFileTMP, "%g 0 %g\n", m_zhed[nIndex + i], m_ztmp[nIndex + i]);
				}

				fclose(pFileTMP);
			}
			pFileTMP = nullptr;
		}

		// �T�C�����搶�p�̍ŏI���ʃt�@�C��
		sprintf(m_fiend, "%s/%s/%s%s", pszPath, FOLDER_SIMON, m_fpfx, EXTENSION_END);
		FILE *pFileEnd = fopen(m_fiend, "wt");
		if (nullptr != pFileEnd)
		{
			for (int z = 0; z < m_ZCount; z++)
			{
				for (int i = 0; i < m_Grains; i++)
				{
					int nIndex = m_ZCount * i + z;
					fprintf(pFileEnd, "%g %g %d %g %g %g\n", m_gpx[i], m_gpy[i], z, m_ymx[nIndex], m_ymy[nIndex], m_ymz[nIndex]);
				}
			}

			fclose(pFileEnd);
		}
		pFileEnd = nullptr;
	}

	// �v�Z�̏�
	// �����^�]�p
	sprintf(szFile, "%s/%s/%s%s", pszPath, FOLDER_SIMON, m_fpfx, EXTENSION_CON);
	FILE *pFileCON = fopen(szFile, "wt");
	if (nullptr != pFileCON)
	{
		fprintf(pFileCON, "%s%s\n", pat, EXTENSION_GP_DAT);
		fprintf(pFileCON, "%s%s\n", pat, EXTENSION_GA_DAT);
		fprintf(pFileCON, "%s\n", pszParameter);

		strftime(szBuffer, sizeof(szBuffer), "%Y%m%d%H%M%S", localtime(&timeStart));
		fprintf(pFileCON, "%s\n", szBuffer);

		time_t timeEnd = time(NULL);
		strftime(szBuffer, sizeof(szBuffer), "%Y%m%d%H%M%S", localtime(&timeEnd));
		fprintf(pFileCON, "%s\n", szBuffer);

		fprintf(pFileCON, "%d\n", m_nht);
		fprintf(pFileCON, "%d\n", m_ntpy);
		fprintf(pFileCON, "%d\n", m_Grains);

		fclose(pFileCON);
	}
	pFileCON = nullptr;

	printf("\nComplete.\n");

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// ������Z���ς��Ƃ�
void CHAMR::myavg(float &amx, float &amy, float &amz)
{
	amx = 0.0;
	amy = 0.0;
	amz = 0.0;

	for (int z = 0; z < m_ZCount; z++)
	{
		amx += m_Yx[z];
		amy += m_Yy[z];
		amz += m_Yz[z];
	}

	amx /= (float)m_ZCount;
	amy /= (float)m_ZCount;
	amz /= (float)m_ZCount;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// ��������
float CHAMR::ran1(int &IDUM)
{
	static int IFF = 0;
	static float R[97];
	static int IX1 = 0, IX2 = 0, IX3 = 0;

	const int M1 = 259200, IA1 = 7141, IC1 = 54773;
	const int M2 = 134456, IA2 = 8121, IC2 = 28411;
	const int M3 = 243000, IA3 = 4561, IC3 = 51349;
	const float RM1 = 3.8580247E-6f, RM2 = 7.4373773E-6f;


	if ((IDUM < 0) || (IFF == 0))
	{
		IFF = 1;
		IX1 = (IC1 - IDUM) % M1;
		IX1 = (IA1 * IX1 + IC1) % M1;
		IX2 = IX1 % M2;
		IX1 = (IA1 * IX1 + IC1) % M1;
		IX3 = IX1 % M3;

		for (int J = 0; J < 97; J++)
		{
			IX1 = (IA1 * IX1 + IC1) % M1;
			IX2 = (IA2 * IX2 + IC2) % M2;
			R[J] = ((float)IX1 + (float)IX2 * RM2) * RM1;
		}

		IDUM = 1;
	}

	IX1 = (IA1 * IX1 + IC1) % M1;
	IX2 = (IA2 * IX2 + IC2) % M2;
	IX3 = (IA3 * IX3 + IC3) % M3;
	int J = (97 * IX3) / M3;
	if ((J > 96) || (J < 0))
		printf("Something Wrong\n");

	float ran1 = R[J];
	R[J] = ((float)IX1 + (float)IX2 * RM2) * RM1;

	return ran1;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// gamma=0.0176 (1/Oe)(1/naosecond)
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
// �K�E�X (���K) ����
float CHAMR::gasdev(int &idum)
{
	static int iset = 0;
	static float gset = 0.0f;
	float gasdev = 0.0f;

	if (idum < 0)
		iset = 0;

	if (iset == 0)
	{
		float v1 = 0.0, v2 = 0.0, rsq = 0.0;
		do
		{
			v1 = 2.0f * ran1(idum) - 1.0f;
			v2 = 2.0f * ran1(idum) - 1.0f;
			rsq = pow(v1, 2.0f) + pow(v2, 2.0f);
		} while ((rsq >= 1.0f) || (rsq == 0.0f));

		float fac = sqrt(-2.0f * log(rsq) / rsq);
		gset = v1 * fac;
		gasdev = v2 * fac;
		iset = 1;
	}
	else
	{
		gasdev = gset;
		iset = 0;
	}

	return gasdev;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// �m���
float CHAMR::rumach()
{
	float U = 1.0f;
	do
	{
		U *= 0.5f;
	} while ((1.0f + U) != 1.0f);

	// 1.19209290e-07
	return U * 2.0f;
}


#ifdef CPU

///////////////////////////////////////////////////////////////////////////////////////////////
// �v�Z���x
void CHAMR::caltemp(float asx, float asy)
{
//	int ksy = (int)(asy + 1.0) - 1;
//	int ksx = (int)(asx + 1.0) - m_kxofs - 1;
	int ksy = (int)asy;
	int ksx = (int)asx - m_kxofs;

	// �������x
	for (int j = 0; j < m_ZCount; j++)
		m_tmpt[j] = 275.0;
//		m_tmpt[j] = 611.0;
//		m_tmpt[j] = 650.0;
//		m_tmpt[j] = 700.0;

	// �v�Z�ΏۂȂ牷�x������B
	if ((asy >= 0.0f) && (ksy < m_ntpy) && (ksx >= 0) && (ksx < m_ntpx))
	{
		for (int z = 0; z < m_ZCount; z++)
		{
			m_tmpt[z] = m_tpf[(m_ZCount * m_ntpy * ksx) + (m_ZCount * ksy) + z];

			if (ksy > 0)
			{
				float tp = m_tmpt[z];
				float tn = m_tpf[(m_ZCount * m_ntpy * ksx) + (m_ZCount * (ksy - 1)) + z];
//				m_tmpt[z] = tn + ((tp - tn) * ((asy - 1.0f - (float)in) / (float)(ksy - in)));
				m_tmpt[z] = tn + ((tp - tn) * (asy - (float)ksy));
			}

			m_hax[z] = m_hay[z] = m_haz[z] = 0.0f;

			// 3D�̏ꍇ1���傫���B
			if ((1 < m_nhfy) && (DUMMY_DOUBLE != m_nps))
			{
				//	nps - (MP_edge - tpf_edge)
				int ksyh = ksy + (int)(m_nps - m_MP_edge - m_tpf_edge);

				if ((0 < ksyh) && (ksyh < m_nhfx))
				{
					m_hax[z] = m_hfx[(m_ZCount * m_nhfy * ksx) + (m_ZCount * ksyh) + z];
					m_hay[z] = m_hfy[(m_ZCount * m_nhfy * ksx) + (m_ZCount * ksyh) + z];
					m_haz[z] = m_hfz[(m_ZCount * m_nhfy * ksx) + (m_ZCount * ksyh) + z];

					if (1 < ksyh)
					{
						float tp = m_hax[z];
						float tn = m_hfx[(m_ZCount * m_nhfy * ksx) + (m_ZCount * (ksyh - 1)) + z];
						m_hax[z] = tn + ((tp - tn) * (asy - (float)ksyh));

						tp = m_hay[z];
						tn = m_hfy[(m_ZCount * m_nhfy * ksx) + (m_ZCount * (ksyh - 1)) + z];
						m_hay[z] = tn + ((tp - tn) * (asy - (float)ksyh));

						tp = m_haz[z];
						tn = m_hfz[(m_ZCount * m_nhfy * ksx) + (m_ZCount * (ksyh - 1)) + z];
						m_haz[z] = tn + ((tp - tn) * (asy - (float)ksyh));
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
// LLG
void CHAMR::fun(float YDx[], float YDy[], float YDz[])
{
	float hex, hey, hez;
	float hanx, hany, hanz;
	float ahx, ahy, ahz;
	float agx, agy, agz;
	float spx, spy, spz;
	float agm;

	m_nLLG++;

	for (int z = 0; z < m_ZCount; z++)
	{
		int zn = z - 1;
		int zp = z + 1;

		if (z == 0)
		{
			// ��ԏ�̑w
			hex = m_hea[z] * m_Yx[zp];
			hey = m_hea[z] * m_Yy[zp];
			hez = m_hea[z] * m_Yz[zp];
		}
		else if(z == (m_ZCount - 1))
		{
			// ��ԉ��̑w
			hex = m_hea[zn] * m_Yx[zn];
			hey = m_hea[zn] * m_Yy[zn];
			hez = m_hea[zn] * m_Yz[zn];
		}
		else
		{
			// ���ԑw
			hex = (m_hea[z] * m_Yx[zp]) + (m_hea[zn] * m_Yx[zn]);
			hey = (m_hea[z] * m_Yy[zp]) + (m_hea[zn] * m_Yy[zn]);
			hez = (m_hea[z] * m_Yz[zp]) + (m_hea[zn] * m_Yz[zn]);
		}

		hanx = 0.0;
		hany = 0.0;
		hanz = m_huk[z] * m_Yz[z];

		ahx = hanx + m_hax[z] + hex;
		ahy = hany + m_hay[z] + hey;
		ahz = hanz + m_haz[z] + hez;

		agx = m_alfg[z] * (ahx + m_htx[z]);
		agy = m_alfg[z] * (ahy + m_hty[z]);
		agz = m_alfg[z] * (ahz + m_htz[z]);

		spx = (ahy * m_Yz[z]) - (ahz * m_Yy[z]);
		spy = (ahz * m_Yx[z]) - (ahx * m_Yz[z]);
		spz = (ahx * m_Yy[z]) - (ahy * m_Yx[z]);

		agm = (agx * m_Yx[z]) + (agy * m_Yy[z]) + (agz * m_Yz[z]);

		YDx[z] = m_amk[z] * (spx + agx - (agm * m_Yx[z]));
		YDy[z] = m_amk[z] * (spy + agy - (agm * m_Yy[z]));
		YDz[z] = m_amk[z] * (spz + agz - (agm * m_Yz[z]));

		//if (0 == z)
		//{
		//	printf("\n%e\t%e\t%e\n", hex, hey, hez);
		//	printf("%e\t%e\t%e\n", hanx, hany, hanz);
		//	printf("%e\t%e\t%e\n", m_hax[z], m_hay[z], m_haz[z]);
		//	printf("%e\t%e\t%e\t%e\t%e\n", m_hea[z], m_huk[z], hanz, m_haz[z], hez);
		//	printf("%e\t%e\t%e\t%e\t%e\n", m_alfg[z], ahz, m_htz[z], ahx, ahy);
		//	printf("%e\t%e\t%e\n", agx, agy, agz);
		//	printf("%e\t%e\t%e\t%e\n", m_amk[z], spz, agz, agm);
		//	printf("%e\t%e\t%e\n", m_Yx[z], m_Yy[z], m_Yz[z]);
		//	printf("%e\t%e\t%e\n\n", YDx[z], YDy[z], YDz[z]);
		//}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
bool CHAMR::slsode(float T, float TOUT, int &ISTATE)
{
	// Declare all other variables.
	float H0, TDIST, TOL, TOLSF, SUM, W0;

	bool bCalculating = true;

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

	if (1 == ISTATE)
	{
		// ����
		if (TOUT == T)
			return true;

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

		m_UROUND = rumach();
		m_TN = T;
		m_JSTART = 0;
		s_NHNIL = 0;
		m_NST = 0;
		s_NSLAST = 0;
		m_HU = 0.0f;
		m_NQU = 0;

		// Initial call to F.  (LF0 points to YH(*,2).) -------------------------

		fun((m_YHx.get() + m_ZCount), (m_YHy.get() + m_ZCount), (m_YHz.get() + m_ZCount));

		// Load the initial value vector in YH. ---------------------------------

		for (int z = 0; z < m_ZCount; z++)
		{
			m_YHx[z] = m_Yx[z];
			m_YHy[z] = m_Yy[z];
			m_YHz[z] = m_Yz[z];
		}

		// Load and invert the EWT array.  (H is temporarily set to 1.0.) -------

		m_NQ = 1;
		m_H = 1.0f;
		sewset();

		for (int z = 0; z < m_ZCount; z++)
		{
			if (0.0f >= m_EWTx[z])		return false;
			if (0.0f >= m_EWTy[z])		return false;
			if (0.0f >= m_EWTz[z])		return false;

			m_EWTx[z] = 1.0f / m_EWTx[z];
			m_EWTy[z] = 1.0f / m_EWTy[z];
			m_EWTz[z] = 1.0f / m_EWTz[z];
		}

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
			TDIST = abs(TOUT - T);
			W0 = std::max<float>(abs(T), abs(TOUT));

			if (TDIST < (2.0f * m_UROUND * W0))
			{
				return false;
			}

			TOL = RTO;

			TOL = std::max<float>(TOL, (100.0f * m_UROUND));
			TOL = std::min<float>(TOL, 0.001f);

			SUM = svnorm((m_YHx.get() + m_ZCount), (m_YHy.get() + m_ZCount), (m_YHz.get() + m_ZCount), m_EWTx.get(), m_EWTy.get(), m_EWTz.get());
			SUM = 1.0f / (TOL * W0 * W0) + TOL * pow(SUM, 2.0f);

			H0 = 1.0f / sqrt(SUM);
			H0 = std::min(H0, TDIST);
			H0 = abs(H0) * ((0.0f <= (TOUT - T)) ? 1.0f: -1.0f);
		}

		// Load H with H0 and scale YH(*,2) by H0. ------------------------------

		m_H = H0;
//		printf("H %e\n", m_H);

		for (int z = 0; z < m_ZCount; z++)
		{
			m_YHx[m_ZCount + z] *= H0;
			m_YHy[m_ZCount + z] *= H0;
			m_YHz[m_ZCount + z] *= H0;
		}
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

		s_NSLAST = m_NST;

		if (((m_TN - TOUT) * m_H) >= 0.0f)
		{
			if (0 != sintdy(TOUT))
			{
				return false;
			}

			bCalculating = false;
		}
	}

	while (bCalculating)
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
		if ((m_NST - s_NSLAST) >= MXSTP0)
		{
			return false;
		}

		sewset();

		for (int z = 0; z < m_ZCount; z++)
		{
			if (0.0f >= m_EWTx[z])	return false;
			if (0.0f >= m_EWTy[z])	return false;
			if (0.0f >= m_EWTz[z])	return false;

			m_EWTx[z] = 1.0f / m_EWTx[z];
			m_EWTy[z] = 1.0f / m_EWTy[z];
			m_EWTz[z] = 1.0f / m_EWTz[z];
		}

// 270
		TOLSF = m_UROUND * svnorm(m_YHx.get(), m_YHy.get(), m_YHz.get(), m_EWTx.get(), m_EWTy.get(), m_EWTz.get());

		if (TOLSF > 1.0f)
		{
			TOLSF = TOLSF * 2.0f;

			return false;
		}

		if ((m_TN + m_H) == m_TN)
		{
			s_NHNIL++;
			if (s_NHNIL <= MXHNL0)
			{
//				printf("%f\t%f\t%f\t\t%d\t%d\n", m_TN, m_H, m_TN, s_NHNIL, MXHNL0);

				printf("SLSODE-  Warning..internal T (=R1) and H (=R2) are\n");
				printf("      such that in the machine, T + H = T on the next step  \n");
				printf("      (H = step size). Solver will continue anyway\n");

				if (s_NHNIL >= MXHNL0)
				{
					printf("SLSODE-  Above warning has been issued I1 times.  \n");
					printf("      It will not be issued again for this problem\n");
				}
			}
		}

		// -----------------------------------------------------------------------
		//  CALL SSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,SPREPJ,SSOLSY)
		// -----------------------------------------------------------------------
//printf("%f\n", m_YHz[0]);
		sstode();
//printf("%f\n", m_Yz[1]);

		if (1 != (1 - m_KFLAG))
		{
			return false;
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
//		printf("\t%d\t%e\t%e\t%e\t%e\t%e\n", m_nLLG, m_TN, TOUT, (m_TN - TOUT), m_H, m_Yx[0]);
		if (((m_TN - TOUT) * m_H) < 0.0)
		{
			continue;
		}

		if (0 != sintdy(TOUT))
			return false;

//printf("%f\n", m_Yz[0]);

		bCalculating = false;
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
	// ISTATE �� 2 �ɐݒ肳��A�I�v�V�����o�͂����[�N�E�A���C�Ƀ��[�h�����B���[�N�E�A���C�Ƀ��[�h�����B
	// -----------------------------------------------------------------------
// 420
	ISTATE = 2;

//printf("%f\n", m_Yz[0]);

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
// SCFODE is called by the integrator routine to set coefficients needed there.
// The coefficients for the current method, as given by the value of METH, are set for all orders and saved.
// The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
// (A smaller value of the maximum order is also allowed.)
// SCFODE is called once at the beginning of the problem, and is not called again unless and until METH is changed.
// 
// The ELCO array contains the basic method coefficients.
// The coefficients el(i), 1 .le. i .le. nq+1, for the method of order nq are stored in ELCO(i,nq).
// They are given by a genetrating polynomial, i.e.,
//     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
// For the implicit Adams methods, l(x) is given by
//     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
// For the BDF methods, l(x) is given by
//     l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
// where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
// 
// The TESCO array contains test constants used for the local error test and the selection of step size and/or order.
// At order nq, TESCO(k,nq) is used for the selection of step size at order nq - 1 if k = 1, at order nq if k = 2, and at order nq + 1 if k = 3.
// 
// SCFODE �͐ϕ����[�`���ŕK�v�ȌW����ݒ肷�邽�߂ɌĂяo����܂��B
// METH �̒l�ɂ���ė^�����錻�݂̃��\�b�h�̌W���́A���ׂĂ̎����ɑ΂��Đݒ肳��A�ۑ������B
// �����őz�肳���ő原���� METH = 1 �Ȃ� 12�AMETH = 2 �Ȃ� 5 �ł��B
// (�ő原����菬�����l���������)�B
// SCFODE �͖��̊J�n���Ɉ�x�Ăяo����AMETH ���ύX����Ȃ�����A�܂��ύX�����܂ōēx�Ăяo����邱�Ƃ͂���܂���B
// 
// ELCO �z��͊�{���\�b�h�̌W�����܂݂܂��B
// ���� nq �̃��\�b�h�̌W�� el(i), 1 .le. i .le. nq+1 �� ELCO(i,nq) �Ɋi�[����܂��B
// �����͈�`�������ŗ^������A
//     l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
// �A�IAdams�@�ł́Al(x)�͎����ŗ^������B
//     dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1), l(-1) = 0.
// BDF�@�ł́Al(x)�͎����ŗ^������B
//     l(x) = (x+1)*(x+2)* ... *(x+nq)/K�ŗ^������A
// ������ K = factorial(nq)*(1 + 1/2 + ... + 1/nq) �ł���B
// 
// TESCO �z��́A�Ǐ��덷����ƃX�e�b�v�T�C�Y�����/�܂��͎����̑I���Ɏg�p����錟��萔���܂ށB
// ���� nq �ɂ����āATESCO(k,nq) �́Ak = 1 �̏ꍇ�͎��� nq - 1 �ŁAk = 2 �̏ꍇ�͎��� nq �ŁAk = 3 �̏ꍇ�͎��� nq + 1 �ŁA�X�e�b�v�T�C�Y�̑I���Ɏg�p�����B
// 
void CHAMR::scfode()
{
	s_ELCO[0][0] = 1.0f;
	s_ELCO[0][1] = 1.0f;
	s_TESCO[0][0] = 0.0f;
	s_TESCO[0][1] = 2.0f;
	s_TESCO[1][0] = 1.0f;
	s_TESCO[11][2] = 0.0f;

	float PC[12];
	PC[0] = 1.0f;

	float RQFAC = 1.0f;

	for (int i = 1; i < 12; i++)
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
		for (int j = 0; j <= NQM1; j++)
		{
			int nIndex = NQP1 - j - 1;
			PC[nIndex] = PC[nIndex - 1] + FNQM1 * PC[nIndex];
		}

		PC[0] *= FNQM1;

		// Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------

		float PINT = PC[0];
		float XPIN = PC[0] / 2.0f;
		float TSIGN = 1.0f;

		for(int j = 1; j <= i; j++)
		{
			TSIGN = -TSIGN;
			PINT += (TSIGN * PC[j] / (float)(j + 1));
			XPIN += (TSIGN * PC[j] / (float)(j + 2));
		}

		// Store coefficients in ELCO and TESCO. --------------------------------

		s_ELCO[i][0] = PINT * RQ1FAC;
		s_ELCO[i][1] = 1.0f;

		for (int j = 1; j <= i; j++)
			s_ELCO[i][j + 1] = RQ1FAC * PC[j] / (float)(j + 1);

		float AGAMQ = RQFAC * XPIN;
		float RAGQ = 1.0f / AGAMQ;
		s_TESCO[i][1] = RAGQ;

		if (NQP1 < 12)
			s_TESCO[NQP1][0] = RAGQ * RQFAC / (float)(NQP1 + 1);

		s_TESCO[NQM1][2] = RAGQ;
	}
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
int CHAMR::sintdy(float T)
{
	// C��1�Œ�̂��ߌv�Z������폜

	int IFLAG = 0;

	float TP = m_TN - m_HU - 100.0f * m_UROUND * (m_TN + m_HU);

	if (((T - TP) * (T - m_TN)) > 0.0f)
	{
// 90
		IFLAG = -2;
		return IFLAG;
	}

	float S = (T - m_TN) / m_H;

//printf("s %f\t%f\t%d\n", m_Yz[0], m_YHz[m_ZCount * m_L], m_L);

// 15
	int nOffset = m_ZCount * m_L;
	for (int z = 0; z < m_ZCount; z++)
	{
		m_Yx[z] = m_YHx[nOffset + z];
		m_Yy[z] = m_YHy[nOffset + z];
		m_Yz[z] = m_YHz[nOffset + z];
	}

//printf("s %f\n", m_Yz[0]);

	for (int i = 0; i <= m_NQ; i++)
	{
// 35
		nOffset = m_ZCount * (m_NQ - i);
		for (int z = 0; z < m_ZCount; z++)
		{
			m_Yx[z] = m_YHx[nOffset + z] + (S * m_Yx[z]);
			m_Yy[z] = m_YHy[nOffset + z] + (S * m_Yy[z]);
			m_Yz[z] = m_YHz[nOffset + z] + (S * m_Yz[z]);
		}
	}

	return IFLAG;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
bool CHAMR::sstode()
{
	int IREDO, IRET, M, NCF, NEWQ;
	float DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP, R, RH, RHDN, RHSM, RHUP, TOLD;

	int nCalculating = 0;

//printf("SSTODE start\n");

	m_KFLAG = 0;
	TOLD = m_TN;
	NCF = 0;
	DELP = 0.0f;

	DEL = 0.0f;

	IRET = 0;

	if (m_JSTART == 0)
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

		m_NQ = 0;
		m_L = 1;
		s_IALTH = 2;
		s_RMAX = 10000.0f;
		m_RC = 0.0f;
		m_EL0 = 1.0f;
		s_CRATE = 0.7f;
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
		scfode();

		nCalculating = 150;
	}

	while (690 > nCalculating)
	{

		if (150 == nCalculating)
		{
// 150
			for (int i = 0; i <= m_L; i++)
			{
				s_EL[i] = s_ELCO[m_NQ][i];
//				printf("%f\n", s_EL[i]);
			}
//			printf("aaaaaaaaaaaaaaaaaaaaa\n");

			m_RC = m_RC * s_EL[0] / m_EL0;
			m_EL0 = s_EL[0];
			s_CONIT = 0.5f / (float)(m_NQ + 2 + 1);

			nCalculating = (2 == IRET) ? 170: 0;
		}

		if (170 == nCalculating)
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
			RH = std::max<float>(RH, 0.0f);

// 175
			RH = std::min<float>(RH, s_RMAX);
			R = 1.0;

			for (int i = 1; i <= m_L; i++)
			{
				R *= RH;

				int nOffset = m_ZCount * i;
				for (int z = 0; z < m_ZCount; z++)
				{
					m_YHx[nOffset + z] *= R;
					m_YHy[nOffset + z] *= R;
					m_YHz[nOffset + z] *= R;
				}
			}

			m_H *= RH;
			m_RC = m_RC * RH;
			s_IALTH = m_L + 1;

//printf("SSTODE 180 H %e\n", m_H);

			if (IREDO == 0)
			{
// 690
				s_RMAX = 10.0;
				nCalculating = 700;
				continue;
			}

			nCalculating = 0;
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
		m_TN += m_H;
//		printf("SSTODE 200 TN %e\t%e\t%d\t%d\n", m_TN, m_H, m_NQ, m_L);
//		printf("SS %e\t%e\t%e\n", m_YHx[0],m_YHx[m_ZCount], m_YHx[m_ZCount*2]);

		for (int i = 0; i <= m_NQ; i++)
		{
			for (int j = (m_NQ - i); j <= m_NQ; j++)
			{
				int nOffset = m_ZCount * (j + 1);
				for (int z = 0; z < m_ZCount; z++)
				{
					m_YHx[nOffset - m_ZCount + z] += m_YHx[nOffset + z];
					m_YHy[nOffset - m_ZCount + z] += m_YHy[nOffset + z];
					m_YHz[nOffset - m_ZCount + z] += m_YHz[nOffset + z];
				}
			}
		}

//		printf("SS %e\t%e\t%e\n", m_YHx[0],m_YHx[m_ZCount], m_YHx[m_ZCount*2]);

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
		for (int z = 0; z < m_ZCount; z++)
		{
			m_Yx[z] = m_YHx[z];
			m_Yy[z] = m_YHy[z];
			m_Yz[z] = m_YHz[z];
		}

		fun(m_SAVFx.get(), m_SAVFy.get(), m_SAVFz.get());

		// -----------------------------------------------------------------------
		// If indicated, the matrix P = I - h*el(1)*J is reevaluated and
		// preprocessed before starting the corrector iteration.  IPUP is set
		// to 0 as an indicator that this has been done.
		// -----------------------------------------------------------------------
		// �w�����ꂽ�ꍇ�AP = I - h*el(1)*J�s�񂪍ĕ]������A�␳��̔������J�n����O�ɑO�������s����B ���ĕ]������, �␳��̔������J�n����O�ɑO���������.
		// IPUP �� �� 0 ���Z�b�g�����B
		// -----------------------------------------------------------------------

		for (int z = 0; z < m_ZCount; z++)
		{
			m_ACORx[z] = 0.0;
			m_ACORy[z] = 0.0;
			m_ACORz[z] = 0.0;
		}


		while (0 == nCalculating)
		{
			//-----------------------------------------------------------------------
			// In the case of functional iteration, update Y directly from
			// the result of the last function evaluation.
			// -----------------------------------------------------------------------
			// �֐������̏ꍇ�A�Ō�̊֐��]������C���璼��Y���X�V����B�𒼐ڍX�V���܂��B
			// -----------------------------------------------------------------------
// 270
//printf("%f\t%f\t%f\t%f\n", m_Yz[0], m_SAVFz[0], m_YHz[m_ZCount], m_ACORz[0]);
//			printf("%e\t%e\t%e\t%e\t%e\n", m_YHz[0], m_YHz[m_ZCount], m_YHz[m_ZCount*2], m_YHz[m_ZCount*3], m_YHz[m_ZCount*4]);
//			printf("A %e\t%e\t%e\t%e\t%e\t%d\n", m_SAVFz[0], m_ACORz[0], m_Yz[0], m_EWTz[0], m_YHz[m_ZCount], s_IALTH);

//			DEL = svnorm(m_Yx.get(), m_Yy.get(), m_Yz.get(), m_EWTx.get(), m_EWTy.get(), m_EWTz.get());
//			printf("DEL\t%f\t%e\t%e\t%e\n", DEL, m_H, m_SAVFz[0], m_YHz[m_ZCount]);

			for (int z = 0; z < m_ZCount; z++)
			{
				m_SAVFx[z] = m_H * m_SAVFx[z] - m_YHx[m_ZCount + z];
				m_SAVFy[z] = m_H * m_SAVFy[z] - m_YHy[m_ZCount + z];
				m_SAVFz[z] = m_H * m_SAVFz[z] - m_YHz[m_ZCount + z];

				m_Yx[z] = m_SAVFx[z] - m_ACORx[z];
				m_Yy[z] = m_SAVFy[z] - m_ACORy[z];
				m_Yz[z] = m_SAVFz[z] - m_ACORz[z];
			}
//printf("%f\t%f\n", m_Yz[0], m_SAVFz[0]);

			DEL = svnorm(m_Yx.get(), m_Yy.get(), m_Yz.get(), m_EWTx.get(), m_EWTy.get(), m_EWTz.get());
//			printf("B %e\t%e\t%e\t%e\t%e\t%e\t%e\n", DEL, m_SAVFz[0], m_ACORz[0], m_Yz[0], m_EWTz[0], m_YHz[m_ZCount], m_H);

			for (int z = 0; z < m_ZCount; z++)
			{
				m_Yx[z] = m_YHx[z] + s_EL[0] * m_SAVFx[z];
				m_Yy[z] = m_YHy[z] + s_EL[0] * m_SAVFy[z];
				m_Yz[z] = m_YHz[z] + s_EL[0] * m_SAVFz[z];

				m_ACORx[z] = m_SAVFx[z];
				m_ACORy[z] = m_SAVFy[z];
				m_ACORz[z] = m_SAVFz[z];
			}
//printf("%f\t%f\t%f\t%f\t%f\n", m_Yz[0], m_YHz[0], s_EL[0], m_SAVFz[0], m_ACORz[0]);

			// -----------------------------------------------------------------------
			// Test for convergence.  If M.gt.0, an estimate of the convergence
			// rate constant is stored in CRATE, and this is used in the test.
			// -----------------------------------------------------------------------
			// �����̃e�X�g�B
			// M.gt.0�̏ꍇ�A�������萔�̐���l��CRATE�Ɋi�[����A�e�X�g�Ɏg�p�����B�̐���l�� CRATE �Ɋi�[����A���ꂪ�e�X�g�Ɏg�p�����B
			// -----------------------------------------------------------------------
// 400
			if (M != 0)
				s_CRATE = std::max<float>((0.2f * s_CRATE), (DEL / DELP));

			DCON = DEL * std::min<float>(1.0f, (1.5f * s_CRATE)) / (s_TESCO[m_NQ][1] * s_CONIT);
//			printf("DCON %e\t%e\t%f\t%f\t%f\t%d\n", DCON, DEL, s_CRATE, s_TESCO[m_NQ][1], s_CONIT, m_NQ);

			if (DCON > 1.0f)
			{
				M++;

				if (M >= MAXCOR)
				{
					nCalculating = 430;
					break;
				}

				if ((M > 1) && (DEL > (2.0f * DELP)))
				{
					nCalculating = 430;
					break;
				}

				DELP = DEL;

				fun(m_SAVFx.get(), m_SAVFy.get(), m_SAVFz.get());

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
			if (M == 0)
				DSM = DEL / s_TESCO[m_NQ][1];

			if (M > 0)
				DSM = svnorm(m_ACORx.get(), m_ACORy.get(), m_ACORz.get(), m_EWTx.get(), m_EWTy.get(), m_EWTz.get()) / s_TESCO[m_NQ][1];

			if (DSM <= 1.0f)
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

				m_KFLAG = 0;
				IREDO = 0;
				m_NST++;
				m_HU = m_H;
				m_NQU = m_NQ;

				for (int i = 0; i <= m_L; i++)
				{
					int nOffset = m_ZCount * i;
					for (int z = 0; z < m_ZCount; z++)
					{
						m_YHx[nOffset + z] += (s_EL[i] * m_ACORx[z]);
						m_YHy[nOffset + z] += (s_EL[i] * m_ACORy[z]);
						m_YHz[nOffset + z] += (s_EL[i] * m_ACORz[z]);
					}
				}

				s_IALTH--;

				if (0 == s_IALTH)
				{
// 520
					RHUP = 0.0f;
					if (m_L != MAXORD)
					{
						int nOffset = m_ZCount * (MAXORD - 1);
						for (int z = 0; z < m_ZCount; z++)
						{
							m_SAVFx[z] = m_ACORx[z] - m_YHx[nOffset + z];
							m_SAVFy[z] = m_ACORy[z] - m_YHy[nOffset + z];
							m_SAVFz[z] = m_ACORz[z] - m_YHz[nOffset + z];
						}

						DUP = svnorm(m_SAVFx.get(), m_SAVFy.get(), m_SAVFz.get(), m_EWTx.get(), m_EWTy.get(), m_EWTz.get()) / s_TESCO[m_NQ][2];
						EXUP = 1.0f / (m_L + 1 + 1);
						RHUP = 1.0f / (1.4f * pow(DUP, EXUP) + 0.0000014f);
					}
				}
				else
				{
					if ((s_IALTH <= 1) && (m_L != MAXORD))
					{
						int nOffset = m_ZCount * (MAXORD - 1);
						for (int z = 0; z < m_ZCount; z++)
						{
							m_YHx[nOffset + z] = m_ACORx[z];
							m_YHy[nOffset + z] = m_ACORy[z];
							m_YHz[nOffset + z] = m_ACORz[z];
						}
					}

					nCalculating = 700;
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
				m_KFLAG--;
				m_TN = TOLD;

				for (int i = 0; i <= m_NQ; i++)
				{
					for (int j = (m_NQ - i); j <= m_NQ; j++)
					{
						int nOffset = m_ZCount * (j + 1);

						for (int z = 0; z < m_ZCount; z++)
						{
							m_YHx[nOffset - m_ZCount + z] -= m_YHx[nOffset + z];
							m_YHy[nOffset - m_ZCount + z] -= m_YHy[nOffset + z];
							m_YHz[nOffset - m_ZCount + z] -= m_YHz[nOffset + z];
						}
					}
				}

				s_RMAX = 2.0f;
				if (abs(m_H) <= 0.0f)
				{
					return false;
				}

				if (m_KFLAG <= -3)
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
					if (m_KFLAG == -10)
					{
						return false;
					}

					RH = 0.1f;
					RH = std::max(0.0f, RH);
					m_H *= RH;

//printf("SSTODE 640 H %e\n", m_H);

					for (int z = 0; z < m_ZCount; z++)
					{
						m_Yx[z] = m_YHx[z];
						m_Yy[z] = m_YHy[z];
						m_Yz[z] = m_YHz[z];
					}

					fun(m_SAVFx.get(), m_SAVFy.get(), m_SAVFz.get());

// 650
					for (int z = 0; z < m_ZCount; z++)
					{
						m_YHx[m_ZCount + z] = m_H * m_SAVFx[z];
						m_YHy[m_ZCount + z] = m_H * m_SAVFy[z];
						m_YHz[m_ZCount + z] = m_H * m_SAVFz[z];
					}

					s_IALTH = 5;

					if (0 == m_NQ)
					{
						// Next 200
						break;
					}

					m_NQ = 0;
					m_L = 1;
					IRET = 3;

					nCalculating = 150;
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
			EXSM = 1.0f / (float)(m_L + 1);
			RHSM = 1.0f / (1.2f * pow(DSM, EXSM) + 0.0000012f);
			RHDN = 0.0f;

			if (0 != m_NQ)
			{
				int nOffset = m_ZCount * m_L;
				DDN = svnorm((m_YHx.get() + nOffset), (m_YHy.get() + nOffset), (m_YHz.get() + nOffset), m_EWTx.get(), m_EWTy.get(), m_EWTz.get()) / s_TESCO[m_NQ][0];

				EXDN = 1.0f / (float)(m_NQ + 1);
				RHDN = 1.0f / (1.3f * pow(DDN, EXDN) + 0.0000013f);
//				printf("SSTODE 540 %f\t%f\t%f\t%f\t%f\t%f\n", RHSM, RHUP, RHDN, DDN, EXDN, RHDN);
			}

// 560
			if ((RHSM < RHUP) && (RHUP > RHDN))
			{
// 590
				NEWQ = m_L;
				RH = RHUP;
//				printf("SSTODE 590 %f\n", RH);

				if (RH < 1.1f)
				{
// 610
					s_IALTH = 3;
					nCalculating = 700;
					break;
				}

				R = s_EL[m_L] / (float)(m_L + 1);

				int nOffset = m_ZCount * (NEWQ + 1);
				for (int z = 0; z < m_ZCount; z++)
				{
					m_YHx[nOffset + z] = m_ACORx[z] * R;
					m_YHy[nOffset + z] = m_ACORy[z] * R;
					m_YHz[nOffset + z] = m_ACORz[z] * R;
				}
			}
			else
			{
// 570
				if (RHSM >= RHUP)
				{
					NEWQ = m_NQ;
					RH = RHSM;
//					printf("SSTODE 570 %f\n", RH);
				}
				else
				{
// 580
					NEWQ = m_NQ - 1;
					RH = RHDN;
//					printf("SSTODE 580 %f\n", RH);

					if ((m_KFLAG < 0) && (RH > 1.0f))
						RH = 1.0;
				}

// 620
//				printf("SSTODE 620 %d\t%d\n", NEWQ, m_NQ);
				if ((m_KFLAG == 0) && (RH < 1.1f))
				{
					s_IALTH = 3;
					nCalculating = 700;
					break;
				}

				if (m_KFLAG <= -2)
					RH = std::min<float>(RH, 0.2f);

				// -----------------------------------------------------------------------
				// If there is a change of order, reset NQ, l, and the coefficients.
				// In any case H is reset according to RH and the YH array is rescaled.
				// Then exit from 690 if the step was OK, or redo the step otherwise.
				// -----------------------------------------------------------------------
				// �I�[�_�[���ύX���ꂽ�ꍇ�́ANQ�Al�A�W�������Z�b�g����B
				// ������̏ꍇ���AH �� RH �ɏ]���ă��Z�b�g����AYH �z��͍ăX�P�[�����O�����B
				// �X�e�b�v�ɖ�肪�Ȃ����690���I�����A�����łȂ���΃X�e�b�v����蒼���B
				// -----------------------------------------------------------------------

				if (NEWQ == m_NQ)
				{
					nCalculating = 170;
					break;
				}
			}

// 630
//			printf("SSTODE 630\n");
			m_NQ = NEWQ;
			m_L = m_NQ + 1;
			IRET = 2;

			nCalculating = 150;
		}

		switch (nCalculating)
		{
		case 430:		break;
		default:		continue;
		}

		nCalculating = 0;

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
		s_RMAX = 2.0f;
		m_TN = TOLD;
//printf("SSTODE 430 TM %e\n", m_TN);

		for (int i = 0; i <= m_NQ; i++)
		{
			for (int j = (m_NQ - i); j <= m_NQ; j++)
			{
				int nOffset = m_ZCount * (j + 1);

				for (int z = 0; z < m_ZCount; z++)
				{
					m_YHx[nOffset - m_ZCount + z] -= m_YHx[nOffset + z];
					m_YHy[nOffset - m_ZCount + z] -= m_YHy[nOffset + z];
					m_YHz[nOffset - m_ZCount + z] -= m_YHz[nOffset + z];
				}
			}
		}

// 445
		if (abs(m_H) < 0.0f)
		{
			return false;
		}

		if (NCF == MXNCF)
		{
			return false;
		}

		RH = 0.25f;
		IREDO = 1;

		nCalculating = 170;
	}

// 700
	R = 1.0f / s_TESCO[m_NQU][1];

	for (int z = 0; z < m_ZCount; z++)
	{
		m_ACORx[z] *= R;
		m_ACORy[z] *= R;
		m_ACORz[z] *= R;
	}

// 720
	m_JSTART = 1;

//printf("%f\n", m_Yz[0]);

	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
void CHAMR::sewset()
{
	for (int z = 0; z < m_ZCount; z++)
	{
		m_EWTx[z] = RTO * abs(m_YHx[z]) + ATO;
		m_EWTy[z] = RTO * abs(m_YHy[z]) + ATO;
		m_EWTz[z] = RTO * abs(m_YHz[z]) + ATO;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
// 
float CHAMR::svnorm(float Vx[], float Vy[], float Vz[], float Wx[], float Wy[], float Wz[])
{
	float SUM = 0.0;
	for (int z = 0; z < m_ZCount; z++)
	{
		SUM += pow((Vx[z] * Wx[z]), 2.0f);
		SUM += pow((Vy[z] * Wy[z]), 2.0f);
		SUM += pow((Vz[z] * Wz[z]), 2.0f);
	}
//	printf("SVNORM SUM %e\n", SUM);
	return sqrt(SUM / (float)(m_ZCount * 3));
}

#endif


