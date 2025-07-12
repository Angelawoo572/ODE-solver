#pragma once

#include "CHAMR.h"
#include "Summary.h"

#include <filesystem>
#include <chrono>


// 
// 異常時のメッセージは↓で囲むようにする
// printf("**********************************************************************\n");
// 


///////////////////////////////////////////////////////////////////////////////////////////////
// 
int main(int argc, char* argv[])
{
	printf("**********************************************************************\n");
	for (int i = 0; i < argc; i++)
		printf("arg%d %s\n", i, argv[i]);
	printf("**********************************************************************\n");
	printf("\n");
	printf("\n");

	// 入力チェック
	if (3 > argc)
	{
		printf("**********************************************************************\n");
		printf("Too few parameters.\n");
		printf("**********************************************************************\n");
		return 0;
	}

	std::chrono::system_clock::time_point time = std::chrono::system_clock::now();
	std::time_t timeStart = std::chrono::system_clock::to_time_t(time);

	const std::filesystem::path strPath = argv[1];
	if (false == std::filesystem::exists(strPath))
	{
		printf("**********************************************************************\n");
		printf("Folder not found. %s\n", argv[1]);
		printf("**********************************************************************\n");
		return 0;
	}

	const std::string strSimonPath = strPath.string() + "/" + FOLDER_SIMON;
	const std::string strResultPath = strPath.string() + "/" + FOLDER_RESULT;

	printf("%s\n", strSimonPath.c_str());
	printf("%s\n", strResultPath.c_str());

	// 結果の保存先を作成
	std::filesystem::create_directory(strSimonPath);
	std::filesystem::create_directory(strResultPath);

	// 始めは計算と集計を同時に遣っていいたが2024/03/21以降は分けて実行するよう変更。
	// この修正はヘッドファイル作成に時間がかかるため、計算中にヘッドファイルを作る時短が目的。
	// そのため、計算・集計を使うには2回動かす必要がある。

	// 計算
	// 入力パラメータが5個以上だと集計だけとなります。
	if (5 > argc)
	{
		// breakを使いたかったので。
		do
		{
			int nLogGrain = 0;
#ifdef GPU_LOG
			if (3 < argc)
				nLogGrain = strtol(argv[3], NULL, 10);
#endif

			// 条件ファイル検索 (*.par)
			std::string strParameterPath("");
			for (const std::filesystem::directory_entry &x1: std::filesystem::recursive_directory_iterator(strPath))
			{
				std::string str = x1.path().filename().string();
//				std::string str = x1.path().extension().string();

				if (0 != str.compare(FILE_PARAMETER))
					continue;

				strParameterPath = str;
				break;
			}
			if (false != strParameterPath.empty())
			{
				printf("**********************************************************************\n");
				printf("Parameter file not found.\n");
				printf("**********************************************************************\n");
				break;
			}

			// 計算
			CHAMR hamr;
			if (false == hamr.Run(strPath.string().c_str(), argv[2], strParameterPath.c_str(), nLogGrain))
			{
				printf("**********************************************************************\n");
				printf("There was an anomaly in the calculation.\n");
				printf("**********************************************************************\n");
			}
		}
		while (0);

		std::chrono::duration<double> timeCalculation = std::chrono::system_clock::now() - time;

		printf("\n\n");
		printf("////////////////////////////////////////////////////\n");
		printf("\tCalculation\t%lf sec\n", timeCalculation.count());
		printf("\t\t\t\t%lf min\n", (timeCalculation.count() /60.0f));
		printf("////////////////////////////////////////////////////\n");
		printf("\n\n");
	}
	else
	{
		// 条件ファイル検索 (*.par)
		std::string strParameterFile("");
		for (const std::filesystem::directory_entry &x1: std::filesystem::recursive_directory_iterator(strPath))
		{
			std::string str = x1.path().extension().string();

			if (0 != str.compare(EXTENSION_PAR))
				continue;

			strParameterFile = strPath.string() + "/" + x1.path().filename().string();
			break;
		}

		// ヘッドファイル
		std::string strReaderPath = strPath.string() + "/" + FILE_HEAD;
		if (4 < argc)
			strReaderPath = strPath.string() + "/" + argv[4];

		if (false == std::filesystem::exists(strReaderPath))
		{
			printf("**********************************************************************\n");
			printf("Reader head not found. %s\n", strReaderPath.c_str());
			printf("**********************************************************************\n");
		}
		else
		{
			// endファイル
			std::vector<std::string> vecEndFile;
			for (const std::filesystem::directory_entry &x1: std::filesystem::recursive_directory_iterator(strSimonPath))
			{
				std::string str = x1.path().filename().string();

				if (std::string::npos == str.find(EXTENSION_END))
					continue;

				if (std::string::npos == str.find(argv[2]))
					continue;

				vecEndFile.push_back(x1.path().string());
			}

			for (const std::string &x: vecEndFile)
			{
				printf("%s\n", x.c_str());

				CSummary Summary;
				Summary.SubMain(strSimonPath.c_str(), strResultPath.c_str(), strPath.string().c_str(), strReaderPath.c_str(), strParameterFile.c_str(), x.c_str(), argv[2]);
			}
		}

		std::chrono::duration<double> timeElaspsed = std::chrono::system_clock::now() - time;

		printf("\n\n");
		printf("////////////////////////////////////////////////////\n");
		printf("\tSummary\t%lf sec\n", timeElaspsed.count());
		printf("\t\t\t%lf min\n\n\n", (timeElaspsed.count() /60.0f));
		printf("////////////////////////////////////////////////////\n");
	}

#ifdef _DEBUG
	system("pause");
#endif
}

