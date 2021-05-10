#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <unistd.h>
#include <cmath>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

struct StrtMesytec
{
	int modID;
	int data[32];
	int modRes;
	int modEC_TS;
};

struct StrtS800
{
	int tS;
	int eC;
	int trig;				// =1: Coincidence (usual case, S800 scintillator provides); =16: Secondary (our ToF or MCP provides)
	int tof[16];			//S800 ToF packet
	int scin[2][2];			//S800 scintillator: first [2]: two PMTs, second [2]: [0] energy; [1] time
	int ionCham[16];		//S800 ion chamber; [16] means 16 segaments
	int crdcCath[2][5][64]; //[2]: two CRDC; [5]: [0]---sample; [1]-->[4]---energy
	int crdcAnode[2][2];	// [0]: energy; [1]: time
	int mesyTDC[16];		//[0]: time from S800; [2]: time from A1900
};

void Root2Count()
{
	if (access("../AnaData/Count", F_OK) != 0)
		system("mkdir ../AnaData/Count");

	int i = 0, j = 0, k = 0, m = 0, n = 0, p = 0;

	const int IdTdcTof[8] = {9, 11, 13, 15, 1, 3, 5, 7};				  //PMT# 1, 2, 3, 4, 5, 6, 7, 8
	const int IdQdcMcp[2][4] = {{11, 15, 9, 13}, {3, 7, 1, 5}};			  // [2]: [0]High gain, [1]Low gain; [4]: [0]TL, [1]TR, [2]BR, [3]BL;
	const int QdcMcpLow[2][4] = {{100, 100, 100, 100}, {20, 20, 20, 20}}; //pedestal values  //***need to check!***

	const string sDet[6] = {"A1900Pla", "MCP", "CRDC1", "S800Pla", "CRDC2", "IonCham"};

	const string sSet[3] = {"Prod", "Ref", "Intm"};

	string sCoinDet[64];
	sCoinDet[0] = "ExpDefault";
	int iCoin = 0;
	for (i = 0; i < 6; i++)
		sCoinDet[++iCoin] = sDet[i];
	for (i = 0; i < 6; i++)
		for (j = i + 1; j < 6; j++)
			sCoinDet[++iCoin] = sDet[i] + "_" + sDet[j];
	for (i = 0; i < 6; i++)
		for (j = i + 1; j < 6; j++)
			for (k = j + 1; k < 6; k++)
				sCoinDet[++iCoin] = sDet[i] + "_" + sDet[j] + "_" + sDet[k];
	for (i = 0; i < 6; i++)
		for (j = i + 1; j < 6; j++)
			for (k = j + 1; k < 6; k++)
				for (m = k + 1; m < 6; m++)
					sCoinDet[++iCoin] = sDet[i] + "_" + sDet[j] + "_" + sDet[k] + "_" + sDet[m];
	for (i = 0; i < 6; i++)
		for (j = i + 1; j < 6; j++)
			for (k = j + 1; k < 6; k++)
				for (m = k + 1; m < 6; m++)
					for (n = m + 1; n < 6; n++)
						sCoinDet[++iCoin] = sDet[i] + "_" + sDet[j] + "_" + sDet[k] + "_" + sDet[m] + "_" + sDet[n];
	sCoinDet[63] = sDet[0] + "_" + sDet[1] + "_" + sDet[2] + "_" + sDet[3] + "_" + sDet[4] + "_" + sDet[5];

	int runMin, runMax, runNum;
	cout << "Input minimum and maximum numbers of run: ";
	cin >> runMin >> runMax;

	string sCount = "../AnaData/Count/count-run-" + to_string(runMin) + "__" + to_string(runMax) + ".root";

	string sRoot = "";
	StrtMesytec mtdc, mqdcMCP;
	StrtS800 s800;
	string setting = "";
	int run = 0;
	int count[64] = {0};

	TFile *fCount = new TFile(sCount.c_str(), "RECREATE");
	TTree *tCount = new TTree("tCount", "tree for count analysis");
	tCount->Branch("setting", &setting);
	tCount->Branch("run", &run, "run/I");
	for (i = 0; i < 64; i++)
		tCount->Branch(sCoinDet[i].c_str(), &count[i], (sCoinDet[i] + "/I").c_str());

	int nDet[6] = {0};
	bool goodDet[6] = {0};

	long long nEnt = 0, iEnt = 0;

	ostringstream ssRun;
	string strRun;
	for (runNum = runMin; runNum <= runMax; runNum++)
	{
		ssRun.str("");
		ssRun << setw(4) << setfill('0') << runNum;
		strRun = ssRun.str();
		sRoot = "../RootData/run-" + strRun + "-00.root";

		printf("\n**********Now converting %s to %s!**********\n\n", sRoot.c_str(), sCount.c_str());

		TFile *fRoot = new TFile(sRoot.c_str());
		if (fRoot->IsZombie())
		{
			cout << "Error in opening " << sRoot << "!\n";
			continue;
		}
		TTree *tData;
		fRoot->GetObject("tData", tData);
		if (!tData)
		{
			cout << "Error read the tree of tData!\n";
			continue;
		}
		memset(&mtdc, 0, sizeof(mtdc));
		memset(&mqdcMCP, 0, sizeof(mqdcMCP));
		memset(&s800, 0, sizeof(s800));

		tData->SetBranchAddress("mtdc", &mtdc);
		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);
		tData->SetBranchAddress("s800", &s800);

		for (i = 0; i < 3; i++)
		{
			string sfSet = "./RunNum_" + sSet[i] + ".dat";
			ifstream fSet(sfSet.c_str());
			string sRead;
			bool isFound = false;
			while (getline(fSet, sRead))
				if (sRead.find(to_string(runNum)) != string::npos)
				{
					isFound = true;
					break;
				}
			if (isFound)
			{
				setting = sSet[i];
				break;
			}
			if (i == 2 && !isFound)
				setting = "other";
		}
		run = runNum;

		memset(count, 0, sizeof(count));

		nEnt = tData->GetEntries();
		// nEnt = 100;
		for (iEnt = 0; iEnt < nEnt; iEnt++)
		{
			tData->GetEntry(iEnt);

			memset(nDet, 0, sizeof(nDet));
			memset(goodDet, 0, sizeof(goodDet));

			for (i = 0; i < 4; i++)
				if (mtdc.data[IdTdcTof[i]] > 1 && mtdc.data[IdTdcTof[i]] < 65535)
					nDet[0]++;
			if (nDet[0] == 4)
				goodDet[0] = true;

			for (i = 0; i < 4; i++)
			{
				j = IdQdcMcp[0][i];
				k = IdQdcMcp[1][i];
				if ((mqdcMCP.data[j] > QdcMcpLow[0][i] && mqdcMCP.data[j] < 3840) || (mqdcMCP.data[k] > QdcMcpLow[1][i] && mqdcMCP.data[k] < 3840))
					nDet[1]++;
			}
			if (nDet[1] == 4)
				goodDet[1] = true;

			for (i = 4; i < 8; i++)
				if (mtdc.data[IdTdcTof[i]] > 1 && mtdc.data[IdTdcTof[i]] < 65535)
					nDet[3]++;
			if (nDet[3] == 4)
				goodDet[3] = true;

			for (i = 0; i < 2; i++)
			{
				p = 2 + 2 * i;
				if (s800.crdcAnode[i][1] > 10 && s800.crdcAnode[i][1] < 4080)
				{
					for (j = 1; j <= 4; j++)
						for (k = 0; k < 64; k++)
						{
							m = (j - 1) * 64 + k;
							if (s800.crdcCath[i][j][k] > 100 && s800.crdcCath[i][j][k] < 990)
								n++;
						}
					if (n > 2)
						nDet[p]++;
				}
				if (nDet[p] == 1)
					goodDet[p] = true;
			}

			for (i = 0; i < 16; i++)
				if (s800.ionCham[i] > 10 && s800.ionCham[i] < 4080)
					nDet[5]++;
			if (nDet[5] > 14)
				goodDet[5] = true;

			count[0]++;
			iCoin = 0;
			for (i = 0; i < 6; i++)
				if (goodDet[i])
					count[++iCoin]++;
			for (i = 0; i < 6; i++)
				for (j = i + 1; j < 6; j++)
					if (goodDet[i] && goodDet[j])
						count[++iCoin]++;
			for (i = 0; i < 6; i++)
				for (j = i + 1; j < 6; j++)
					for (k = j + 1; k < 6; k++)
						if (goodDet[i] && goodDet[j] && goodDet[k])
							count[++iCoin]++;
			for (i = 0; i < 6; i++)
				for (j = i + 1; j < 6; j++)
					for (k = j + 1; k < 6; k++)
						for (m = k + 1; m < 6; m++)
							if (goodDet[i] && goodDet[j] && goodDet[k] && goodDet[m])
								count[++iCoin]++;
			for (i = 0; i < 6; i++)
				for (j = i + 1; j < 6; j++)
					for (k = j + 1; k < 6; k++)
						for (m = k + 1; m < 6; m++)
							for (n = m + 1; n < 6; n++)
								if (goodDet[i] && goodDet[j] && goodDet[k] && goodDet[m] && goodDet[n])
									count[++iCoin]++;
			count[63]++;
		}
		tCount->Fill();
		fRoot->Close();
		delete fRoot;
	}
	fCount->cd();
	tCount->Write();
	delete tCount;
	fCount->Close();
	delete fCount;
}

#ifndef __CINT__
/*
void StandaloneApplication(int argc, char** argv)
{
	Root2Count();
}

int main(int argc, char** argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	app.Run();
	return 0;
}
*/
int main()
{
	Root2Count();
	return 0;
}
#endif