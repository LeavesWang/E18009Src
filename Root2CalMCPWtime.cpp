#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <unistd.h>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TStyle.h"

using namespace std;

struct StrtMesytec
{
	int modID;
	int data[32];
	int modRes;
	int modEC_TS;
};

void Root2CalMCPWtime()
{
	const int QdcUp = 3840;
	int QdcLow[32] = {0}; //pedestal values
	QdcLow[9] = 300;
	QdcLow[11] = 120;
	QdcLow[13] = 160;
	QdcLow[15] = 180;
	double TdcLow[32] = {0}, TdcUp[32] = {0};
	TdcLow[1] = 16E3;
	TdcLow[3] = 16E3;
	TdcLow[5] = 16E3;
	TdcLow[7] = 16E3;
	TdcUp[1] = 30E3;
	TdcUp[3] = 30E3;
	TdcUp[5] = 30E3;
	TdcUp[7] = 30E3;

	string sRoot, sCalMcp;
	StrtMesytec mtdc, mqdcMCP;

	int i = 0, j = 0, k = 0;

	int run = 0;
	double ampTR = 0, ampTL = 0, ampBL = 0, ampBR = 0;
	double xRaw = 0, yRaw = 0, tdcMCP = 0, timeMCP2S800FP = 0;

	long long iEnt = 0, nEnt = 0;

	int runMin, runMax, runNum;
	cout << "Input minimum and maximum numbers of run: ";
	cin >> runMin >> runMax;

	if (access("../AnaData", F_OK) != 0)
		system("mkdir ../AnaData");
	sCalMcp = "../AnaData/calMCPwTime-run-" + to_string(runMin) + "--" + to_string(runMax) + ".root";
	TFile *fCalMcp = new TFile(sCalMcp.c_str(), "RECREATE");
	TTree *tCalMcp = new TTree("tCalMcpWTime", "tree for calibrating MCP with time");
	tCalMcp->Branch("run", &run, "run/I");
	tCalMcp->Branch("ampTR", &ampTR, "ampTR/D");
	tCalMcp->Branch("ampTL", &ampTL, "ampTL/D");
	tCalMcp->Branch("ampBR", &ampBR, "ampBR/D");
	tCalMcp->Branch("ampBL", &ampBL, "ampBL/D");
	tCalMcp->Branch("xRaw", &xRaw, "xRaw/D");
	tCalMcp->Branch("yRaw", &yRaw, "yRaw/D");
	tCalMcp->Branch("tdcMCP", &tdcMCP, "tdcMCP/D");
	tCalMcp->Branch("timeMCP2S800FP", &timeMCP2S800FP, "timeMCP2S800FP/D");

	ostringstream ssRun;
	for (runNum = runMin; runNum <= runMax; runNum++)
	{
		ssRun.str("");
		ssRun << setw(4) << setfill('0') << runNum;

		sRoot = "../RootData/run-" + ssRun.str() + "-00.root";
		printf("\n**********Now converting %s to %s!**********\n\n", sRoot.c_str(), sCalMcp.c_str());

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
		memset(&mqdcMCP, 0, sizeof(mqdcMCP));
		memset(&mtdc, 0, sizeof(mtdc));
		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);
		tData->SetBranchAddress("mtdc", &mtdc);

		nEnt = tData->GetEntries();
		for (iEnt = 0; iEnt < nEnt; iEnt++)
		{
			tData->GetEntry(iEnt);
			TRandom3 r(0);
			run = runNum;
			ampTR = 0;
			ampTL = 0;
			ampBL = 0;
			ampBR = 0;
			xRaw = 0;
			yRaw = 0;
			tdcMCP = 0;
			timeMCP2S800FP = 0;

			if (mqdcMCP.data[9] > QdcLow[9] && mqdcMCP.data[9] < QdcUp)
				ampBR = mqdcMCP.data[9] - QdcLow[9] + r.Uniform(-0.5, 0.5);
			if (mqdcMCP.data[11] > QdcLow[11] && mqdcMCP.data[11] < QdcUp)
				ampTL = mqdcMCP.data[11] - QdcLow[11] + r.Uniform(-0.5, 0.5);
			if (mqdcMCP.data[13] > QdcLow[13] && mqdcMCP.data[13] < QdcUp)
				ampBL = mqdcMCP.data[13] - QdcLow[13] + r.Uniform(-0.5, 0.5);
			if (mqdcMCP.data[15] > QdcLow[15] && mqdcMCP.data[15] < QdcUp)
				ampTR = mqdcMCP.data[15] - QdcLow[15] + r.Uniform(-0.5, 0.5);

			if (ampBR > 0 && ampTL > 0 && ampBL > 0 && ampTR > 0)
			{
				xRaw = (ampBR + ampTR - ampTL - ampBL) / (ampBR + ampTR + ampTL + ampBL);
				yRaw = (ampTL + ampTR - ampBR - ampBL) / (ampBR + ampTR + ampTL + ampBL);

				if (mtdc.data[14] > 1 && mtdc.data[14] < 65536)
				{
					tdcMCP = mtdc.data[14];
					if (mtdc.data[1] > TdcLow[1] && mtdc.data[1] < TdcUp[1] && mtdc.data[3] > TdcLow[3] && mtdc.data[3] < TdcUp[3] && mtdc.data[5] > TdcLow[5] && mtdc.data[5] < TdcUp[5] && mtdc.data[7] > TdcLow[7] && mtdc.data[7] < TdcUp[7])
						timeMCP2S800FP = 3.90625 * (mtdc.data[1] + r.Uniform(-0.5, 0.5) + mtdc.data[3] + r.Uniform(-0.5, 0.5) + mtdc.data[5] + r.Uniform(-0.5, 0.5) + mtdc.data[7] + r.Uniform(-0.5, 0.5)) / 4 - 3.90625 * (mtdc.data[14] + r.Uniform(-0.5, 0.5));
				}
				tCalMcp->Fill();
			}
		}
		fRoot->Close();
	}
	fCalMcp->cd();
	tCalMcp->Write();
	fCalMcp->Close();
}

#ifndef __CINT__
int main()
{
	Root2CalMCPWtime();
	return 0;
}
#endif
