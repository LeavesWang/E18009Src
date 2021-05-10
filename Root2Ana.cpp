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
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TH1S.h"
#include "TH1C.h"

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

const int AdcIChLow = 10, AdcIChUp = 4080;
// const double CalICh[16][2] = {{0, 1.33}, {0, 3.36}, {0, 1}, {0, 1.36}, {0, 1.4}, {0, 1.36}, {0, 1.32}, {0, 1.37}, {0, 1.3}, {0, 1.15}, {0, 1.45}, {0, 3.47}, {0, 1.32}, {0, 1.37}, {0, 3.11}, {0, 1.34}}; //need to check ch2

const int AdcTofLow = 500, AdcTofUp = 7680;
const int IdAdcClkTof[8] = {8, 10, 12, 14, 0, 2, 4, 6};																				//PMT# 1, 2, 3, 4, 5, 6, 7, 8
const double CalAdcClk[8] = {6.438028734, 6.63393242, 6.664039922, 6.447331894, 6.456309852, 6.65084734, 6.564673297, 6.494280891}; //unit: ps/ch
const int IdAdcTof[4] = {16, 18, 20, 24};																							//PMT# 5&1, 6&2, 7&3, 8&4
const double CalAdc[4] = {13.29158516, 13.06948277, 1, 1};																			//unit: ps/ch  //need to calibrate ADC-20 and ADC-24

const double CalTdc = 3.90625;						 //ps/ch
const int IdTdcTof[8] = {9, 11, 13, 15, 1, 3, 5, 7}; //PMT# 1, 2, 3, 4, 5, 6, 7, 8
const int TdcTofLow[8] = {41000, 42000, 43000, 37000, 18000, 18000, 18000, 18000};
const int TdcTofUp[8] = {61000, 62000, 63000, 57000, 28000, 28000, 28000, 28000};

const int IdQdcTof[8] = {1, 3, 5, 6, 8, 10, 12, 15};			   //PMT# 1, 2, 3, 4, 5, 6, 7, 8
const int QdcTofLow[8] = {100, 100, 100, 100, 100, 100, 100, 100}; //pedestal values for PMT# 1, 2, 3, 4, 5, 6, 7, 8

const int IdQdcMcp[2][4] = {{11, 15, 9, 13}, {3, 7, 1, 5}};			  // [2]: [0]High gain, [1]Low gain; [4]: [0]TL, [1]TR, [2]BR, [3]BL;
const int QdcMcpLow[2][4] = {{120, 180, 300, 160}, {20, 20, 20, 29}}; //pedestal values
const int QdcMcpUp = 3840;
const double ParXMcpCorrTof = 0.04; // unit ns/mm, xMCP-correction parameter for ToF

// const double CalXMcp[2][10] = {{0.733129, 26.744387, -0.091781, 1.043661, 0.047598, 9.192684, 2.637526, -0.929438, 2.056948, 0.576781}, {0.802060, 26.063777, -0.897100, 1.296354, 1.163047, 11.688516, 3.208674, -1.230582, -2.736673, 3.004569}};	   //[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*x*x*x+[6]*y*y*y+[7]*x*y+[8]*x*x*y+[9]*x*y*y  //two gain settings: [0] high gain setting; [1] comb-gain setting
// const double CalYMcp[2][10] = {{3.652901, 19.180574, 1.578795, -1.716251, 0.330541, 11.410052, -0.641449, -0.958885, 0.507911, 5.328422}, {3.727687, 18.762661, -0.510623, -1.588110, -0.511162, 10.227921, -1.138502, 0.227536, 0.858179, 4.114189}}; //[0]+[1]*y+[2]*y*y+[3]*x+[4]*x*x+[5]*y*y*y+[6]*x*x*x+[7]*x*y+[8]*y*y*x+[9]*x*x*y  //two gain settings: [0] high gain setting; [1] comb-gain setting
const double CalXMcp[2][10] = {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}}; //[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*x*x*x+[6]*y*y*y+[7]*x*y+[8]*x*x*y+[9]*x*y*y //for raw pos
const double CalYMcp[2][10] = {{0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}}; //[0]+[1]*y+[2]*y*y+[3]*x+[4]*x*x+[5]*y*y*y+[6]*x*x*x+[7]*x*y+[8]*y*y*x+[9]*x*x*y //for raw pos

const double CalTof[4][2] = {{500, 0.001}, {500, 0.001}, {500, 0.001}, {500, 0.001}}; //unit ns, ns/ps
const double Brho0 = 3.8806;														  //Tm  need to check
const double Disp = 106.8368;														  // mm/%  need to check
const double LOF = 60.763;															  //m  need to check
const double CalZ[4][2] = {{0, 1}, {0, 1}, {0, 1}, {0, 1}};

const string sSet[3] = {"Prod", "Ref", "Intm"};

void Root2Ana()
{
	if (access("../AnaData", F_OK) != 0)
		system("mkdir ../AnaData");

	int i = 0, j = 0, k = 0, m = 0, n = 0, p = 0;

	string sRoot, sAna;
	StrtMesytec madc, mtdc, cqdcTOF, mqdcMCP;
	StrtS800 s800;
	string setting = "";
	int run = 0;
	double tof[3] = {0}; //[3]: [0] TAC+ADC+clock; [1] TAC+ADC; [2]: regular CFD+TDC
	bool goodElc[3] = {0};
	double tD[3][8][8] = {0};
	double xPlaT[3][2] = {0}; //[2]: [0] for A1900 plastic from time info, [1] for S800 plastic from time info
	double yPlaT[3][2] = {0};
	double egyPMT[2][4] = {0}; //energies in 4 PMTs and each plastic
	double egyPla[2][2] = {0}; //first [2]: [0] is energy in plastic at A1900, [1] is energy in plastic at S800; second [2]: [0]arithmetic mean of four PMTs' Qs; [1]geometric mean of dour PMTs' Qs
	double xPlaQ[2] = {0};	   //[2]: [0] for A1900 plastic from amp info, [1] for S800 plastic from amp info
	double yPlaQ[2] = {0};
	double xMCP[2] = {0};		   //two gain settings: [0] high gain setting; [1] comb-gain setting
	double yMCP[2] = {0};		   //two gain settings: [0] high gain setting; [1] comb-gain setting
	double tofCorr[3] = {0};	   //Momentum-corrected ToF; [3]: [0] TAC+ADC+clock; [1] TAC+ADC; [2]: regular CFD+TDC
	int nGoodICh = 0;			   //count of good channels of IC
	bool goodICh[16] = {0};		   //if echa channel works well
	double delE0 = 0;			   // Average Delta E of 16 channels
	double delEICh[16] = {0};	   // Delta E of each channel
	bool goodCRDC = 0;			   //if CRDC works well
	double xCrdc[2][2] = {0};	   //first [2]: two CRDCs; second [2]: [0] is from gravity center; [1] is from Gaussian fit
	double xEgyCrdc[2][256] = {0}; //256 energies for 256 cathode pads
	double yCrdc[2] = {0};		   //y is only considered from electron's drift time
	double brho[2] = {0};
	double Z[3] = {0};
	double AoQ[3][2] = {0};

	long long iEnt = 0, nEnt = 0;

	double calICh[16][2] = {0};
	string sDraw = "", sCut = "";

	int nQdcTof[2] = {0};

	bool isCombGain = false, isCRDC = false;
	int nGoodMcp[2] = {0};
	double valQdcMcp[2][4] = {0};
	double xMcpRaw = 0, yMcpRaw = 0;

	int nGoodElc = 0, iElc = 0, numTime[3][2] = {0};
	double tPMT[8] = {0}, timeDet[2] = {0};

	int runMin, runMax, runNum;
	cout << "Input minimum and maximum numbers of run: ";
	cin >> runMin >> runMax;

	ostringstream ssRun;
	for (runNum = runMin; runNum <= runMax; runNum++)
	{
		ssRun.str("");
		ssRun << setw(4) << setfill('0') << runNum;
		sRoot = "../RootData/run-" + ssRun.str() + "-00.root";
		sAna = "../AnaData/ana-run-" + to_string(runNum) + ".root";
		printf("\n**********Now converting %s to %s!**********\n\n", sRoot.c_str(), sAna.c_str());

		TFile *fAna = new TFile(sAna.c_str(), "RECREATE");
		TTree *tAna = new TTree("tAna", "tree for data analysis");

		tAna->Branch("setting", &setting);
		tAna->Branch("run", &run, "run/I");
		tAna->Branch("tof", tof, "tof[3]/D");
		tAna->Branch("goodElc", goodElc, "goodElc[3]/O");
		tAna->Branch("tD", tD, "tD[3][8][8]/D");
		tAna->Branch("xPlaT", xPlaT, "xPlaT[3][2]/D");
		tAna->Branch("yPlaT", yPlaT, "yPlaT[3][2]/D");
		tAna->Branch("egyPMT", egyPMT, "egyPMT[2][4]/D");
		tAna->Branch("egyPla", egyPla, "egyPla[2][2]/D");
		tAna->Branch("xPlaQ", xPlaQ, "xPlaQ[2]/D");
		tAna->Branch("yPlaQ", yPlaQ, "yPlaQ[2]/D");
		tAna->Branch("xMCP", xMCP, "xMCP[2]/D");
		tAna->Branch("yMCP", yMCP, "yMCP[2]/D");
		tAna->Branch("tofCorr", tof, "tofCorr[3]/D");
		tAna->Branch("nGoodICh", &nGoodICh, "nGoodICh/I");
		tAna->Branch("goodICh", goodICh, "goodICh[16]/O");
		tAna->Branch("delE0", &delE0, "delE0/D");
		tAna->Branch("delEICh", delEICh, "delEICh[16]/D");
		tAna->Branch("goodCRDC", &goodCRDC, "goodCRDC/O");
		tAna->Branch("xCrdc", xCrdc, "xCrdc[2][2]/D");
		tAna->Branch("xEgyCrdc", xEgyCrdc, "xEgyCrdc[2][256]/D");
		tAna->Branch("yCrdc", yCrdc, "yCrdc[2]/D");
		tAna->Branch("brho", brho, "brho[2]/D");
		tAna->Branch("Z", Z, "Z[3]/D");
		tAna->Branch("AoQ", AoQ, "AoQ[3][2]/D");

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
		memset(&madc, 0, sizeof(madc));
		memset(&mtdc, 0, sizeof(mtdc));
		memset(&cqdcTOF, 0, sizeof(cqdcTOF));
		memset(&mqdcMCP, 0, sizeof(mqdcMCP));
		memset(&s800, 0, sizeof(s800));

		tData->SetBranchAddress("madc", &madc);
		tData->SetBranchAddress("mtdc", &mtdc);
		tData->SetBranchAddress("cqdcTOF", &cqdcTOF);
		tData->SetBranchAddress("mqdcMCP", &mqdcMCP);
		tData->SetBranchAddress("s800", &s800);

		tData->SetEstimate(-1);
		nEnt = tData->GetEntries();
		sDraw = "s800.ionCham[8]";
		sCut = "s800.ionCham[8]>" + to_string(AdcIChLow) + " && s800.ionCham[8]<" + to_string(AdcIChUp);
		for (j = 0; j < 16; j++)
			sDraw += ":s800.ionCham[" + to_string(j) + "]";
		int nData = tData->Draw(sDraw.c_str(), sCut.c_str(), "goff");
		TGraph *grIC[16];
		double *valIC[16];
		double *valIC8 = tData->GetV1();
		for (j = 0; j < 16; j++)
		{
			valIC[j] = tData->GetVal(j + 1);
			grIC[j] = new TGraph();
			n = 0;
			for (i = 0; i < nData; i++)
				if (valIC[j][i] > AdcIChLow && valIC[j][i] < AdcIChUp)
					grIC[j]->SetPoint(n++, valIC[j][i], valIC8[i]);
			if (n < 2)
				continue;
			TFitResultPtr fitResIC = grIC[j]->Fit("pol1", "SQ");
			int fitStIC = fitResIC;
			if (fitStIC != 0 && fitStIC != 4000)
				continue;
			calICh[j][0] = fitResIC->Parameter(0);
			calICh[j][1] = fitResIC->Parameter(1);
		}
		for (j = 0; j < 16; j++)
			if (!grIC[j])
				delete grIC[j];

		double mcpGainMat[8][2];
		for (i = 0; i < 8; i++)
		{
			mcpGainMat[i][0] = 0;
			mcpGainMat[i][1] = 1;
		}
		if (isCombGain)
		{
			double *errCh = new double[nEnt];
			for (iEnt = 0; iEnt < nEnt; iEnt++)
				errCh[iEnt] = 0.5;
			// string strCanv="gainMat_"+to_string(runNum);
			// TCanvas *cGain=new TCanvas(strCanv.c_str(), strCanv.c_str());
			// cGain->Divide(2,2);
			TGraphErrors *gr[4];
			for (i = 0; i < 4; i++)
			{
				j = IdQdcMcp[0][i];
				m = IdQdcMcp[1][i];
				// cGain->cd(i+1);
				string sCut = "mqdcMCP.data[" + to_string(j) + "]>" + to_string(QdcMcpLow[0][i]) + "&&mqdcMCP.data[" + to_string(j) + "]<3840&&mqdcMCP.data[" + to_string(m) + "]>" + to_string(QdcMcpLow[1][i]) + "&&mqdcMCP.data[" + to_string(m) + "]<3840";
				string sDraw = "mqdcMCP.data[" + to_string(j) + "]-" + to_string(QdcMcpLow[0][i]) + ":mqdcMCP.data[" + to_string(m) + "]-" + to_string(QdcMcpLow[1][i]);

				long long nData = tData->Draw(sDraw.c_str(), sCut.c_str(), "goff");
				if (nData < 2)
					continue;
				double *highGain = tData->GetV1();
				double *lowGain = tData->GetV2();
				gr[i] = new TGraphErrors(nData, lowGain, highGain, errCh, errCh);
				// gr[i]->Draw("AP");
				// gr[i]->SetTitle(("high"+to_string(j)+"_vs_low"+to_string(m)).c_str());
				TFitResultPtr fitRes = gr[i]->Fit("pol1", "SQ");
				int fitSt = fitRes;
				if (fitSt != 0 && fitSt != 4000)
					continue;
				mcpGainMat[j][0] = fitRes->Parameter(0);
				mcpGainMat[j][1] = fitRes->Parameter(1);
				// printf("%f %f\n",mcpGainMat[j][0],mcpGainMat[j][1]);
			}
			delete[] errCh;
			// cGain->SaveAs(("/home/kailong/ExpData/Jul2018/Graphs/Charts/"+strCanv+".png").c_str());
			// cGain->Close();
			for (i = 0; i < 4; i++)
				if (!gr[i])
					delete gr[i];
			// delete cGain;
		}

		for (iEnt = 0; iEnt < nEnt; iEnt++)
		// for (iEnt = 0; iEnt < 100; iEnt++)
		{
			tData->GetEntry(iEnt);
			if (s800.trig == 1 || s800.trig == 16)
			{
				memset(tof, 0, sizeof(tof));
				memset(goodElc, 0, sizeof(goodElc));
				memset(tD, 0, sizeof(tD));
				memset(xPlaT, 0, sizeof(xPlaT));
				memset(yPlaT, 0, sizeof(yPlaT));
				memset(egyPMT, 0, sizeof(egyPMT));
				memset(egyPla, 0, sizeof(egyPla));
				memset(xPlaQ, 0, sizeof(xPlaQ));
				memset(yPlaQ, 0, sizeof(yPlaQ));
				memset(xMCP, 0, sizeof(xMCP));
				memset(yMCP, 0, sizeof(yMCP));
				memset(tofCorr, 0, sizeof(tofCorr));
				nGoodICh = 0;
				memset(goodICh, 0, sizeof(goodICh));
				delE0 = 0;
				memset(delEICh, 0, sizeof(delEICh));
				goodCRDC = false;
				memset(xCrdc, 0, sizeof(xCrdc));
				memset(xEgyCrdc, 0, sizeof(xEgyCrdc));
				memset(yCrdc, 0, sizeof(yCrdc));
				memset(brho, 0, sizeof(brho));
				memset(Z, 0, sizeof(Z));
				memset(AoQ, 0, sizeof(AoQ));

				TRandom3 r(0);

				for (i = 0; i < 16; i++)
					if (s800.ionCham[i] > AdcIChLow && s800.ionCham[i] < AdcIChUp)
					{
						goodICh[i] = true;
						nGoodICh++;

						delEICh[i] = calICh[i][0] + calICh[i][1] * (s800.ionCham[i] + r.Uniform(-0.5, 0.5));
						delE0 += delEICh[i];
					}
				if (nGoodICh > 0)
					delE0 /= nGoodICh;

				if (isCRDC)
					for (m = 0; m < 2; m++)
						if (s800.crdcAnode[m][0] > 0 && s800.crdcAnode[m][1] > 0)
						{
							//calculate y from the drift time to anode
							yCrdc[m] = s800.crdcAnode[m][1];

							TH1S hx("hx", "hx", 256, 0, 256);
							TH1C hc("hc", "hc", 256, 0, 256);
							for (i = 1; i <= 4; i++)
								for (j = 0; j < 64; j++)
								{
									k = (i - 1) * 64 + j;
									if (k > 90 && k < 130 && s800.crdcCath[m][i][j] > 100 && s800.crdcCath[m][i][j] < 990)
									{
										xEgyCrdc[m][k] = s800.crdcCath[m][i][j];
										hx.Fill(k, s800.crdcCath[m][i][j]);
										hc.Fill(k);
									}
								}
							if (hx.GetEntries() > 2)
							{
								double maxPad = hx.GetBinCenter(hx.GetMaximumBin());
								hx.SetAxisRange(maxPad - 10, maxPad + 10, "X");
								hc.SetAxisRange(maxPad - 10, maxPad + 10, "X");
								if (hc.Integral() > 2)
								{
									//calculate x from gravity center
									xCrdc[m][0] = hx.GetMean();
									//calculate x from Gaussian fit
									TFitResultPtr fitX = hx.Fit("gaus", "S0Q");
									int fitSt = fitX;
									if (fitSt != 0 && fitX->Parameter(1) > 0 && fitX->Parameter(1) < 256)
										xCrdc[m][1] = fitX->Parameter(1);

									goodCRDC = true;
								}
							}
						}

				memset(nGoodMcp, 0, sizeof(nGoodMcp));
				memset(valQdcMcp, 0, sizeof(valQdcMcp));
				for (i = 0; i < 4; i++)
				{
					j = IdQdcMcp[0][i];
					if (mqdcMCP.data[j] > QdcMcpLow[0][i] && mqdcMCP.data[j] < QdcMcpUp)
					{
						nGoodMcp[0]++;
						valQdcMcp[0][i] = mqdcMCP.data[j] - QdcMcpLow[0][i] + r.Uniform(-0.5, 0.5);

						nGoodMcp[1]++;
						valQdcMcp[1][i] = valQdcMcp[0][i];
					}
					m = IdQdcMcp[1][i];
					if (mqdcMCP.data[j] >= QdcMcpUp && mqdcMCP.data[m] > QdcMcpLow[1][i] && mqdcMCP.data[m] < QdcMcpUp)
					{

						valQdcMcp[1][i] = mcpGainMat[m][0] + mcpGainMat[m][1] * (mqdcMCP.data[m] - QdcMcpLow[1][i] + r.Uniform(-0.5, 0.5));
						if (valQdcMcp[1][i] > QdcMcpUp - QdcMcpLow[0][i])
							nGoodMcp[1]++;
					}
				}
				for (i = 0; i < 2; i++)
					if (nGoodMcp[i] == 4)
					{
						xMcpRaw = (valQdcMcp[i][1] + valQdcMcp[i][2] - valQdcMcp[i][0] - valQdcMcp[i][3]) / (valQdcMcp[i][0] + valQdcMcp[i][1] + valQdcMcp[i][2] + valQdcMcp[i][3]);

						yMcpRaw = (valQdcMcp[i][0] + valQdcMcp[i][1] - valQdcMcp[i][2] - valQdcMcp[i][3]) / (valQdcMcp[i][0] + valQdcMcp[i][1] + valQdcMcp[i][2] + valQdcMcp[i][3]);

						xMCP[i] = CalXMcp[i][0] + CalXMcp[i][1] * xMcpRaw + CalXMcp[i][2] * pow(xMcpRaw, 2) + CalXMcp[i][3] * yMcpRaw + CalXMcp[i][4] * pow(yMcpRaw, 2) + CalXMcp[i][5] * pow(xMcpRaw, 3) + CalXMcp[i][6] * pow(yMcpRaw, 3) + CalXMcp[i][7] * xMcpRaw * yMcpRaw + CalXMcp[i][8] * pow(xMcpRaw, 2) * yMcpRaw + CalXMcp[i][9] * xMcpRaw * pow(yMcpRaw, 2);

						yMCP[i] = CalYMcp[i][0] + CalYMcp[i][1] * yMcpRaw + CalYMcp[i][2] * pow(yMcpRaw, 2) + CalYMcp[i][3] * xMcpRaw + CalYMcp[i][4] * pow(xMcpRaw, 2) + CalYMcp[i][5] * pow(yMcpRaw, 3) + CalYMcp[i][6] * pow(xMcpRaw, 3) + CalYMcp[i][7] * yMcpRaw * xMcpRaw + CalYMcp[i][8] * pow(yMcpRaw, 2) * xMcpRaw + CalYMcp[i][9] * yMcpRaw * pow(xMcpRaw, 2);

						brho[i] = Brho0 * (1 + xMCP[i] / Disp / 100);
					}

				memset(nQdcTof, 0, sizeof(nQdcTof));
				for (i = 0; i < 8; i++)
				{
					k = IdQdcTof[i];
					if (cqdcTOF.data[k] > QdcTofLow[i])
					{
						j = i % 4;
						n = i / 4;
						nQdcTof[n]++;
						egyPMT[n][j] = cqdcTOF.data[k] + r.Uniform(-0.5, 0.5) - QdcTofLow[i];
					}
				}
				if (nQdcTof[0] == 4)
				{
					xPlaQ[0] = log(egyPMT[0][1] / egyPMT[0][3]); // need to confirm
					yPlaQ[0] = log(egyPMT[0][2] / egyPMT[0][0]);
					egyPla[0][0] = (egyPMT[0][0] + egyPMT[0][1] + egyPMT[0][2] + egyPMT[0][3]) / 4;
					egyPla[0][1] = pow(egyPMT[0][0] * egyPMT[0][1] * egyPMT[0][2] * egyPMT[0][3], 0.25);
				}
				if (nQdcTof[1] == 4)
				{
					xPlaQ[1] = log(egyPMT[1][1] / egyPMT[1][3]); // need to confirm
					yPlaQ[1] = log(egyPMT[1][0] / egyPMT[1][2]);
					egyPla[1][0] = (egyPMT[1][0] + egyPMT[1][1] + egyPMT[1][2] + egyPMT[1][3]) / 4;
					egyPla[1][1] = pow(egyPMT[1][0] * egyPMT[1][1] * egyPMT[1][2] * egyPMT[1][3], 0.25);
				}

				if (nGoodICh > 8 && nQdcTof[0] == 4 && nQdcTof[1] == 4 && nGoodMcp[0] == 4)
				{
					nGoodElc = 0;
					memset(numTime, 0, sizeof(numTime));
					for (iElc = 0; iElc < 3; iElc++)
					{
						memset(tPMT, 0, sizeof(tPMT));
						memset(timeDet, 0, sizeof(timeDet));
						if (iElc != 1) //For TAC+ADC+clock or regularCFD+TDC
						{
							for (i = 0; i < 8; i++)
							{
								k = i / 4;
								if (iElc == 0) //For TAC+ADC+clock
								{
									j = IdAdcClkTof[i];
									if (madc.data[j] > AdcTofLow && madc.data[j] < AdcTofUp)
									{
										numTime[0][k]++;
										tPMT[i] = (CalAdcClk[i] + 0) * (madc.data[j] + r.Uniform(-0.5, 0.5));
									}
								}
								if (iElc == 2) //For regularCFD+TDC
								{
									j = IdTdcTof[i];
									if (mtdc.data[j] > TdcTofLow[i] && mtdc.data[j] < TdcTofUp[i])
									{
										numTime[2][k]++;
										tPMT[i] = (CalTdc + 0) * (mtdc.data[j] + r.Uniform(-0.5, 0.5));
									}
								}
								timeDet[k] += tPMT[i];
							}
							for (i = 0; i < 7; i++)
								for (j = i + 1; j < 8; j++)
									if (abs(tPMT[i]) > 0 && abs(tPMT[j]) > 0)
									{
										tD[iElc][i][j] = tPMT[i] - tPMT[j];
										tD[iElc][j][i] = -tD[iElc][i][j];
									}
							if (numTime[iElc][0] == 4 && numTime[iElc][1] == 4)
							{
								goodElc[iElc] = true;
								nGoodElc++;
								tof[iElc] = timeDet[1] / 4 - timeDet[0] / 4;

								xPlaT[iElc][0] = tD[iElc][3][1];
								yPlaT[iElc][0] = tD[iElc][0][2];

								xPlaT[iElc][1] = tD[iElc][7][5];
								yPlaT[iElc][1] = tD[iElc][6][4];
							}
						}
						if (iElc == 1) //For TAC+ADC
						{
							for (i = 0; i < 4; i++)
							{
								j = IdAdcTof[i];
								m = i;
								p = i + 4;
								if (madc.data[j] > AdcTofLow && madc.data[j] < AdcTofUp)
								{
									numTime[1][0]++;
									numTime[1][1]++;
									tD[1][m][p] = (CalAdc[i] + 0) * (madc.data[j] + r.Uniform(-0.5, 0.5));
									tD[1][p][m] = -tD[1][p][m];
									tof[1] += tD[1][p][m];
								}
							}
							if (numTime[1][0] == 4)
							{
								tof[1] /= 4;
								goodElc[iElc] = true;
								nGoodElc++;
							}
						}
					}

					if (nGoodElc > 1)
					{
						for (iElc = 0; iElc < 3; iElc++)
							if (goodElc[iElc])
							{
								tof[iElc] = CalTof[iElc][0] + CalTof[iElc][1] * tof[iElc];

								tofCorr[iElc] = tof[iElc] - ParXMcpCorrTof * tof[iElc];

								double beta = LOF / tof[iElc] / 0.299792458;
								if (beta > 0 && beta < 1)
								{
									double gamma = 1 / sqrt(1 - beta * beta);
									Z[iElc] = sqrt(delE0 / (log(5930 / (1 / beta / beta - 1)) / beta / beta - 1));
									Z[iElc] = CalZ[iElc][0] + CalZ[iElc][1] * Z[iElc];

									for (k = 0; k < 2; k++)
										AoQ[iElc][k] = brho[k] / beta / gamma * 0.32184043;
								}
							}
						tAna->Fill();
					}
				}
			}
		}
		fRoot->Close();
		delete fRoot;

		fAna->cd();
		tAna->Write();
		fAna->Close();
		delete fAna;
	}
}

#ifndef __CINT__
/*
void StandaloneApplication(int argc, char** argv)
{
	Root2Ana();
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
	Root2Ana();
	return 0;
}
#endif