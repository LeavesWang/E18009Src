#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <unistd.h>

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

const int IdCrdcPad[2][64] = {{30, 31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21, 18, 19, 16, 17, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 33, 32, 35, 34, 37, 36, 39, 38, 41, 40, 43, 42, 45, 44, 47, 46, 49, 48, 51, 50, 53, 52, 55, 54, 57, 56, 59, 58, 61, 60, 63, 62}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32}}; //[2]: [0]for connector=0 or 2; [1] for connector=1 or 3

void Evt2Root()
{
	string sEvt, sRoot;
	unsigned char ringHead[4];
	unsigned int evtWord;
	int run[2];
	StrtMesytec madc, mtdc, cqdcTOF, mqdcMCP;
	StrtS800 s800;
	int ringSize = 0, ringType = 0, bHSig = 0, srcID = 0;
	int nWords = 0, iCh = 0, iConnect = 0, iPad = 0, iCRDC = 0, iHodo = 0;
	int i = 0, j = 0;
	int iCtrlCRDC = 0, iEvt = 0;
	unsigned char buf[8192];
	unsigned char *pBuf;
	int runMin = 0, runMax = 0, runNum = 0;
	cout << "Input minimum and maximum numbers of run: ";
	cin >> runMin >> runMax;
	ostringstream ssRun1, ssRun2, ssIEvt;
	for (runNum = runMin; runNum <= runMax; runNum++)
	{
		ssRun1.str("");
		ssRun1 << runNum;
		ssRun2.str("");
		ssRun2 << setw(4) << setfill('0') << runNum;
		for (iEvt = 0; iEvt <= 10; iEvt++)
		{
			ssIEvt.str("");
			ssIEvt << setw(2) << setfill('0') << iEvt;

			sEvt = "../EvtData/run" + ssRun1.str() + "/run-" + ssRun2.str() + "-" + ssIEvt.str() + ".evt";
			if (access("../RootData", F_OK) != 0)
				system("mkdir ../RootData");
			sRoot = "../RootData/run-" + ssRun2.str() + "-" + ssIEvt.str() + ".root";
			printf("\n**********Now converting %s to %s!**********\n\n", sEvt.c_str(), sRoot.c_str());

			ifstream fEvt(sEvt.c_str(), ios::binary);
			if (!fEvt.is_open())
			{
				cout << "Error in Opening " << sEvt << endl;
				continue;
			}
			TFile *fRoot = new TFile(sRoot.c_str(), "RECREATE");
			TTree *tData = new TTree("tData", "tree of data");
			tData->Branch("run", run, "run[2]/I"); //[0]: run number;  [1]: sub-number of evt
			tData->Branch("madc", &madc, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("mtdc", &mtdc, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("cqdcTOF", &cqdcTOF, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("mqdcMCP", &mqdcMCP, "modID/I:data[32]/I:modRes/I:modEC_TS/I");
			tData->Branch("s800", &s800, "tS/I:eC/I:trig/I:tof[16]/I:scin[2][2]/I:ionCham[16]/I:crdcCath[2][5][64]/I:crdcAnode[2][2]/I:mesyTDC[16]/I");
			j = 0;
			run[0] = runNum;
			run[1] = iEvt;
			while (!fEvt.eof())
			{
				j++;
				//cout<<"j= "<<j<<", ";
				memset(&madc, 0, sizeof(madc));
				memset(&mtdc, 0, sizeof(mtdc));
				memset(&cqdcTOF, 0, sizeof(cqdcTOF));
				memset(&mqdcMCP, 0, sizeof(mqdcMCP));
				memset(&s800, 0, sizeof(s800));

				fEvt.read((char *)ringHead, sizeof(ringHead));
				ringSize = ringHead[0] | ringHead[1] << 8 | ringHead[2] << 16 | ringHead[3] << 24;
				if (ringSize > 8192)
					continue;

				fEvt.read((char *)buf, ringSize - 4);
				pBuf = buf;

				ringType = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
				pBuf += 4;

				if (ringType == 30)
				{
					// cout<<"ringSize= "<<ringSize<<endl;
					// cout<<"ringType= "<<ringType<<endl;
					bHSig = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					pBuf += 4;
					if (bHSig == 20)
						pBuf += 16;
					pBuf += 12;
					srcID = *pBuf | *(pBuf + 1) << 8;
					if (srcID == 2)
						continue;

					pBuf += 42;

					//Decode ADC data
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
					{
						pBuf += 4;
						madc.modID = (evtWord >> 16) & 0xFF;

						madc.modRes = (evtWord >> 12) & 0x7;
						nWords = evtWord & 0xFFF;

						//printf("ADC Header: %x\n", evtWord);

						for (i = 1; i <= nWords - 1; i++)
						{
							evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
							pBuf += 4;
							iCh = (evtWord >> 16) & 0x1F;
							madc.data[iCh] = evtWord & 0x1FFF;
							// printf("adc[%d]=%d\n", iCh, madc.data[iCh]);
						}
						evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
						pBuf += 4;
						madc.modEC_TS = evtWord & 0x3FFFFFFF;
					}
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
						continue;
					pBuf += 4;

					//Decode TDC data
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
					{
						pBuf += 4;

						mtdc.modID = (evtWord >> 16) & 0xFF;

						mtdc.modRes = (evtWord >> 12) & 0xF;
						nWords = evtWord & 0xFFF;

						//printf("TDC Header: %x\n", evtWord);

						for (i = 1; i <= nWords - 1; i++)
						{
							evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
							pBuf += 4;
							iCh = (evtWord >> 16) & 0x3F;
							if (iCh >= 0 && iCh < 32)
							{
								mtdc.data[iCh] = evtWord & 0xFFFF;
								// printf("tdc[%d]=%d\n", iCh, mtdc.data[iCh]);
							}
						}
						evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
						pBuf += 4;
						mtdc.modEC_TS = evtWord & 0x3FFFFFFF;
						//printf("TDC End: %x\n", evtWord);
					}
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
						continue;
					pBuf += 4;

					//Decode Mesytec QDC data of MCP
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
					{
						pBuf += 4;

						mqdcMCP.modID = (evtWord >> 16) & 0xFF;

						nWords = evtWord & 0xFFF;

						//printf("QDC Header: %x\n", evtWord);

						for (i = 1; i <= nWords - 1; i++)
						{
							evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
							pBuf += 4;
							iCh = (evtWord >> 16) & 0x1F;
							mqdcMCP.data[iCh] = evtWord & 0xFFF;
							//printf("qdcMCP[%d]=%d\n", iCh, mqdcMCP.data[iCh]);
						}
						evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
						pBuf += 4;
						mqdcMCP.modEC_TS = evtWord & 0x3FFFFFFF;
						//printf("QDC End: %x\n", evtWord);
					}
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
						continue;
					pBuf += 4;

					//Decode CAEN V792N QDC data of ToF
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
					{
						pBuf += 4;

						cqdcTOF.modID = (evtWord >> 24) & 0x7;

						nWords = (evtWord >> 8) & 0x3F;

						//printf("QDC Header: %x\n", evtWord);

						for (i = 1; i <= nWords; i++)
						{
							evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
							pBuf += 4;
							iCh = (evtWord >> 17) & 0xF;
							cqdcTOF.data[iCh] = evtWord & 0xFFF;
							//printf("qdcTOF[%d]=%d\n", iCh, cqdcTOF.data[iCh]);
						}
						evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
						pBuf += 4;
						cqdcTOF.modEC_TS = evtWord & 0xFFFFFF;
						//printf("QDC End: %x\n", evtWord);
					}
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					if (evtWord != 0xFFFFFFFF)
						continue;

					pBuf += 4;

					if (pBuf - buf >= ringSize - 4)
						continue;

					pBuf += 60;
					if (pBuf - buf >= ringSize - 4)
						continue;
					evtWord = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					s800.tS = evtWord & 0x3FFFFFFF;

					pBuf += 12;
					if (pBuf - buf >= ringSize - 4)
						continue;
					s800.eC = *pBuf | *(pBuf + 1) << 8 | *(pBuf + 2) << 16 | *(pBuf + 3) << 24;
					// cout<<"s800.eC="<<s800.eC<<endl;

					pBuf += 6;
					if (pBuf - buf >= ringSize - 4)
						continue;
					nWords = *pBuf | *(pBuf + 1) << 8;
					// cout<<"triggernWords="<<nWords<<endl;
					pBuf += 4;
					if (pBuf - buf >= ringSize - 4)
						continue;
					s800.trig = *pBuf | *(pBuf + 1) << 8;
					pBuf += (2 * nWords - 4); //skip trigger packet 0X5801
					if (pBuf - buf >= ringSize - 4)
						continue;

					//S800 ToF packet
					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += 4; //skip 0X000* and 0X5802
					// cout<<"nWords="<<nWords<<endl;
					if (pBuf - buf >= ringSize - 4)
						continue;
					for (i = 0; i < nWords - 2; i++)
					{
						evtWord = *pBuf | *(pBuf + 1) << 8;
						pBuf += 2; //go to next word
						if (pBuf - buf >= ringSize - 4)
							break;
						iCh = (evtWord >> 12) & 0xF;
						s800.tof[iCh] = evtWord & 0xFFF;
						// cout<<"s800.tof[iCh]="<<s800.tof[iCh]<<endl;
					}

					//S800 scintillator packet
					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += 4; //skip 0X000* and 0X5810
					if (pBuf - buf >= ringSize - 4)
						continue;
					for (i = 0; i < nWords - 2; i++)
					{
						evtWord = *pBuf | *(pBuf + 1) << 8;
						pBuf += 2; //go to next word
						if (pBuf - buf >= ringSize - 4)
							break;
						iCh = (evtWord >> 12) & 0xF;
						s800.scin[iCh][i % 2] = evtWord & 0xFFF;
						// cout<<"s800.scin[iCh][i%2]="<<s800.scin[iCh][i%2]<<endl;
					}

					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += 8; //skip 0X00**, 0X5820, 0X00**, 0X5821
					// cout<<"ionchambernWords="<<nWords<<endl;
					if (pBuf - buf >= ringSize - 4)
						continue;
					for (i = 0; i < nWords - 4; i++)
					{
						evtWord = *pBuf | *(pBuf + 1) << 8;
						pBuf += 2; //go to next word
						if (pBuf - buf >= ringSize - 4)
							break;
						iCh = (evtWord >> 12) & 0xF;
						s800.ionCham[iCh] = evtWord & 0xFFF;
					}

					for (iCRDC = 0; iCRDC < 2; iCRDC++)
					{
						pBuf += 6; //skip head of CRDC packet
						if (pBuf - buf >= ringSize - 4)
							break;
						nWords = *pBuf | *(pBuf + 1) << 8;
						pBuf += 6;
						if (pBuf - buf >= ringSize - 4)
							break;
						// cout<<"CRDCnWords="<<nWords<<endl;
						iCh = -1;
						for (i = 0; i < nWords - 3; i++)
						{
							evtWord = *pBuf | *(pBuf + 1) << 8;
							pBuf += 2;
							if (pBuf - buf >= ringSize - 4)
								break;
							iCtrlCRDC = (evtWord >> 15) & 0x1;
							if (iCtrlCRDC == 1)
							{
								iCh = evtWord & 0x3F;
								s800.crdcCath[iCRDC][0][iCh] = (evtWord >> 6) & 0x1FF;
							}
							if (iCtrlCRDC == 0 && iCh != -1)
							{
								iConnect = 1 + ((evtWord >> 10) & 0x3);
								iPad = IdCrdcPad[(iConnect - 1) % 2][iCh];
								s800.crdcCath[iCRDC][iConnect][iPad] = evtWord & 0x3FF;
							}
						}

						pBuf += 4; //skip the head of CRDC Anode sub-packet
						if (pBuf - buf >= ringSize - 4)
							break;
						s800.crdcAnode[iCRDC][0] = *pBuf | *(pBuf + 1) << 8;
						pBuf += 2;
						if (pBuf - buf >= ringSize - 4)
							break;
						s800.crdcAnode[iCRDC][1] = *pBuf | *(pBuf + 1) << 8;
						pBuf += 2;
						if (pBuf - buf >= ringSize - 4)
							break;
					}
					// S800 hodo packet 0X58b0
					for (iHodo = 0; iHodo < 2; iHodo++)
					{
						nWords = *pBuf | *(pBuf + 1) << 8;
						pBuf += (2 * nWords); //skip hodo packet 0X5870
						if (pBuf - buf >= ringSize - 4)
							break;
					}
					pBuf += 12; //skip hodo packet coincidence register
					if (pBuf - buf >= ringSize - 4)
						continue;

					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += (2 * nWords); //skip TPPACs packet 0X5870
					if (pBuf - buf >= ringSize - 4)
						continue;

					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += (2 * nWords); //skip OBJECT PIN packet 0X58a0
					if (pBuf - buf >= ringSize - 4)
						continue;

					//focal plane PIN packet 0X5805
					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += (2 * nWords); //skip focal plane PIN packet 0X5805
					if (pBuf - buf >= ringSize - 4)
						continue;

					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += (2 * nWords); //skip Galotte packet 0X58d0
					if (pBuf - buf >= ringSize - 4)
						continue;

					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += (2 * nWords); //skip LaBr packet 0X58E0
					if (pBuf - buf >= ringSize - 4)
						continue;

					//process MesytecTDC data: 0X58F0
					nWords = *pBuf | *(pBuf + 1) << 8;
					pBuf += 4;
					if (pBuf - buf >= ringSize - 4)
						continue;
					// cout<<"MTDCnWords="<<nWords<<endl;
					for (i = 0; i < nWords - 2; i++)
					{
						evtWord = *pBuf | *(pBuf + 1) << 8;
						pBuf += 2;
						if (pBuf - buf >= ringSize - 4)
							break;
						iCh = evtWord & 0xF;
						s800.mesyTDC[iCh] = *pBuf | *(pBuf + 1) << 8;
						pBuf += 2;
						if (pBuf - buf >= ringSize - 4)
							break;
					}

					cout << runNum << "--" << j << ":  madc.modEC_TS=" << madc.modEC_TS << ",  s800.tS=" << s800.tS << ",  DeltaTS=" << madc.modEC_TS - s800.tS << "  s800.trig=" << s800.trig << endl;
					tData->Fill();
				}
			}
			fRoot->Write();
			fRoot->Close();
			fEvt.close();
		}
	}
}

#ifndef __CINT__
/*
void StandaloneApplication(int argc, char** argv)
{
	Evt2Root();
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
	Evt2Root();
	return 0;
}
#endif