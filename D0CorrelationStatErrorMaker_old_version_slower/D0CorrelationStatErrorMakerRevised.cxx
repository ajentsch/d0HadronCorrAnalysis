/////ROOT LIBRARIES///////////////
#include "TROOT.h"
#include "TSystem.h" 
#include "TApplication.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <iomanip>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TPad.h>
#include <TMath.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TLeaf.h>
#include <TLatex.h>
#include "TError.h"
gErrorIgnoreLevel = 4000;  //default is kInfo

///////C++ Libraries/////////////////
#include <iostream>
#include <iomanip>
#include <math.h>
#include <map>
#include <vector>
#include <utility>
#include <climits>
#include <sstream>
#include <string>
#include <fstream>
#include <TAttMarker.h>
#include <algorithm>
#include <memory>
#include <iterator>
#include <ctype.h>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>
#include <istream>

/****
    
    Code to calculate statistical errors for D0-Hadron Correlations
    Author: Alex Jentsch -- alex.jentsch@utexas.edu
    
    Inputs: D0-Hadron singles distributions on eta and phi
    
    Outputs: Histograms of errors
    
    Last update: 9/6/2017
    
    
*****/


using namespace std;

 


int D0CorrelationStatErrorMakerRevised(){

    //constants///////////////////////////////////////
    int NUM_ETA_CORR_BINS = 9;
    int NUM_PHI_CORR_BINS = 12;
    int NUM_ETA_PRIM_BINS = 9;
    int NUM_PHI_PRIM_BINS = 24;
    
    double delEtaCorrBinCenters[10];    //NUM_ETA_CORR_BINS + 1 for same indexing
    double delPhiCorrBinCenters[13];    //NUM_PHI_CORR_BINS + 1 for same indexing
    
    double phiBinShift = TMath::Pi()/12.0;
    double phiBinShiftPrimary = TMath::Pi()/24.0;
    double NMix      = 5.0;
    double NEv[4][5][10][16];
    
    NEv[0][2][6][11] = 300;
    NEv[1][2][6][11] = 225;
    NEv[2][2][6][11] = 150;
    
    //TH2D* zPhiHist[4][5][10][16];    //this will be zPhi histogram   
    //TH2D* zEtaHist[4][5][10][16];    //this will be zEta histogram
    TH2D* mBarHist[4][5][10][16];  //this will be mBar histogram
    TH2D* nBarHist[4][5][10][16];    //this will be nBar histogram
    TH2D* errorHists[4][5][10][16];   //final errors stored here
    
    ////////////////////////////////////////////////
    
    TFile * inputFile = new TFile("final_stat_error_input_file.root");
    TFile * outputFile = new TFile("output_errors_Revised.root", "RECREATE");
    
    TString vzLabel          = "_VzBin_";
    TString multLabel        = "_CentBin_";
    TString ptBinLabel[5]    = {"0", "1", "2", "3", "4"};
    TString vzBinLabel[10]   = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    TString multBinLabel[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    
    TString zEtaLabel = "Sibling__US_Signal_Region_etaD0_vs_etaH_PtBin_";
    TString zPhiLabel = "Sibling__US_Signal_Region_phiD0_vs_phiH_PtBin_";
    TString mBarHistLabel = "Sibling__US_Signal_Region_phiD0_vs_etaD0_PtBin_";
    TString nBarHistLabel = "Sibling__US_Signal_Region_phiH_vs_etaH_PtBin_";
    TString eventCountsLabel = "event_counts_per_Vz_Cent_bin_US_Signal_Region__PtBin_";
    TString errorHistLabel = "ErrorHistogram";
    
    
    TH2D* zEtaHist = new TH2D("", "", NUM_ETA_PRIM_BINS, -1, 1, NUM_ETA_PRIM_BINS, -1, 1);
    TH2D* zPhiHist = new TH2D("", "", NUM_PHI_PRIM_BINS, -TMath::Pi()+phiBinShiftPrimary, TMath::Pi()+phiBinShiftPrimary, NUM_PHI_PRIM_BINS, -TMath::Pi()+phiBinShiftPrimary, TMath::Pi()+phiBinShiftPrimary);  
    TH2D* eventCounts = new TH2D("","",  16, 0, 16, 10, 0, 10);   
    
    TString str1;
    
    for(int band = 1; band < 2; band++){
        for(int ptBin = 0; ptBin < 5; ptBin++){
            for(int vzBin = 0; vzBin < 10; vzBin++){
                for(int multBin = 0; multBin < 16; multBin++){
                
                    
                    str1 = mBarHistLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    mBarHist[band][ptBin][vzBin][multBin] = (TH2D*) inputFile->Get(str1);
                    mBarHist[band][ptBin][vzBin][multBin]->SetDirectory(0);
                    
                    str1 = nBarHistLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    nBarHist[band][ptBin][vzBin][multBin] = (TH2D*) inputFile->Get(str1);
                    nBarHist[band][ptBin][vzBin][multBin]->SetDirectory(0);
                    
                    str1 = errorHistLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    errorHists[band][ptBin][vzBin][multBin] = new TH2D(str1, str1, NUM_ETA_CORR_BINS, -2, 2, NUM_PHI_CORR_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                
                }
            }
        }
    }
    
    for(int band = 1; band < 2; band++){
        for(int ptBin = 0; ptBin < 5; ptBin++){
        
            str1 = eventCountsLabel + ptBinLabel[ptBin];
            eventCounts = (TH2D*) inputFile->Get(str1);
            eventCounts->SetDirectory(0);
        
            for(int vzBin = 0; vzBin < 10; vzBin++){
                for(int multBin = 0; multBin < 16; multBin++){
                
                    
                    NEv[band][ptBin][vzBin][multBin] = eventCounts->GetBinContent(multBin+1,vzBin+1);
                    
                }
            }
        }
    }
    //calculated values/////////////////////////////
    
  
    
        //*************BEGIN CALCULATION OF Z VALUES FOR EACH BIN**********************//
        //calculate z-value matrix - corresponds to eta-eta and phi-phi histograms
        
        double etaZValues[10][10];
        double phiZValues[25][25];  //1 more than needed to keep indexing the same as histo-bins
        
        bool oneOrOneHalf = false; // false -- 1, true -- 1/2
        
        for(int y = 1; y < 10; y++){
            for(int x = 1; x < 10; x++){
            
                if(!oneOrOneHalf) { etaZValues[x][y] = 1.0; oneOrOneHalf = true; continue; }
                if(oneOrOneHalf)  { etaZValues[x][y] =  .5; oneOrOneHalf = false; continue; }
            }
        }
      
        oneOrOneHalf = false;
        
        for(int y = 1; y < 25; y++){         //z-values for phi
            for(int x = 1; x < 25; x++){
            
                if(!oneOrOneHalf) { phiZValues[x][y] = 1.0; oneOrOneHalf = true; continue; }
                if(oneOrOneHalf)  { phiZValues[x][y] =  .5; oneOrOneHalf = false; continue; }
            }
            if(oneOrOneHalf) { oneOrOneHalf = false; }
            else oneOrOneHalf = true;
        }
        
       
        //*************END CALCULATION OF Z VALUES FOR EACH BIN**********************//
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    
    double etaZValues[10][10];
    double phiZValues[25][25];  
    
    //Here we can fill the mBar and nBar arrays to store the occupancies. OR WE CAN JUST GRAB THEM DIRECTLY FROM THE HISTO    
    
    double zEtaDBinCenter, zEtaHBinCenter, delEtaStripe, delEtaStripeMirror;
    double zPhiDBinCenter, zPhiHBinCenter, delPhiStripe, delPhiStripeMirror;
   
    double term1 = 0;
    double term2 = 0;
    double denominator = 0;
    
    double delEtaValue;
    double delPhiValue;
    
    double NMix      = 5.0;
    double NEv[4][5][10][16];
    double mBar[10][25];
    double nBar[10][25];
    double mBarPrime;
    double nBarPrime;
    double finalErrorValueOnDelEtaDelPhi;
    bool nineOne = false;
    //double NEv
    
	int etaDBinMin, etaHBinMax, phiDBinMin, phiHBinMax, delPhiBinSymm;
	
    //TH1D * binCountEtaZ = new TH1D("check", "check", 9, -2, 2);
    
	cout << endl;
	
    ofstream testing;
    
    testing.open("testOutput.txt");
    
    for(int band = 1; band < 2; band++){
	for(int ptBin = 1; ptBin < 2; ptBin++){
	for(int vzBin = 0; vzBin < 1; vzBin++){
	for(int multBin = 5; multBin < 6; multBin++){
    
        //fill mBar and nBar arrays
        for(int i = 1; i < 10; i++){
            for(int j = 1; j < 25; j++){
            
                mBar[i][j] = mBarHist[band][ptBin][vzBin][multBin]->GetBinContent(i, j)/NEv[band][ptBin][vzBin][multBin];
                nBar[i][j] = nBarHist[band][ptBin][vzBin][multBin]->GetBinContent(i, j)/NEv[band][ptBin][vzBin][multBin];
            }
        }
		
        
        for(int delEtaBin = 2; delEtaBin < 3; delEtaBin++){
		
			if(delEtaBin == 5) {continue;}
		
                for(int delPhiBin = 1; delPhiBin < 2; delPhiBin++){
                    
                    if(delEtaBin == 1 || delEtaBin == 9){
                    
                        if(delPhiBin == 3 || delPhiBin == 9) {continue;} //skips delEta = 0, delPhi = 0, delPhi = Pi
                                
                        testing << "delEtabin: " << delEtaBin << "   delPhiBin: " << delPhiBin << endl;		
                                
                        etaDBinMin = getStartingEtaDBin(delEtaBin);
                        phiDBinMin = getStartingPhiDBin(delPhiBin);
                                    
                        etaHBinMax = getStartingEtaHBin(delEtaBin);
                        phiHBinMax = getStartingPhiHBin(delPhiBin);
                                    
                     //Special case for outer corners of eta-eta space           
                        for(int i = 0; i < 3; i++){
                            for(int etaDBin = etaDBinMin+i; etaDBin < 10; etaDBin++){
                                for(int etaHBin = 1; etaHBin < (etaHBinMax+1)-i; etaHBin++){
                                
                                    if(etaDBin == 9 && etaHBin == 1 && nineOne) {continue;}
                                    if(etaDBin == 9 && etaHBin == 1) { nineOne = true; }
                                    if(delEtaBin == 1 && delEtaBin == 9 && etaDBin == 8 && etaHBin == 2) { continue; }
                                    for(int j = 0; j < 3; j++){
                                        
                                        for(int phiDBin = phiDBinMin+j; phiDBin < 25; phiDBin++){
                                            for(int phiHBin = 1; phiHBin < (phiHBinMax+1)-j; phiHBin++){
                                
                                                testing << "(etaD, etaH): " << etaDBin << " , " << etaHBin << "   (phiD, phiH): " << phiDBin << " , " << phiHBin << endl;
                                
                                                if(mBar[etaDBin][phiDBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaDBin][etaHBin]*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*phiZValues[phiDBin][phiHBin]*nBar[etaHBin][phiHBin]*mBar[etaDBin][phiDBin]*(1-mBar[etaDBin][phiDBin]);
                                
                                                //symmeterize -- swap x and y
                                               
                                                if(mBar[etaHBin][phiHBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaHBin][etaDBin]*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*phiZValues[phiHBin][phiDBin]*nBar[etaDBin][phiDBin]*mBar[etaHBin][phiHBin]*(1-mBar[etaHBin][phiHBin]);
                                
                                                for(int g = 0; g < 3; g++){
                                
                                                    for(int etaDBinPrime = etaDBinMin+g; etaDBinPrime < 10; etaDBinPrime++){
                                                        for(int etaHBinPrime = 1; etaHBinPrime < (etaHBinMax+1)-g; etaHBinPrime++){
                                                        
                                                            for(int h = 0; h < 3; h++){
                                                                if(etaDBinPrime == 8 && etaHBinPrime == 2) { continue; }
                                                            
                                                                for(int phiDBinPrime = phiDBinMin+h; phiDBinPrime < 25; phiDBinPrime++){
                                                                    for(int phiHBinPrime = 1; phiHBinPrime < (phiHBinMax+1)-h; phiHBinPrime++){
                                        
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*etaZValues[etaDBinPrime][etaHBin]*phiZValues[phiDBinPrime][phiHBin]*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin]*nBar[etaHBinPrime][phiHBinPrime];
                                                                         
                                                                        //symmeterize -- swap x and y 
                                                                         
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*etaZValues[etaHBinPrime][etaDBin]*phiZValues[phiHBinPrime][phiDBin]*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin]*nBar[etaDBinPrime][phiDBinPrime];
                                                                         
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    
                        delPhiBinSymm = getSymmetricPhiBin(delPhiBin);
                        phiDBinMin = getStartingPhiDBin(delPhiBinSymm);
                        phiHBinMax = getStartingPhiHBin(delPhiBinSymm);
                        etaDBinMin = getStartingEtaDBin(delEtaBin);
                        etaHBinMax = getStartingEtaHBin(delEtaBin);
                        
                        nineOne = false;
                                
                        for(int a = 0; a < 3; a++){
                            for(int etaDBin = etaDBinMin+a; etaDBin < 10; etaDBin++){
                                for(int etaHBin = 1; etaHBin < (etaHBinMax+1)-a; etaHBin++){
                                
                                    if(etaDBin == 9 && etaHBin == 1 && nineOne) {continue;}
                                    if(etaDBin == 9 && etaHBin == 1) { nineOne = true; }
                                    if(delEtaBin == 1 && delEtaBin == 9 && etaDBin == 8 && etaHBin == 2) {continue;}
                                
                                    for(int b = 0; b < 3; b++){
                                        for(int phiDBin = phiDBinMin+b; phiDBin < 25; phiDBin++){
                                            for(int phiHBin = 1; phiHBin < (phiHBinMax+1)-b; phiHBin++){
                                
                                                testing << "(etaD, etaH): " << etaDBin << " , " << etaHBin << "   (phiD, phiH): " << phiDBin << " , " << phiHBin << endl;
                                
                                                if(mBar[etaDBin][phiDBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaDBin][etaHBin]*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*phiZValues[phiDBin][phiHBin]*nBar[etaHBin][phiHBin]*mBar[etaDBin][phiDBin]*(1-mBar[etaDBin][phiDBin]);
                                
                                                //symmeterize -- swap x and y
                                               
                                                if(mBar[etaHBin][phiHBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaHBin][etaDBin]*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*phiZValues[phiHBin][phiDBin]*nBar[etaDBin][phiDBin]*mBar[etaHBin][phiHBin]*(1-mBar[etaHBin][phiHBin]);
                                                        
                                                for(int c = 0; c < 3; c++){
                                
                                                    for(int etaDBinPrime = etaDBinMin+c; etaDBinPrime < 10; etaDBinPrime++){
                                                        for(int etaHBinPrime = 1; etaHBinPrime < (etaHBinMax+1)-c; etaHBinPrime++){
                                                        
                                                            for(int d = 0; d < 3; d++){
                                                                if(etaDBinPrime == 8 && etaHBinPrime == 2) { continue; }
                                                            
                                                                for(int phiDBinPrime = phiDBinMin+d; phiDBinPrime < 25; phiDBinPrime++){
                                                                    for(int phiHBinPrime = 1; phiHBinPrime < (phiHBinMax+1)-d; phiHBinPrime++){
                                        
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*etaZValues[etaDBinPrime][etaHBin]*phiZValues[phiDBinPrime][phiHBin]*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin]*nBar[etaHBinPrime][phiHBinPrime];
                                                                         
                                                                        //symmeterize -- swap x and y 
                                                                         
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*etaZValues[etaHBinPrime][etaDBin]*phiZValues[phiHBinPrime][phiDBin]*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin]*nBar[etaDBinPrime][phiDBinPrime];
                                                                         
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                            
                        cout << "term1: " << term1 << ".   term2:  " << term2 << ".    denom: " << denominator << endl;
                        
                        finalErrorValueOnDelEtaDelPhi = TMath::Sqrt((term1+term2)/denominator);
                        cout << "Final Error: "  << finalErrorValueOnDelEtaDelPhi << endl << endl;
                        errorHists[band][ptBin][vzBin][multBin]->Fill(delEtaValue, delPhiValue, finalErrorValueOnDelEtaDelPhi);
                        testing << endl << endl;
                        term1 = 0;
                        term2 = 0;
                        denominator = 0;

                        //This concludes the loop over a delEta,delPhi Pair                            
                 
                }
                
                if(delEtaBin == 2 || delEtaBin == 3 || delEtaBin == 4 || delEtaBin == 6 || delEtaBin == 7 || delEtaBin == 8){
                    
                        if(delPhiBin == 3 || delPhiBin == 9) {continue;} //skips delEta = 0, delPhi = 0, delPhi = Pi
                                
                        testing << "delEtabin: " << delEtaBin << "   delPhiBin: " << delPhiBin << endl;		
                                
                        etaDBinMin = getStartingEtaDBin(delEtaBin);
                        phiDBinMin = getStartingPhiDBin(delPhiBin);
                                    
                        etaHBinMax = getStartingEtaHBin(delEtaBin);
                        phiHBinMax = getStartingPhiHBin(delPhiBin);
                                    
                              
                        for(int i = 0; i < 3; i++){
                            for(int etaHBin = 1; etaHBin < (etaHBinMax+1)-i; etaHBin++){
                                for(int etaDBin = etaDBinMin+i; etaDBin < etaDBin+3 && etaDBin < 10; etaDBin++){
                                
                                
                                   
                                    for(int j = 0; j < 3; j++){
                                        for(int phiHBin = 1; phiHBin < (phiHBinMax+1)-j; phiHBin++){
                                            for(int phiDBin = phiDBinMin+j; phiDBin < 25; phiDBin++){
                                            
                                
                                                testing << "(etaD, etaH): " << etaDBin << " , " << etaHBin << "   (phiD, phiH): " << phiDBin << " , " << phiHBin << endl;
                                
                                                if(mBar[etaDBin][phiDBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaDBin][etaHBin]*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*phiZValues[phiDBin][phiHBin]*nBar[etaHBin][phiHBin]*mBar[etaDBin][phiDBin]*(1-mBar[etaDBin][phiDBin]);
                                
                                                //symmeterize -- swap x and y
                                               
                                                if(mBar[etaHBin][phiHBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaHBin][etaDBin]*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*phiZValues[phiHBin][phiDBin]*nBar[etaDBin][phiDBin]*mBar[etaHBin][phiHBin]*(1-mBar[etaHBin][phiHBin]);
                                
                                                for(int g = 0; g < 3; g++){
                                                    for(int etaHBinPrime = 1; etaHBinPrime < (etaHBinMax+1)-g; etaHBinPrime++){
                                                        for(int etaDBinPrime = etaDBinMin+g; etaDBinPrime < etaDBinPrime+3 && etaDBinPrime < 10; etaDBinPrime++){
                                                        
                                                            for(int h = 0; h < 3; h++){
                                                                
                                                                for(int phiHBinPrime = 1; phiHBinPrime < (phiHBinMax+1)-h; phiHBinPrime++){
                                                                    for(int phiDBinPrime = phiDBinMin+h; phiDBinPrime < 25; phiDBinPrime++){
                                                                    
                                        
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*etaZValues[etaDBinPrime][etaHBin]*phiZValues[phiDBinPrime][phiHBin]*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin]*nBar[etaHBinPrime][phiHBinPrime];
                                                                         
                                                                        //symmeterize -- swap x and y 
                                                                         
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*etaZValues[etaHBinPrime][etaDBin]*phiZValues[phiHBinPrime][phiDBin]*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin]*nBar[etaDBinPrime][phiDBinPrime];
                                                                        
                                                                         
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    
                        delPhiBinSymm = getSymmetricPhiBin(delPhiBin);
                        phiDBinMin = getStartingPhiDBin(delPhiBinSymm);
                        phiHBinMax = getStartingPhiHBin(delPhiBinSymm);
                        etaDBinMin = getStartingEtaDBin(delEtaBin);
                        etaHBinMax = getStartingEtaHBin(delEtaBin);
                        
                       
                        for(int a = 0; a < 3; a++){
                            for(int etaHBin = 1; etaHBin < (etaHBinMax+1)-a; etaHBin++){
                                for(int etaDBin = etaDBinMin+a; etaDBin < etaDBin+3 && etaDBin < 10; etaDBin++){
                                
                                    for(int b = 0; b < 3; b++){
                                        for(int phiHBin = 1; phiHBin < (phiHBinMax+1)-b; phiHBin++){
                                            for(int phiDBin = phiDBinMin+b; phiDBin < 25; phiDBin++){
                                            
                                
                                                testing << "(etaD, etaH): " << etaDBin << " , " << etaHBin << "   (phiD, phiH): " << phiDBin << " , " << phiHBin << endl;
                                
                                                if(mBar[etaDBin][phiDBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaDBin][etaHBin]*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*phiZValues[phiDBin][phiHBin]*nBar[etaHBin][phiHBin]*mBar[etaDBin][phiDBin]*(1-mBar[etaDBin][phiDBin]);
                                
                                                //symmeterize -- swap x and y
                                               
                                                if(mBar[etaHBin][phiHBin] == 0) { continue; }
                                                term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaHBin][etaDBin]*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*phiZValues[phiHBin][phiDBin]*nBar[etaDBin][phiDBin]*mBar[etaHBin][phiHBin]*(1-mBar[etaHBin][phiHBin]);
                                
                                                for(int c = 0; c < 3; c++){
                                                    for(int etaHBinPrime = 1; etaHBinPrime < (etaHBinMax+1)-c; etaHBinPrime++){
                                                        for(int etaDBinPrime = etaDBinMin+c; etaDBinPrime < etaDBinPrime+3 && etaDBinPrime < 10; etaDBinPrime++){
                                                        
                                                            for(int d = 0; d < 3; d++){
                                                                
                                                                for(int phiHBinPrime = 1; phiHBinPrime < (phiHBinMax+1)-d; phiHBinPrime++){
                                                                    for(int phiDBinPrime = phiDBinMin+d; phiDBinPrime < 25; phiDBinPrime++){
                                                                    
                                        
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*etaZValues[etaDBinPrime][etaHBin]*phiZValues[phiDBinPrime][phiHBin]*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaDBin][phiDBin]*mBar[etaDBinPrime][phiDBinPrime]*nBar[etaHBin][phiHBin]*nBar[etaHBinPrime][phiHBinPrime];
                                                                         
                                                                        //symmeterize -- swap x and y 
                                                                         
                                                                        term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaHBin][etaDBin]*phiZValues[phiHBin][phiDBin]*etaZValues[etaHBinPrime][etaDBin]*phiZValues[phiHBinPrime][phiDBin]*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin];
                                                                        denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar[etaHBin][phiHBin]*mBar[etaHBinPrime][phiHBinPrime]*nBar[etaDBin][phiDBin]*nBar[etaDBinPrime][phiDBinPrime];
                                                                         
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                            
                        cout << "term1: " << term1 << ".   term2:  " << term2 << ".    denom: " << denominator << endl;
                        
                        finalErrorValueOnDelEtaDelPhi = TMath::Sqrt((term1+term2)/denominator);
                        cout << "Final Error: "  << finalErrorValueOnDelEtaDelPhi << endl << endl;
                        errorHists[band][ptBin][vzBin][multBin]->Fill(delEtaValue, delPhiValue, finalErrorValueOnDelEtaDelPhi);
                        testing << endl << endl;
                        term1 = 0;
                        term2 = 0;
                        denominator = 0;

                        //This concludes the loop over a delEta,delPhi Pair                            
                 
                }
                 
            }
        }                        
    
                    errorHists[band][ptBin][vzBin][multBin]->Write();
    
	}
	}
	}  
    }        
    //algotrithm setup up so that each delEta, delPhi bin in correlation space gives the band(s) in primary space.
    
    
    
    //textfile.close();
    testing.close();
    outputFile->Close();
    inputFile->Close();

}




	
int getStartingPhiDBin(int delPhiBin){ //PhiD bin is "x" in my notebook notation

	if(delPhiBin == 2)  { return 2;  }  //if we get bin 2, call function to 
    if(delPhiBin == 4)  { return 22; }
    
	if(delPhiBin == 1)  { return 4;  }
	if(delPhiBin == 5)  { return 20; }
	
	if(delPhiBin == 6)  { return 6;  }
	if(delPhiBin == 12) { return 18; }
	
	if(delPhiBin == 7)  { return 8;  }
	if(delPhiBin == 11) { return 16; }
	
	if(delPhiBin == 8)  { return 10; }
	if(delPhiBin == 10) { return 14; }
	
	if(delPhiBin == 3)  { return 2;  } //special case because of end corner bins -- be careful

	if(delPhiBin == 9)  { return 12; }
	
	else return -1;
	
}

int getStartingPhiHBin(int delPhiBin){ //PhiH bin is "y" in my notebook notation

	if(delPhiBin == 2)  { return 23;  }  //if we get bin 2, call function to 
    if(delPhiBin == 4)  { return 3; }
    
	if(delPhiBin == 1)  { return 21;  }
	if(delPhiBin == 5)  { return 5; }
	
	if(delPhiBin == 6)  { return 19;  }
	if(delPhiBin == 12) { return 7; }
	
	if(delPhiBin == 7)  { return 17;  }
	if(delPhiBin == 11) { return 9; }
	
	if(delPhiBin == 8)  { return 15; }
	if(delPhiBin == 10) { return 11; }
	
	if(delPhiBin == 3)  { return 2;  } //special case because of end corner bins -- be careful

	if(delPhiBin == 9)  { return 13; }
	
	else return -1;
	
}

int getStartingEtaDBin(int delEtaBin){

	if(delEtaBin == 4 || delEtaBin == 6) { return 2; }
	if(delEtaBin == 3 || delEtaBin == 7) { return 4; }
	if(delEtaBin == 2 || delEtaBin == 8) { return 6; }
	if(delEtaBin == 1 || delEtaBin == 9) { return 8; }
	
	else return -1;
	
}

int getStartingEtaHBin(int delEtaBin){

	if(delEtaBin == 4 || delEtaBin == 6) { return 8; }
	if(delEtaBin == 3 || delEtaBin == 7) { return 6; }
	if(delEtaBin == 2 || delEtaBin == 8) { return 4; }
	if(delEtaBin == 1 || delEtaBin == 9) { return 2; }
	
	else return -1;
	
}

int getSymmetricDelEtaBin(int delEtaBin){

	if( delEtaBin == 1) { return 9; }
	if( delEtaBin == 2) { return 8; }
	if( delEtaBin == 3) { return 7; }
	if( delEtaBin == 4) { return 6; }
	//if( delEtaBin == 5) { return 5; }
	if( delEtaBin == 6) { return 4; }
	if( delEtaBin == 7) { return 3; }
	if( delEtaBin == 8) { return 2; }
	if( delEtaBin == 9) { return 1; }
	
	else return -1;
	
}

int getSymmetricPhiBin(int delPhiBin){

	if(delPhiBin == 2)  { return 4;  }  //if we get bin 2, call function to 
    if(delPhiBin == 4)  { return 2; }
    
	if(delPhiBin == 1)  { return 5;  }
	if(delPhiBin == 5)  { return 1; }
	
	if(delPhiBin == 6)  { return 12;  }
	if(delPhiBin == 12) { return 6; }
	
	if(delPhiBin == 7)  { return 11;  }
	if(delPhiBin == 11) { return 7; }
	
	if(delPhiBin == 8)  { return 10; }
	if(delPhiBin == 10) { return 8; }
	
	//if(delPhiBin == 3)  { return 3;  } //special case because of end corner bins -- be careful

	//if(delPhiBin == 9)  { return 9; }
	
	else return -1;
}