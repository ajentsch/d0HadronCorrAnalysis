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

//threshold = .00001

bool isZero(double value, double threshold){
    return value >= -threshold && value <= threshold;
}

bool isDiffZero(double a, double b, double threshold){
    
    double value = a - b;
    
    return value >= -threshold && value <= threshold;
}

bool isSumLessThanZero(double a, double b, double threshold){

    double value = a + b;
    
    if(value < -threshold) { return true; }
    
    else return false;
}    

bool isDiffGreaterThanZero(double a, double b, double threshold){

    double value = a - b;
    
    if(value > threshold) { return true; }
    
    else return false;
}    

bool isDiffGreaterThanOrEqualToZero(double a, double b, double threshold){

    double value = a - b;
    
    if(value >= threshold) { return true; }
    
    else return false;
}    


int D0CorrelationStatErrorMaker(){

    //constants///////////////////////////////////////
    int NUM_ETA_CORR_BINS = 9;
    int NUM_PHI_CORR_BINS = 12;
    int NUM_ETA_PRIM_BINS = 9;
    int NUM_PHI_PRIM_BINS = 24;
    
    double delEtaCorrBinCenters[10];    //NUM_ETA_CORR_BINS + 1 for same indexing
    double delPhiCorrBinCenters[13];    //NUM_PHI_CORR_BINS + 1 for same indexing
    double delPhiZBinCenters[25];    //NUM_PHI_CORR_BINS + 1 for same indexing
    
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
    
    TFile * inputFile = new TFile("best_output_for_stat_errors.root");
    TFile * outputFile = new TFile("output_errors_correct_6.root", "RECREATE");
    
    TString vzLabel          = "_VzBin_";
    TString multLabel        = "_CentBin_";
    TString siblingLabel     = "Sibling__";
    TString bandLabel[3]     = {"SBR_", "US_Signal_Region_", "SBL_"};
    TString ptBinLabel[5]    = {"0", "1", "2", "3", "4"};
    TString vzBinLabel[10]   = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    TString multBinLabel[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    
    //TString zEtaLabel = "Signal_Region_etaD0_vs_etaH_PtBin_";
    //TString zPhiLabel = "Signal_Region_phiD0_vs_phiH_PtBin_";
    TString mBarHistLabel = "phiD0_vs_etaD0_PtBin_";
    TString nBarHistLabel = "phiH_vs_etaH_PtBin_";
    TString eventCountsLabel = "event_counts_per_Vz_Cent_bin_US_Signal_Region__PtBin_";
    TString errorHistLabel = "ErrorHistogram";
    
    
    TH2D* zEtaHist = new TH2D("", "", NUM_ETA_PRIM_BINS, -1, 1, NUM_ETA_PRIM_BINS, -1, 1);
    //TH2D* zPhiHist = new TH2D("", "", NUM_PHI_PRIM_BINS, -TMath::Pi()+phiBinShiftPrimary, TMath::Pi()+phiBinShiftPrimary, NUM_PHI_PRIM_BINS, -TMath::Pi()+phiBinShiftPrimary, TMath::Pi()+phiBinShiftPrimary);  
    TH2D* eventCounts = new TH2D("","",  16, 0, 16, 10, 0, 10);   
    
    TString str1;
    
    for(int band = 0; band < 3; band++){
        for(int ptBin = 0; ptBin < 5; ptBin++){
            for(int vzBin = 0; vzBin < 10; vzBin++){
                for(int multBin = 0; multBin < 16; multBin++){
                
                    
                    str1 = siblingLabel + bandLabel[band] + mBarHistLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    mBarHist[band][ptBin][vzBin][multBin] = (TH2D*) inputFile->Get(str1);
                    mBarHist[band][ptBin][vzBin][multBin]->SetDirectory(0);
                    
                    str1 = siblingLabel + bandLabel[band] + nBarHistLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    nBarHist[band][ptBin][vzBin][multBin] = (TH2D*) inputFile->Get(str1);
                    nBarHist[band][ptBin][vzBin][multBin]->SetDirectory(0);
                    
                    str1 = errorHistLabel + ptBinLabel[ptBin] + vzLabel + vzBinLabel[vzBin] + multLabel + multBinLabel[multBin];
                    errorHists[band][ptBin][vzBin][multBin] = new TH2D(str1, str1, NUM_ETA_CORR_BINS, -2, 2, NUM_PHI_CORR_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                
                }
            }
        }
    }
    
    for(int band = 0; band < 3; band++){
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
    
    
    double delEtaBinWidthCorr = 4.0/NUM_ETA_CORR_BINS;
    double firstDelEtaBin = -2+.5*(delEtaBinWidthCorr); 
    
    for(int i = 1; i < NUM_ETA_CORR_BINS+1; i++){
    
        delEtaCorrBinCenters[i] = firstDelEtaBin + (i-1)*delEtaBinWidthCorr;
        cout << delEtaCorrBinCenters[i] << endl;
        
    }
    
    cout << endl;
    
    double delPhiBinWidthCorr = TMath::TwoPi()/NUM_PHI_CORR_BINS;
    double firstDelPhiBin = -8.0;
    double firstZPhiBin = -22.0;
    
    for(int i = 1; i < NUM_PHI_CORR_BINS+1; i++){
    
        delPhiCorrBinCenters[i] = firstDelPhiBin + (i-1.0)*(4.0);
        //if(isZero(delPhiCorrBinCenters[i], .00001)) { delPhiCorrBinCenters[i] = 0.0; } 

        cout << delPhiCorrBinCenters[i] << endl;
        
    }
    
    for(int i = 1; i < NUM_PHI_PRIM_BINS+1; i++){
    
        delPhiZBinCenters[i] = firstZPhiBin + (i-1.0)*(2.0);
        //if(isZero(delPhiCorrBinCenters[i], .00001)) { delPhiCorrBinCenters[i] = 0.0; } 

        cout << delPhiZBinCenters[i] << endl;
        
    }
    
    
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
    double mBar;
    double nBar;
    double mBarPrime;
    double nBarPrime;
    double finalDelPhiValue;
    double finalErrorValueOnDelEtaDelPhi;
    //double NEv
    
    TH1D * binCountEtaZ = new TH1D("check", "check", 9, -2, 2);
    
    int numPhiBins  = 0;
    
    for(int band = 1; band < 2; band++){
        for(int ptBin = 2; ptBin < 3; ptBin++){
            for(int vzBin = 5; vzBin < 6; vzBin++){
                for(int multBin = 10; multBin < 11; multBin++){
    
                   // for(int delEtaBin = 1; delEtaBin < NUM_ETA_CORR_BINS+1; delEtaBin++){
                       // for(int delPhiBin = 1; delPhiBin < NUM_PHI_CORR_BINS+1; delPhiBin++){
                    
                    for(int delEtaBin = 1; delEtaBin < 10; delEtaBin++){
                        for(int delPhiBin = 1; delPhiBin < 13; delPhiBin++){
                    
                            //ofstream textfile;
                            //textfile.open("output.txt");
                    
                            //for(int delEtaBin = 1; delEtaBin < 2; delEtaBin++){
                            //for(int delPhiBin = 12; delPhiBin < 13; delPhiBin++){
                    
                            delEtaValue = delEtaCorrBinCenters[delEtaBin];            //this is the actual delEta value -- gives the "stripe(s)"
                            delPhiValue = delPhiCorrBinCenters[delPhiBin];            //this is the actual delPhi value -- gives the "stripe(s)"
                    
                            cout << "delEtaBin & value: " << delEtaBin << ", " << delEtaValue << endl;
                            cout << "delPhiBin & value: " << delPhiBin << ", " << delPhiValue << endl;
                            
                            for(int etaDBin = 1; etaDBin < NUM_ETA_PRIM_BINS+1; etaDBin++){
                                for(int etaHBin = 1; etaHBin < NUM_ETA_PRIM_BINS+1; etaHBin++){
                                
                                    if(numPhiBins != 0) {cout << "numPhiBins: "<< numPhiBins << endl;}
                                
                                    numPhiBins = 0;
                                
                                    for(int phiDBin = 1; phiDBin < NUM_PHI_PRIM_BINS+1; phiDBin++){
                                        for(int phiHBin = 1; phiHBin < NUM_PHI_PRIM_BINS+1; phiHBin++){
                                        
                                        
                                            mBar = mBarHist[band][ptBin][vzBin][multBin]->GetBinContent(etaDBin, phiDBin)/NEv[band][ptBin][vzBin][multBin];
                                            nBar = nBarHist[band][ptBin][vzBin][multBin]->GetBinContent(etaHBin, phiHBin)/NEv[band][ptBin][vzBin][multBin];
                                            
                                            //cout << mBar << "     " << nBar << endl;
                                            
                                            if(mBar == 0 || nBar == 0) { continue; }
                                            
                                            zEtaDBinCenter = zEtaHist->GetXaxis()->GetBinCenter(etaDBin);  //gets eta value for D bin center
                                            zEtaHBinCenter = zEtaHist->GetYaxis()->GetBinCenter(etaHBin);  //gets eta value for H bin center
                                            delEtaStripe = TMath::Abs(zEtaDBinCenter - zEtaHBinCenter);            //gets delEta from those values
                                            delEtaStripeMirror = -delEtaStripe;             //gets the other stripe for the symmetric delEta bin
                                            
                                               
                                            if(!(isDiffZero(delEtaStripe, delEtaValue, .0001)) && !(isDiffZero((delEtaStripe - .5*delEtaBinWidthCorr), delEtaValue, .0001)) && !(isDiffZero((delEtaStripe + .5*delEtaBinWidthCorr), delEtaValue, .0001)) &&
                                               !(isDiffZero(delEtaStripeMirror, delEtaValue, .0001)) && !(isDiffZero((delEtaStripeMirror - .5*delEtaBinWidthCorr), delEtaValue, .0001)) && !(isDiffZero((delEtaStripeMirror + .5*delEtaBinWidthCorr), delEtaValue, .0001))) {continue;}   
                                               
                                            zPhiDBinCenter = delPhiZBinCenters[phiDBin];  //gets phi value for D bin center
                                            zPhiHBinCenter = delPhiZBinCenters[phiHBin];  //gets phi value for H bin center
                                            delPhiStripe = zPhiDBinCenter - zPhiHBinCenter;  //gets absolute delPhi from those values
                                            
                                            //cout << "1: " << delPhiStripe << endl;
                                            
                                            if(delPhiStripe < -22.0) { delPhiStripe = delPhiStripe + 48.0; } //shifts [-2pi, 2pi] -> [-pi, pi]
                                            else if(delPhiStripe > 26.0) { delPhiStripe = delPhiStripe - 48.0; } //shifts [-2pi, 2pi] -> [-pi, pi]
                                            
                                            delPhiStripe = TMath::Abs(delPhiStripe);
                                            
                                            //cout <<  "2: " << delPhiStripe << endl;
                                            
                                            if((delPhiValue == 36.0) && (delPhiStripe >= 38.0)) { delPhiStripe = delPhiStripe - 48.0; } 
                                            else if(delPhiValue != 36.0 && delPhiStripe > 38.0) { delPhiStripe = delPhiStripe - 48.0; }
                                            
                                            if(delPhiStripe == 0.0) { delPhiStripeMirror = 0.0; }
                                            else delPhiStripeMirror = -delPhiStripe;
                                            
                                            //cout << "3: " << delPhiStripeMirror << endl;
                                            
                                            if((delPhiValue == 36.0) && (delPhiStripeMirror <= -10.0)) { delPhiStripeMirror = delPhiStripeMirror + 48.0; } 
                                            else if(delPhiValue != 36.0 && delPhiStripeMirror < -10.0) { delPhiStripeMirror = delPhiStripeMirror + 48.0; } 
                                           
                                            if(!(delPhiStripe == delPhiValue) && !((delPhiStripe - 2.0) == delPhiValue) && !(((delPhiStripe + 2.0) == delPhiValue)) &&
                                               !(delPhiStripeMirror == delPhiValue) && !((delPhiStripeMirror - 2.0) == delPhiValue) && !(((delPhiStripeMirror + 2.0) == delPhiValue)) ) {continue;}
                                            
                                            //cout << "eta-eta Bin: " << etaDBin << "  , " << etaHBin << endl;
                                            //textfile << "phi-phi Bin: " << phiDBin << "  , " << phiHBin << endl;
                                            numPhiBins++;
                                           
                                            term1 = term1 + NEv[band][ptBin][vzBin][multBin]*NMix*(1+NMix)*etaZValues[etaDBin][etaHBin]*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*phiZValues[phiDBin][phiHBin]*nBar*mBar*(1-mBar);
                    
                                            binCountEtaZ->Fill(delEtaBin);
                    
                                            for(int etaDBinPrime = 1; etaDBinPrime < NUM_ETA_PRIM_BINS+1; etaDBinPrime++){
                                                for(int etaHBinPrime = 1; etaHBinPrime < NUM_ETA_PRIM_BINS+1; etaHBinPrime++){
                                                    for(int phiDBinPrime = 1; phiDBinPrime < NUM_PHI_PRIM_BINS+1; phiDBinPrime++){
                                                        for(int phiHBinPrime = 1; phiHBinPrime < NUM_PHI_PRIM_BINS+1; phiHBinPrime++){
                    
                                                            mBarPrime = mBarHist[band][ptBin][vzBin][multBin]->GetBinContent(etaDBinPrime, phiDBinPrime)/NEv[band][ptBin][vzBin][multBin];
                                                            nBarPrime = nBarHist[band][ptBin][vzBin][multBin]->GetBinContent(etaHBinPrime, phiHBinPrime)/NEv[band][ptBin][vzBin][multBin];
            
                                                            if(mBarPrime == 0 || nBarPrime == 0) { continue; }
                    
                                                            zEtaDBinCenter = zEtaHist->GetXaxis()->GetBinCenter(etaDBinPrime);  //gets eta value for D bin center
                                                            zEtaHBinCenter = zEtaHist->GetYaxis()->GetBinCenter(etaHBinPrime);  //gets eta value for H bin center
                                                            delEtaStripe = TMath::Abs(zEtaDBinCenter - zEtaHBinCenter);            //gets delEta from those values
                                                            delEtaStripeMirror = -delEtaStripe;             //gets the other stripe for the symmetric delEta bin
                                                            
                                                            if(!(isDiffZero(delEtaStripe, delEtaValue, .0001)) && !(isDiffZero((delEtaStripe - .5*delEtaBinWidthCorr), delEtaValue, .0001)) && !(isDiffZero((delEtaStripe + .5*delEtaBinWidthCorr), delEtaValue, .0001)) &&
                                                               !(isDiffZero(delEtaStripeMirror, delEtaValue, .0001)) && !(isDiffZero((delEtaStripeMirror - .5*delEtaBinWidthCorr), delEtaValue, .0001)) && !(isDiffZero((delEtaStripeMirror + .5*delEtaBinWidthCorr), delEtaValue, .0001))) {continue;}                                    
                                                            
                                                            zPhiDBinCenter = delPhiZBinCenters[phiDBinPrime];  //gets phi value for D bin center
                                                            zPhiHBinCenter = delPhiZBinCenters[phiHBinPrime];  //gets phi value for H bin center
                                                            delPhiStripe = zPhiDBinCenter - zPhiHBinCenter;  //gets absolute delPhi from those values
                                                            
                                                            if(delPhiStripe < -24.0) { delPhiStripe = delPhiStripe + 50.0; } //shifts [-2pi, 2pi] -> [-pi, pi]
                                                            else if(delPhiStripe > 24.0) { delPhiStripe = delPhiStripe - 46.0; } //shifts [-2pi, 2pi] -> [-pi, pi]
                                                            
                                                            delPhiStripe = TMath::Abs(delPhiStripe);
                                                            
                                                            if((delPhiValue == 36.0) && (delPhiStripe >= 38.0)) { delPhiStripe = delPhiStripe - 48.0; } 
                                                            else if(delPhiValue != 36.0 && delPhiStripe > 38.0) { delPhiStripe = delPhiStripe - 48.0; }
                                                            
                                                            if(delPhiStripe == 0.0) { delPhiStripeMirror = 0.0; }
                                                            else delPhiStripeMirror = -delPhiStripe;
                                                            
                                                            if((delPhiValue == 36.0) && (delPhiStripeMirror <= -10.0)) { delPhiStripeMirror = delPhiStripeMirror + 48.0; } 
                                                            else if(delPhiValue != 36.0 && delPhiStripeMirror < -10.0) { delPhiStripeMirror = delPhiStripeMirror + 48.0; } 
                                                           
                                                            if(!(delPhiStripe == delPhiValue) && !((delPhiStripe - 2.0) == delPhiValue) && !(((delPhiStripe + 2.0) == delPhiValue)) &&
                                                               !(delPhiStripeMirror == delPhiValue) && !((delPhiStripeMirror - 2.0) == delPhiValue) && !(((delPhiStripeMirror + 2.0) == delPhiValue)) ) {continue;}
                                                            
                                                            
                                                            term2 = term2 + 2*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*etaZValues[etaDBin][etaHBin]*phiZValues[phiDBin][phiHBin]*etaZValues[etaDBinPrime][etaHBinPrime]*phiZValues[phiDBinPrime][phiHBinPrime]*mBar*mBarPrime*nBar;
                                                            
                                                            denominator = denominator + NEv[band][ptBin][vzBin][multBin]*NEv[band][ptBin][vzBin][multBin]*NMix*NMix*mBar*mBarPrime*nBar*nBarPrime;
                                                            
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
                            finalDelPhiValue = delPhiValue*(TMath::Pi()/24.0);
                            errorHists[band][ptBin][vzBin][multBin]->Fill(delEtaValue, finalDelPhiValue, finalErrorValueOnDelEtaDelPhi);
                            
                            term1 = 0;
                            term2 = 0;
                            denominator = 0;

                            //This concludes the loop over a delEta,delPhi Pair                            
                            
                        }
                    }                        
    
                    errorHists[band][ptBin][vzBin][multBin]->Write();
    
                }
            }
        }  
    }        
    //algotrithm setup up so that each delEta, delPhi bin in correlation space gives the band(s) in primary space.
    
    binCountEtaZ->Write();
    
    //textfile.close();
    outputFile->Close();
    inputFile->Close();

}



