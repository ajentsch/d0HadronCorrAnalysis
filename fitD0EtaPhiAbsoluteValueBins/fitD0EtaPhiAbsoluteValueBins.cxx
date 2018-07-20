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
//gErrorIgnoreLevel = 4000;  //default is kInfo

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



using namespace std;


int fitD0EtaPhiAbsoluteValueBins(){

    TFile * inputFile = new TFile("d0HadronCorrMaker_OUTPUT_pt_2_10_GeV_nominal_cuts.root");
    TFile * outputFile = new TFile("ABS_FITTING_OUTPUT_d0HadronCorrMaker_OUTPUT_pt_2_10_GeV_nominal_cuts.root", "RECREATE");
    
    double phiBinShift = (TMath::Pi()/12.0);
    int NUM_PT_BINS = 3;
    double alpha = .05; //confidence interval for p-value calculation
	
	double NUM_SYMM_ETA_BINS = 13;
	double NUM_SYMM_PHI_BINS = 12;
	int NUM_ABS_VAL_ETA_BINS = 6;
	int NUM_ABS_VAL_PHI_BINS = 7;
	double LOWER_ETA_BIN = -0.153846;
	double UPPER_ETA_BIN = 1.69+0.153846;
    
    bool OFFSET_ONLY = false;
    bool OFFSET_PLUS_DIPOLE = false;
    bool OFFSET_PLUS_DIPOLE_PLUS_QUADRUPOLE = false;
    bool FULL_FIT_WITH_DIPOLE = false;                      //mid-central
    bool FULL_FIT_NO_DIPOLE = false;
    bool FULL_FIT_WITH_AWAY_SIDE_GAUSS = false;
    bool FULL_FIT_WITH_DIPOLE_AND_ETA_GAUSSIAN = false;    //Lanny suggests for central bin 5-15-2018
    bool DOUBLE_2D_GAUSS = false;
    bool OFFSET_QUAD_JETPEAK_AS_1D_GAUSS = true; //Lanny suggests for peripheral bin
    
    cout << endl;
    
    if(OFFSET_ONLY) { cout << "Using option: OFFSET ONLY (28 - 1) DF" << endl; }
    else if(OFFSET_PLUS_DIPOLE) { cout << "Using option: OFFSET PLUS DIPOLE ONLY (28 - 2) DF" << endl; }
    else if(OFFSET_PLUS_DIPOLE_PLUS_QUADRUPOLE) { cout << "Using option: OFFSET PLUS DIPOLE PLUS QUADRUPOLE ONLY (28 - 3) DF" << endl; }
    else if(FULL_FIT_NO_DIPOLE) { cout << "Using option: FULL FIT NO DIPOLE (28 - 5) DF" << endl; }
    else if(FULL_FIT_WITH_DIPOLE) { cout << "Using option: FULL FIT USING DIPOLE (28 - 6) DF" << endl; }
    else if(FULL_FIT_WITH_AWAY_SIDE_GAUSS) { cout << "Using option: FULL FIT USING AS GAUSSIAN (28 - 7) DF" << endl; }
    else if(FULL_FIT_WITH_DIPOLE_AND_ETA_GAUSSIAN) { cout << "Using option: FULL FIT USING DIPOLE & ETA_DEPENDENT GAUSSIAN (35 - 7) DF" << endl; }
    else if(DOUBLE_2D_GAUSS) { cout << "Using option: FULL FIT COMPLICATED JET SHAPE (28 - 10) DF" << endl; }
    else if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { cout << "Using option: Lanny Fit for Peripheral (35 - 7) = 21 DF " << endl; }
    
    cout << endl;
    
    
    TH1D * fitParameterOffset[5];
    TH1D * fitParameterDipole[5];
    TH1D * fitParameterQuad[5];
    TH1D * fitParameterJetAmp[5];
    TH1D * fitParameterSigEta[5];
    TH1D * fitParameterSigPhi[5];
    
    TH1D * jetVolumePlot[5];
    
    TH1D * fitParameterASGaussWidth[5];
    TH1D * fitParameterASGaussAmp[5];
    
    TString offsetLabel = "Offset_A0_";
    TString dipoleLabel = "Dipole_AD_";
    TString quadLabel   = "Quadrupole_AQ_";
    TString jetAmpLabel = "JetAmp_A1_";
    TString sigEtaLabel = "SigmaEta_";
    TString sigPhiLabel = "SigmaPhi_";
    TString ASGaussWidthLabel = "AwaySide_Gauss_Width";
    TString ASGaussAmpLabel = "AwaySide_Gauss_Amp";
    
    TH2D * inputCorrHist[5][3];
    TH1D * inputCorrHistPhiProjection[5][3];
    TH1D * inputCorrHistPhiProjectionFit[5][3];

    TH2D * absCorrHist[5][3];
    TH2D * absCorrHistFit[5][3];
    TH2D * absCorrHistRes[5][3];
    TH2D * chiSquareHist[5][3];
    TString str1;
    TString str2;
    TF2 * myfit[5][3];
    
    double nBarPrime[3] = {74.4, 351, 1027};
    
    //TString corrHistLabel   = "FullSubtractedCorr_SideBand__"; 
	TString corrHistLabel   = "Full_DStar_Correction_FullSubtractedCorr_SideBand__"; 
    TString ptBinLabel[5]   = {"PtBin_0", "PtBin_1", "PtBin_2", "PtBin_3", "PtBin_4"};
    TString centBinLabel[3] = {"Peripheral", "MidCentral", "Central"};
    TString absValueBin1     = "AbsoulteValue_";
    TString absValueBin2     = "#||{#Delta#eta}, #||{#Delta#phi}_";
    TString fitFunctionLabel = "FitFunction";
    TString fitLabel         = "_Fit";
    TString resLabel         = "_Residual";
    TString chiSquareLabel   = "_ChiSquare";
    TString pngLabel         = ".png";
    
    //read-in final, symmeterized correlations
    
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        for(int centBin = 0; centBin < 3; centBin++){
        
            str1 = corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin];
            inputCorrHist[ptBin][centBin] = (TH2D*) inputFile->Get(str1);
            inputCorrHist[ptBin][centBin]->SetDirectory(0);
        }
    }
    
    TString histogramProjectionLabel = "phiProjection_";
    
    TCanvas * projectionCan = new TCanvas("", "", 3300, 1100);
    projectionCan->Divide(3,1);


    
    
    //create histograms for absolute value data
    
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        for(int centBin = 0; centBin < 3; centBin++){
    
            str1 = absValueBin1 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin];
            str2 = absValueBin2 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin];
        
            absCorrHist[ptBin][centBin] = new TH2D(str1, str2, NUM_ABS_VAL_ETA_BINS, LOWER_ETA_BIN, UPPER_ETA_BIN, NUM_ABS_VAL_PHI_BINS, -phiBinShift, TMath::Pi()+phiBinShift);  //changed binning to direct copy
            
            str1 = absValueBin1 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin] + fitLabel;
            str2 = absValueBin2 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin] + fitLabel;
            
            absCorrHistFit[ptBin][centBin] = new TH2D(str1, str2, NUM_ABS_VAL_ETA_BINS, LOWER_ETA_BIN, UPPER_ETA_BIN, NUM_ABS_VAL_PHI_BINS, -phiBinShift, TMath::Pi()+phiBinShift);
            
            str1 = absValueBin1 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin] + resLabel;
            str2 = absValueBin2 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin] + resLabel;
            
            absCorrHistRes[ptBin][centBin] = new TH2D(str1, str2, NUM_ABS_VAL_ETA_BINS, LOWER_ETA_BIN, UPPER_ETA_BIN, NUM_ABS_VAL_PHI_BINS, -phiBinShift, TMath::Pi()+phiBinShift);
            
            str1 = absValueBin1 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin] + chiSquareLabel;
            str2 = absValueBin2 + corrHistLabel + ptBinLabel[ptBin] + centBinLabel[centBin] + chiSquareLabel;
            
            chiSquareHist[ptBin][centBin] = new TH2D(str1, str2, NUM_ABS_VAL_ETA_BINS, LOWER_ETA_BIN, UPPER_ETA_BIN, NUM_ABS_VAL_PHI_BINS, -phiBinShift, TMath::Pi()+phiBinShift);
            
            absCorrHist[ptBin][centBin]->SetStats(0);
            absCorrHistFit[ptBin][centBin]->SetStats(0);
            absCorrHistRes[ptBin][centBin]->SetStats(0);
            chiSquareHist[ptBin][centBin]->SetStats(0);
            //absCorrHist[ptBin][centBin]->Write();

        }
    }

    //parse symmeterized histograms to fill absolute value histograms
    
    double sqrt2 = TMath::Sqrt(2);
    
    
    
    int corrEtaBin = 0;
    int corrPhiBin = 0;
    
    
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        for(int centBin = 0; centBin < 3; centBin++){
            for(int etaBin = 1; etaBin < NUM_ABS_VAL_ETA_BINS+1; etaBin++){
                for(int phiBin = 1; phiBin < NUM_ABS_VAL_PHI_BINS+1; phiBin++){
                
                    corrEtaBin = getEtaBin(etaBin);
                    corrPhiBin = getPhiBin(phiBin);
               
                    if(centBin == 1 && etaBin == 1 && phiBin == 1){absCorrHist[ptBin][centBin]->SetBinError(etaBin, phiBin, 1.2*inputCorrHist[ptBin][centBin]->GetBinError(corrEtaBin,corrPhiBin));}
                    else absCorrHist[ptBin][centBin]->SetBinError(etaBin, phiBin, inputCorrHist[ptBin][centBin]->GetBinError(corrEtaBin,corrPhiBin));
                    
                    absCorrHist[ptBin][centBin]->SetBinContent(etaBin, phiBin, inputCorrHist[ptBin][centBin]->GetBinContent(corrEtaBin, corrPhiBin));
                    
                    
                }
            }
            
            absCorrHist[ptBin][centBin]->Write();

        }
    } 
    
  
    for(int k = 2; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){ 
 
            str1 = fitFunctionLabel + ptBinLabel[k] + centBinLabel[i];
            if(OFFSET_ONLY) {myfit[k][i] = new TF2(str1,"[0]", LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift); }
            else if(OFFSET_PLUS_DIPOLE) {myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y)",  LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift); }
            else if(OFFSET_PLUS_DIPOLE_PLUS_QUADRUPOLE) {myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y)", 0, UPPER_ETA_BIN, 0 , TMath::Pi()); }
            else if(FULL_FIT_NO_DIPOLE) {myfit[k][i] = new TF2(str1,"[0] + [1]*2*cos(2*y) + [2]*(exp(-0.5*(y*y)/([3]*[3]))+exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4]))", LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift); }
            else if(FULL_FIT_WITH_DIPOLE) {myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift); }
            else if(FULL_FIT_WITH_DIPOLE_AND_ETA_GAUSSIAN) {myfit[k][i] = new TF2(str1,"[0] + [1]*2*cos(2*y) + [2]*(exp(-0.5*(y*y)/([3]*[3]))+exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4])) + [5]*cos(y)*exp(-0.5*(x*x)/([6]*[6]))", LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift); }
            else if(FULL_FIT_WITH_AWAY_SIDE_GAUSS) {myfit[k][i] = new TF2(str1,"[0] + [1]*(exp(-0.5*((y-3.14)*(y-3.14))/([2]*[2]))+exp(-0.5*((y+3.14)*(y+3.14))/([2]*[2]))) + [3]*2*cos(2*y) + [4]*(exp(-0.5*(y*y)/([5]*[5]))+exp(-0.5*((y-6.28)*(y-6.28))/([5]*[5])))*exp(-0.5*(x*x)/([6]*[6]))",  LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift);}
            else if(DOUBLE_2D_GAUSS) {myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + ([3]*exp(-0.5*(x*x)/([4]*[4]))+([5]*exp(-0.5*(x*x)/([6]*[6]))))*(exp(-0.5*(y*y)/([7]*[7]))+exp(-0.5*((y-6.28)*(y-6.28))/([7]*[7]))) + [8]*exp(-0.5*(x*x)/([9]*[9]))", LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift); }
            else if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { myfit[k][i] = new TF2(str1,"[0] + 2*[1]*cos(2*y) + [2]*(exp(-0.5*(y*y)/([3]*[3])) + exp(-0.5*((y-6.283185)*(y-6.283185))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4])) + [5]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([6]*[6]))+exp(-0.5*((y+3.14159)*(y+3.14159))/([6]*[6])))", LOWER_ETA_BIN, UPPER_ETA_BIN, -phiBinShift, TMath::Pi()+phiBinShift); }
		}
	}
    
    if(OFFSET_ONLY){
    
        //myfit[2][0]->SetParLimits(0, -.04, .04);    //offset A0
        //myfit[2][1]->SetParLimits(0, -.04, .04);    //offset A0
        //myfit[2][2]->SetParLimits(0, -.04, .04);    //offset A0
        
    }
    
    else if(OFFSET_PLUS_DIPOLE){
    
        myfit[2][0]->SetParLimits(0, -.04, .04);    //offset A0
        myfit[2][0]->SetParLimits(1, -.04, .04);    //Dipole AD
        
        myfit[2][1]->SetParLimits(0, -.04, .04);    //offset A0
        myfit[2][1]->SetParLimits(1, -.04, .04);    //Dipole AD
        
        myfit[2][2]->SetParLimits(0, -.04, .04);    //offset A0
        myfit[2][2]->SetParLimits(1, -.04, .04);    //Dipole AD
    }
    
    else if(OFFSET_PLUS_DIPOLE_PLUS_QUADRUPOLE){
    
        myfit[2][0]->SetParLimits(0, -.04, .04);    //offset A0
        myfit[2][0]->SetParLimits(1, -.04, .04);    //Dipole AD
        myfit[2][0]->SetParLimits(2, -.04, .04);    //quadrupole AQ
        
        myfit[2][1]->SetParLimits(0, -.04, .04);    //offset A0
        myfit[2][1]->SetParLimits(1, -.04, .04);    //Dipole AD
        myfit[2][1]->SetParLimits(2, -.04, .04);    //quadrupole AQ
        
        myfit[2][2]->SetParLimits(0, -.04, .04);    //offset A0
        myfit[2][2]->SetParLimits(1, -.04, .04);    //Dipole AD
        myfit[2][2]->SetParLimits(2, -.04, .04);    //quadrupole AQ
    }
    
    else if(FULL_FIT_WITH_AWAY_SIDE_GAUSS){ 
        
        //                         A0      AS Gauss amp  AS Gauss width     AQ    A1   sPhi sEta  
        myfit[2][0]->SetParameters(0.0,   .04,            1.0,             .011,  .11,  .3, .3);
        
        myfit[2][0]->SetParLimits(0, -0.01, 0.01);    //offset A0
        myfit[2][0]->SetParLimits(1, 0.01,  .05);    //Away-side gaussian amp
        myfit[2][0]->SetParLimits(2, .6,  1.2);    //Away-side gaussian width
        myfit[2][0]->SetParLimits(3, 0.009,  .012);    //v2     AQ
        myfit[2][0]->SetParLimits(4, 0.08, .13);    //A1 jet amp
        myfit[2][0]->SetParLimits(5, 0.2, .4);     //phi width
        myfit[2][0]->SetParLimits(6, 0.2, .4);    //eta width
        
        //                         A0      AS Gauss amp  AS Gauss width     AQ    A1   sPhi sEta  
        myfit[2][1]->SetParameters(-.012,   .018,            1.5,          .007,  .04,  .45, .55);
        
        myfit[2][1]->SetParLimits(0, -.015, -.010);    //offset A0
        myfit[2][1]->SetParLimits(1, 0.015,  0.02);    //Away-side gaussian amp
        myfit[2][1]->SetParLimits(2, 1.2,  1.7);    //Away-side gaussian width
        myfit[2][1]->SetParLimits(3, 0.006,  .0082);    //v2     AQ
        myfit[2][1]->SetParLimits(4, 0.035, .05);    //A1 jet amp
        myfit[2][1]->SetParLimits(5, 0.2, .6);     //phi width
        myfit[2][1]->SetParLimits(6, 0.4, .6);    //eta width
        
        //                         A0      AS Gauss amp  AS Gauss width     AQ    A1   sPhi sEta  
        myfit[2][2]->SetParameters(-.01,   .02,            1.0,             .007,  .04,  .43, .6);
        
        myfit[2][2]->SetParLimits(0, -.02, .01);    //offset A0
        myfit[2][2]->SetParLimits(1, 0.0,  0.04);    //Away-side gaussian amp
        myfit[2][2]->SetParLimits(2, .5,  1.3);    //Away-side gaussian width
        myfit[2][2]->SetParLimits(3, 0.0,  .02);    //v2     AQ
        myfit[2][2]->SetParLimits(4, 0.0, .04);    //A1 jet amp
        myfit[2][2]->SetParLimits(5, 0.3, .8);     //phi width
        myfit[2][2]->SetParLimits(6, 0.3, .8);    //eta width 
    }
    
    else if(FULL_FIT_NO_DIPOLE){
    
     //                             A0     AQ    A1   sPhi sEta  
        myfit[2][0]->SetParameters(-0.01, .02, .03,  .2, .2);
    
        myfit[2][0]->SetParLimits(0, -.05, 0.0);    //offset A0
        myfit[2][0]->SetParLimits(1, 0.01,  .07);    //v2     AQ
        myfit[2][0]->SetParLimits(2, 0.02, .1);    //A1 jet amp
        myfit[2][0]->SetParLimits(3, 0.0, 1.0);     //phi width
        myfit[2][0]->SetParLimits(4, 0.0, 1.0);    //eta width
        
         //                         A0    AQ    A1   sPhi sEta  
        myfit[2][1]->SetParameters(-0.01, .02, .01,  .2, .2);
        
        myfit[2][1]->SetParLimits(0, -0.02, 0.0);    //offset A0
        myfit[2][1]->SetParLimits(1, 0.01,  .04);    //v2     AQ
        myfit[2][1]->SetParLimits(2, 0.02, .04);    //A1 jet amp
        myfit[2][1]->SetParLimits(3, 0.1, 1.2);     //phi width
        myfit[2][1]->SetParLimits(4, 0.1, 1.2);    //eta width
        
        //                         A0        AQ    A1   sPhi sEta  
        myfit[2][2]->SetParameters(-.0007, .003, .02,  .32, .37);
        
        myfit[2][2]->SetParLimits(0, -.003, 0.0);    //offset A0
        myfit[2][2]->SetParLimits(1, 0.0,  .015);    //v2     AQ
        myfit[2][2]->SetParLimits(2, 0.0, .05);    //A1 jet amp
        myfit[2][2]->SetParLimits(3, 0.1, 1.2);     //phi width
        myfit[2][2]->SetParLimits(4, 0.1, 1.2);    //eta width 
    }
    
    else if(FULL_FIT_WITH_DIPOLE){
    
       //                             A0   AD   AQ    A1   sPhi sEta  
        myfit[2][0]->SetParameters(-0.001, -.01, .02, .03,  .2, .4);
    
        myfit[2][0]->SetParLimits(0, -.01, 0.0);    //offset A0
        myfit[2][0]->SetParLimits(1, -.05,  0.0);    //    AD
        myfit[2][0]->SetParLimits(2, 0.0,  .05);    //v2     AQ
        myfit[2][0]->SetParLimits(3, 0.01, .09);    //A1 jet amp
        myfit[2][0]->SetParLimits(4, 0.0, .8);     //phi width
        myfit[2][0]->SetParLimits(5, 0.0, 1.0);    //eta width
        
         //                         A0    AD    AQ    A1   sPhi sEta  
        myfit[2][1]->SetParLimits(0, -.02, 0.0);    //offset A0
        //myfit[2][1]->SetParLimits(1, -.03, -0.001);    //v1     AD
        //myfit[2][1]->SetParLimits(2, 0.009,  .015);    //v2     AQ
        //myfit[2][1]->FixParameter(1, -.0177424); //AD     --smaller eta range
        //myfit[2][1]->FixParameter(2, .00605651); //v2 AQ  --smaller eta range
		
		myfit[2][1]->FixParameter(1, -.0165570); //AD
        myfit[2][1]->FixParameter(2, .00647402); //v2 AQ
        myfit[2][1]->SetParLimits(3, 0.02, .08);    //A1 jet amp
        myfit[2][1]->SetParLimits(4, 0.2, 1.0);     //phi width
        myfit[2][1]->SetParLimits(5, .5, 4.3);    //eta width
        
        //                         A0       AD     AQ    A1   sPhi sEta  
        myfit[2][2]->SetParameters(-.0007, -.001, .003, .02,  .32, .37);
        
        myfit[2][2]->SetParLimits(0, -.01, 0.0);    //offset A0
        myfit[2][2]->SetParLimits(1, -.01,  0.0);    //    AD
        myfit[2][2]->SetParLimits(2, 0.0,  .015);    //v2     AQ
        myfit[2][2]->SetParLimits(3, 0.0, .05);    //A1 jet amp
        myfit[2][2]->SetParLimits(4, 0.1, 1.2);     //phi width
        myfit[2][2]->SetParLimits(5, 0.1, 1.2);    //eta width 
        
        
    }
    
    else if(FULL_FIT_WITH_DIPOLE_AND_ETA_GAUSSIAN){ //central fitting
    
    /* double A0       = -.04;
    double AD       = -.03;
    double AQ       = .015;
    double A1       = .13;
    double A1SigPhi = .2;
    double A1SigEta = .2;
    double etaPeak  = .04;
    double etaWidth = .8; */
    
        myfit[2][0]->SetParLimits(0, -.02, .01);    //offset A0
        myfit[2][0]->SetParLimits(1, 0.002,  0.004);    //v2     AQ
        myfit[2][0]->SetParLimits(2, 0.0, .01);    //A1 jet amp
        myfit[2][0]->SetParLimits(3, 0.2, .6);     //phi width
        myfit[2][0]->SetParLimits(4, 1.5, 2.5);    //eta width 
        myfit[2][0]->SetParLimits(5, -.05,  0.0);    //As dipole
        myfit[2][0]->SetParLimits(6, 0.3, .6);     //AS Eta Width
        
        myfit[2][1]->SetParLimits(0, -.02, .01);    //offset A0
        myfit[2][1]->SetParLimits(1, 0.002,  0.004);    //v2     AQ
        myfit[2][1]->SetParLimits(2, 0.0, .01);    //A1 jet amp
        myfit[2][1]->SetParLimits(3, 0.2, .6);     //phi width
        myfit[2][1]->SetParLimits(4, 1.5, 2.5);    //eta width 
        myfit[2][1]->SetParLimits(5, -.05,  0.0);    //As dipole
        myfit[2][1]->SetParLimits(6, 0.3, .6);     //AS Eta Width
        
        //                         A0      AD     AQ    A1   sPhi sEta  etaGPeak  etaGwidth
        //myfit[2][2]->SetParameters(-.001, -.01, .01, .025, .3, .3, .02, .03);
        
         myfit[2][2]->SetParameters(-.001, 0.0,  .04,     .75,    1.32,     -.0145       ,   1.06      );
        
        myfit[2][2]->SetParLimits(0, -.015, 0.0);    //offset  
        myfit[2][2]->SetParLimits(1, 0.0, .003);    //quadrupole  
        myfit[2][2]->SetParLimits(2, 0.008,  .19);   //jet amp
        myfit[2][2]->SetParLimits(3, 0.3, .9);    //jet phi width
        myfit[2][2]->SetParLimits(4, .5, 1.9);    //jet eta width 
        myfit[2][2]->SetParLimits(5, -0.03, 0.0);    //AS dipole
        myfit[2][2]->SetParLimits(6, 0.05, 3.5);    //AS eta width
        
    }
    
    else if(DOUBLE_2D_GAUSS){
    
    /* double A0       = -.04;
    double AD       = -.03;
    double AQ       = .015;
    double A1       = .13;
    double A1SigPhi = .2;
    double A1SigEta = .2;
    double etaPeak  = .04;
    double etaWidth = .8; */
    
        myfit[2][0]->SetParameters(-.001, -.005, .005, .09, .25, .09, .25, .25, .01, .9);
    
       //myfit[2][0]->SetParLimits(0, -.02, 0.0);    //offset A0
        myfit[2][0]->SetParLimits(1, -.009,  0.0);    //v1     AD
        myfit[2][0]->SetParLimits(2, 0.0,  .05);    //v2     AQ
        myfit[2][0]->SetParLimits(2, 0.0,  .05);
        myfit[2][0]->SetParLimits(3, 0.03, .07);    //A1 jet amp 1
        myfit[2][0]->SetParLimits(4, 0.1, .4);    //jet eta width 1
        myfit[2][0]->SetParLimits(5, 0.03, .07);    //A1 jet amp 2
        myfit[2][0]->SetParLimits(6, 0.1, .4);    //jet eta width 2
        myfit[2][0]->SetParLimits(7, 0.1, .4);     //phi width 
        //myfit[2][0]->SetParLimits(8, 0.0, .02);     //eta_Gauss_peak
        myfit[2][0]->SetParLimits(9, 1.0, 2.5);    //eta_Gauss width
        
        myfit[2][1]->SetParameters(-.01, -.01, .03, .06, .3, .06, .3, .3, .01, 1.2);
       
        //myfit[2][1]->SetParLimits(0, -.02, 0.0);    //offset A0
        myfit[2][1]->SetParLimits(1, -.009,  0.0);    //v1     AD
        myfit[2][1]->SetParLimits(2, 0.0,  .05);    //v2     AQ
        myfit[2][1]->SetParLimits(3, 0.03, .07);    //A1 jet amp 1
        myfit[2][1]->SetParLimits(4, 0.1, .4);    //jet eta width 1
        myfit[2][1]->SetParLimits(5, 0.03, .07);    //A1 jet amp 2
        myfit[2][1]->SetParLimits(6, 0.1, .4);    //jet eta width 2
        myfit[2][1]->SetParLimits(7, 0.1, .4);     //phi width 
        //myfit[2][1]->SetParLimits(8, 0.0, .02);     //eta_Gauss_peak
        myfit[2][1]->SetParLimits(9, 1.0, 2.5);    //eta_Gauss width
        
        
        myfit[2][2]->SetParameters(-.001, -.005, .01, .025, .3, .025, .3, .3, .02, .03);
        
        //myfit[2][2]->SetParLimits(0, -.02, 0.0);    //offset A0
        myfit[2][2]->SetParLimits(1, -.009,  0.0);    //v1     AD
        myfit[2][2]->SetParLimits(2, 0.0,  .05);    //v2     AQ
        myfit[2][2]->SetParLimits(3, 0.03, .07);    //A1 jet amp 1
        myfit[2][2]->SetParLimits(4, 0.1, .4);    //jet eta width 1
        myfit[2][2]->SetParLimits(5, 0.03, .07);    //A1 jet amp 2
        myfit[2][2]->SetParLimits(6, 0.1, .4);    //jet eta width 2
        myfit[2][2]->SetParLimits(7, 0.1, .4);     //phi width 
        //myfit[2][2]->SetParLimits(8, 0.0, .02);     //eta_Gauss_peak
        myfit[2][2]->SetParLimits(9, 1.0, 2.5);    //eta_Gauss width
    }
  
    if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS){  //peripheral 
  
        myfit[2][0]->SetParameters(-0.011, .004, .11, .3, .3, .028, .52);
        
        myfit[2][0]->SetParLimits(0, -.015, -0.002);    //offset A0
        myfit[2][0]->SetParLimits(1, 0.002,  0.005);    //v2     AQ
        myfit[2][0]->SetParLimits(2, 0.01, .14);    //A1 jet amp
        myfit[2][0]->SetParLimits(3, 0.2, .4);     //phi width
        myfit[2][0]->SetParLimits(4, 0.24, .4);    //eta width 
        myfit[2][0]->SetParLimits(5, .01,  0.032);    //AS amp
        myfit[2][0]->SetParLimits(6, 0.48, .7);     //AS Phi Width
        
        myfit[2][1]->SetParLimits(0, -.02, .01);    //offset A0
        myfit[2][1]->SetParLimits(1, 0.002,  0.004);    //v2     AQ
        myfit[2][1]->SetParLimits(2, 0.0, .01);    //A1 jet amp
        myfit[2][1]->SetParLimits(3, 0.2, .6);     //phi width
        myfit[2][1]->SetParLimits(4, 1.5, 2.5);    //eta width 
        myfit[2][1]->SetParLimits(5, -.05,  0.0);    //AS amp
        myfit[2][1]->SetParLimits(6, 0.3, .6);     //AS Phi Width
        
        myfit[2][2]->SetParLimits(0, -.02, .01);    //offset A0
        myfit[2][2]->SetParLimits(1, 0.002,  0.004);    //v2     AQ
        myfit[2][2]->SetParLimits(2, 0.0, .01);    //A1 jet amp
        myfit[2][2]->SetParLimits(3, 0.2, .6);     //phi width
        myfit[2][2]->SetParLimits(4, 1.5, 2.5);    //eta width 
        myfit[2][2]->SetParLimits(5, -.05,  0.0);    //AS amp
        myfit[2][2]->SetParLimits(6, 0.3, .6);     //AS Phi Width
        
    }
  
    cout << endl << "BEGIN FITTING CORRELATIONS HERE" << endl << endl;
    
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        for(int centBin = 0; centBin < 3; centBin++){
    
            absCorrHist[ptBin][centBin]->Fit(myfit[ptBin][centBin], "R0");
        
            //cout << myfit[ptBin][centBin]->GetChisquare() << endl;
        
        }
    }
    
    
    cout << endl;
    //CALCULATE CHI-SQUARE HERE BY HAND
    
    double chiSquare = 0;
    double chiSquareFull = 0;
    double tempData = 0;
    double tempModel = 0;
    double tempError = 0;
    double tempFinal = 0;
    double etaBinCenter = 0;
    double phiBinCenter = 0;
    
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        for(int centBin = 0; centBin < 3; centBin++){
        
            chiSquareFull = 0;
            chiSquare = 0;
            tempData = 0;
            tempModel = 0;
            tempError = 0;
            tempFinal = 0;
			
			double quarterEtaBin     = (4.0/NUM_SYMM_ETA_BINS)*.25;
			double quarterPhiBinLow  = (TMath::TwoPi()/NUM_SYMM_PHI_BINS)*.25;
			double quarterPhiBinHigh = TMath::Pi() - quarterPhiBinLow ;
	
			for(int etaBin = 1; etaBin < NUM_ABS_VAL_ETA_BINS+1; etaBin++){
                for(int phiBin = 1; phiBin < NUM_ABS_VAL_PHI_BINS+1; phiBin++){
                
                    etaBinCenter = absCorrHist[ptBin][centBin]->GetXaxis()->GetBinCenter(etaBin);
                    phiBinCenter = absCorrHist[ptBin][centBin]->GetYaxis()->GetBinCenter(phiBin);
                    //if(TMath::Abs(etaBinCenter) < .005) { etaBinCenter = 0.0; }
                    //if(TMath::Abs(phiBinCenter) < .005) { phiBinCenter = 0.0; }
                    
                    tempData = absCorrHist[ptBin][centBin]->GetBinContent(etaBin, phiBin);
                    
                    //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Eval(myfit[k][i]);
                    
                    if((etaBin == 1 && phiBin == 1)) { tempModel = myfit[ptBin][centBin]->Eval(quarterEtaBin, quarterPhiBinLow, 0, 0); }
                    else if((etaBin == 1 && phiBin == 9)) { tempModel = myfit[ptBin][centBin]->Eval(quarterEtaBin, quarterPhiBinHigh, 0, 0); }
                    else if((etaBin == 1 && phiBin != 1)) { tempModel = myfit[ptBin][centBin]->Eval(quarterEtaBin, phiBinCenter, 0, 0); }
                    else if((etaBin != 1 && phiBin == 1)) { tempModel = myfit[ptBin][centBin]->Eval(etaBinCenter, quarterPhiBinLow, 0, 0); }
                    else if((etaBin != 1 && phiBin == 9)) { tempModel = myfit[ptBin][centBin]->Eval(etaBinCenter, quarterPhiBinHigh, 0, 0); }
                    else tempModel = myfit[ptBin][centBin]->Eval(etaBinCenter, phiBinCenter, 0, 0);
                    
                    tempError = absCorrHist[ptBin][centBin]->GetBinError(etaBin, phiBin);
                    
                    tempFinal = (tempData - tempModel)/tempError;
                    
                    //cout << tempData << "  " << tempModel << " " << tempError << endl;
                    
                    chiSquare = (tempFinal*tempFinal);
                    
                    chiSquareHist[ptBin][centBin]->SetBinContent(etaBin, phiBin, chiSquare);
                
                    chiSquareFull = chiSquareFull + chiSquare;
                }
            }
            
			//#bins = 45 for 13 and 16
			//"    " = 28 for 9 and 12
			//#parameters for full-fit+dipole = 6
			
            cout << "root chi square: " << myfit[ptBin][centBin]->GetChisquare() << endl;
            if(FULL_FIT_WITH_DIPOLE) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(29.0) << endl; }
            else if(FULL_FIT_NO_DIPOLE) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(23.0) << endl; }
            else if(OFFSET_ONLY) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(27.0) << endl; }
            else if(OFFSET_PLUS_DIPOLE) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(26.0) << endl; }
            else if(OFFSET_PLUS_DIPOLE_PLUS_QUADRUPOLE) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(25.0) << endl; }
            else if(FULL_FIT_WITH_AWAY_SIDE_GAUSS) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(21.0) << endl; }
            else if(FULL_FIT_WITH_DIPOLE_AND_ETA_GAUSSIAN) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(28.0) << endl; }
            else if(DOUBLE_2D_GAUSS) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(18.0) << endl; }
            else if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { cout << "root chi square p.d.f. : " << myfit[ptBin][centBin]->GetChisquare()/(28.0) << endl; }
            
            cout << "hand-calculated Chi-Square: " << chiSquareFull << endl;
            if(FULL_FIT_WITH_DIPOLE) { cout << "hand-calculated Chi-Square p.d.f.: " << (chiSquareFull)/(29.0) << endl; }
            else if(FULL_FIT_NO_DIPOLE) { cout << "hand-calculated Chi-Square p.d.f.: " << (chiSquareFull)/(23.0) << endl; }
            else if(OFFSET_ONLY) { cout << "root chi square p.d.f. : " << (chiSquareFull)/(27.0) << endl; }
            else if(OFFSET_PLUS_DIPOLE) { cout << "root chi square p.d.f. : " << (chiSquareFull)/(26.0) << endl; }
            else if(OFFSET_PLUS_DIPOLE_PLUS_QUADRUPOLE) { cout << "root chi square p.d.f. : " << (chiSquareFull)/(25.0) << endl; }
            else if(FULL_FIT_WITH_AWAY_SIDE_GAUSS) { cout << "hand-calculated Chi-Square p.d.f.: " << (chiSquareFull)/(21.0) << endl; }
            else if(FULL_FIT_WITH_DIPOLE_AND_ETA_GAUSSIAN) { cout << "hand-calculated Chi-Square p.d.f.: " << (chiSquareFull)/(28.0) << endl; }
            else if(DOUBLE_2D_GAUSS) { cout << "hand-calculated Chi-Square p.d.f.: " << (chiSquareFull)/(18.0) << endl; }
            else if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { cout << "hand-calculated Chi-Square p.d.f.: " << (chiSquareFull)/(28.0) << endl; }
            cout << endl;
        }
    }
    
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        for(int centBin = 0; centBin < 3; centBin++){ 
        
            absCorrHistFit[ptBin][centBin]->Eval(myfit[ptBin][centBin]);
            absCorrHistRes[ptBin][centBin]->Add(absCorrHist[ptBin][centBin], absCorrHistFit[ptBin][centBin], 1, -1);
            
            absCorrHistRes[ptBin][centBin]->SetMaximum(absCorrHist[ptBin][centBin]->GetMaximum());
        
            absCorrHistFit[ptBin][centBin]->Write();
            absCorrHistRes[ptBin][centBin]->Write();
        
        
        }
    }    
    
    TString outputFitsLabel = "outputFits_";
    
    for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        
        TCanvas * can = new TCanvas("", "", 4400, 3300);
        can->Divide(4,3);
        int padNum = 1;
  
    
        for(int centBin = 0; centBin < 3; centBin++){ 
        
            can->cd(padNum);
            absCorrHist[ptBin][centBin]->Draw("SURF1");
            
            padNum++;
            can->cd(padNum);
            absCorrHistFit[ptBin][centBin]->Draw("SURF1");
            
            padNum++;
            can->cd(padNum);
            absCorrHistRes[ptBin][centBin]->Draw("SURF1");
            
            padNum++;
            can->cd(padNum);
            chiSquareHist[ptBin][centBin]->Draw("SURF1");
            
            padNum++;
        }
        
        str1 = outputFitsLabel + ptBinLabel[ptBin] + pngLabel;
        can->SaveAs(str1);
    }    
  
    
    
    TString fitParametersLabel = "FitParameters_";
    
    
    if(FULL_FIT_WITH_DIPOLE){
        
        TCanvas * canFitParameters = new TCanvas("", "", 3300, 2200);
        canFitParameters->Divide(3,2);
        
        for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        
            str1 = offsetLabel + ptBinLabel[ptBin];
            fitParameterOffset[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            str1 = dipoleLabel + ptBinLabel[ptBin];
            fitParameterDipole[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            str1 = quadLabel + ptBinLabel[ptBin];
            fitParameterQuad[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            str1 = jetAmpLabel + ptBinLabel[ptBin];
            fitParameterJetAmp[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            str1 = sigEtaLabel + ptBinLabel[ptBin];
            fitParameterSigEta[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            str1 = sigPhiLabel + ptBinLabel[ptBin];
            fitParameterSigPhi[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            
            fitParameterOffset[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterOffset[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterOffset[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            fitParameterDipole[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterDipole[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%"); 
            fitParameterDipole[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%"); 
            fitParameterQuad[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterQuad[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterQuad[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");            
            fitParameterJetAmp[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterJetAmp[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterJetAmp[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            fitParameterSigEta[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterSigEta[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterSigEta[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            fitParameterSigPhi[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterSigPhi[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterSigPhi[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            
            fitParameterOffset[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(0));
            fitParameterOffset[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(0));
            fitParameterOffset[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(0));
            fitParameterOffset[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(0));
            fitParameterOffset[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(0));
            fitParameterOffset[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(0));
            
            fitParameterDipole[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(1));
            fitParameterDipole[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(1));
            fitParameterDipole[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(1));
            fitParameterDipole[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(1));
            fitParameterDipole[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(1));
            fitParameterDipole[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(1));
            
            fitParameterQuad[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(2));
            fitParameterQuad[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(2));
            fitParameterQuad[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(2));
            fitParameterQuad[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(2));
            fitParameterQuad[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(2));
            fitParameterQuad[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(2));
            
            fitParameterJetAmp[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(3));
            fitParameterJetAmp[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(3));
            fitParameterJetAmp[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(3));
            fitParameterJetAmp[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(3));
            fitParameterJetAmp[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(3));
            fitParameterJetAmp[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(3));
            
            fitParameterSigEta[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(5));
            fitParameterSigEta[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(5));
            fitParameterSigEta[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(5));
            fitParameterSigEta[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(5));
            fitParameterSigEta[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(5));
            fitParameterSigEta[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(5));
            
            fitParameterSigPhi[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(4));
            fitParameterSigPhi[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(4));
            fitParameterSigPhi[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(4));
            fitParameterSigPhi[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(4));
            fitParameterSigPhi[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(4));
            fitParameterSigPhi[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(4));
            
            canFitParameters->cd(1);
            fitParameterOffset[ptBin]->SetMarkerStyle(20);
            fitParameterOffset[ptBin]->Draw("P");
            canFitParameters->cd(2);
            fitParameterDipole[ptBin]->SetMarkerStyle(20);
            fitParameterDipole[ptBin]->Draw("p");
            canFitParameters->cd(3);
            fitParameterQuad[ptBin]->SetMarkerStyle(20);
            fitParameterQuad[ptBin]->Draw("p");
            canFitParameters->cd(4);
            fitParameterJetAmp[ptBin]->SetMarkerStyle(20);
            fitParameterJetAmp[ptBin]->Draw("p");
            canFitParameters->cd(5);
            fitParameterSigEta[ptBin]->SetMarkerStyle(20);
            fitParameterSigEta[ptBin]->Draw("p");
            canFitParameters->cd(6);
            fitParameterSigPhi[ptBin]->SetMarkerStyle(20);
            fitParameterSigPhi[ptBin]->Draw("p");
        
            str1 = fitParametersLabel + ptBinLabel[ptBin] + pngLabel;
        
            canFitParameters->SaveAs(str1);
        }
    }
  
    if(FULL_FIT_WITH_AWAY_SIDE_GAUSS){
    
        TCanvas * canFitParameters = new TCanvas("", "", 4400, 2200);
        canFitParameters->Divide(4,2);
  
  
        for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
       
            offsetLabel = offsetLabel + ptBinLabel[ptBin];
            dipoleLabel = dipoleLabel + ptBinLabel[ptBin];
            quadLabel   = quadLabel + ptBinLabel[ptBin];
            jetAmpLabel = jetAmpLabel + ptBinLabel[ptBin];
            sigEtaLabel = sigEtaLabel + ptBinLabel[ptBin];
            sigPhiLabel = sigPhiLabel + ptBinLabel[ptBin];
    
            fitParameterOffset[ptBin] = new TH1D(offsetLabel, offsetLabel, 3, 0, 3);
            fitParameterASGaussWidth[ptBin] = new TH1D(ASGaussWidthLabel, ASGaussWidthLabel, 3, 0, 3);
            fitParameterASGaussAmp[ptBin] = new TH1D(ASGaussAmpLabel, ASGaussAmpLabel, 3, 0, 3);
            fitParameterQuad[ptBin]   = new TH1D(quadLabel, quadLabel, 3, 0, 3);
            fitParameterJetAmp[ptBin] = new TH1D(jetAmpLabel, jetAmpLabel, 3, 0, 3);
            fitParameterSigEta[ptBin] = new TH1D(sigEtaLabel, sigEtaLabel, 3, 0, 3);
            fitParameterSigPhi[ptBin] = new TH1D(sigPhiLabel, sigPhiLabel, 3, 0, 3);
            
            fitParameterOffset[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterOffset[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterOffset[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            fitParameterASGaussWidth[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterASGaussWidth[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%"); 
            fitParameterASGaussWidth[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%"); 
            fitParameterASGaussAmp[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterASGaussAmp[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%"); 
            fitParameterASGaussAmp[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%"); 
            fitParameterQuad[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterQuad[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterQuad[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");            
            fitParameterJetAmp[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterJetAmp[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterJetAmp[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            fitParameterSigEta[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterSigEta[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterSigEta[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            fitParameterSigPhi[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            fitParameterSigPhi[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            fitParameterSigPhi[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
            
            fitParameterOffset[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(0));
            fitParameterOffset[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(0));
            fitParameterOffset[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(0));
            fitParameterOffset[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(0));
            fitParameterOffset[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(0));
            fitParameterOffset[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(0));
            
            fitParameterASGaussAmp[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(1));
            fitParameterASGaussAmp[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(1));
            fitParameterASGaussAmp[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(1));
            fitParameterASGaussAmp[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(1));
            fitParameterASGaussAmp[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(1));
            fitParameterASGaussAmp[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(1));
            
            fitParameterASGaussWidth[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(2));
            fitParameterASGaussWidth[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(2));
            fitParameterASGaussWidth[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(2));
            fitParameterASGaussWidth[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(2));
            fitParameterASGaussWidth[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(2));
            fitParameterASGaussWidth[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(2));
            
            fitParameterQuad[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(3));
            fitParameterQuad[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(3));
            fitParameterQuad[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(3));
            fitParameterQuad[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(3));
            fitParameterQuad[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(3));
            fitParameterQuad[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(3));
            
            fitParameterJetAmp[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(4));
            fitParameterJetAmp[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(4));
            fitParameterJetAmp[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(4));
            fitParameterJetAmp[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(4));
            fitParameterJetAmp[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(4));
            fitParameterJetAmp[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(4));
            
            fitParameterSigEta[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(6));
            fitParameterSigEta[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(6));
            fitParameterSigEta[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(6));
            fitParameterSigEta[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(6));
            fitParameterSigEta[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(6));
            fitParameterSigEta[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(6));
            
            fitParameterSigPhi[ptBin]->SetBinContent(1, myfit[ptBin][0]->GetParameter(5));
            fitParameterSigPhi[ptBin]->SetBinContent(2, myfit[ptBin][1]->GetParameter(5));
            fitParameterSigPhi[ptBin]->SetBinContent(3, myfit[ptBin][2]->GetParameter(5));
            fitParameterSigPhi[ptBin]->SetBinError(1, myfit[ptBin][0]->GetParError(5));
            fitParameterSigPhi[ptBin]->SetBinError(2, myfit[ptBin][1]->GetParError(5));
            fitParameterSigPhi[ptBin]->SetBinError(3, myfit[ptBin][2]->GetParError(5));
            
            canFitParameters->cd(1);
            fitParameterOffset[ptBin]->Draw();
            canFitParameters->cd(2);
            fitParameterASGaussAmp[ptBin]->Draw();
            canFitParameters->cd(3);
            fitParameterASGaussWidth[ptBin]->Draw();
            canFitParameters->cd(4);
            fitParameterQuad[ptBin]->Draw();
            canFitParameters->cd(5);
            fitParameterJetAmp[ptBin]->Draw();
            canFitParameters->cd(6);
            fitParameterSigEta[ptBin]->Draw();
            canFitParameters->cd(7);
            fitParameterSigPhi[ptBin]->Draw();
            
            
            canFitParameters->SaveAs("outputfitParameters.png");
        
        }
    }
    
  
  
    
    
    double jetVolumeYield;
    double jetVolumeYieldError;
   
   
    TString jetVolumeYieldLabel = "Jet per-trigger yield ";
    TCanvas * jetYieldCanvas = new TCanvas("", "", 3300, 1100);
    jetYieldCanvas->Divide(3,1);
        
   
    if(FULL_FIT_WITH_DIPOLE) { 
        for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        
            str1 = jetVolumeYieldLabel + ptBinLabel[ptBin];
            jetVolumePlot[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
        
            jetYieldCanvas->cd(ptBin-1);
        
            for(int centBin = 0; centBin < 3; centBin++){
            
                
                jetVolumeYield = (myfit[ptBin][centBin]->GetParameter(3)*myfit[ptBin][centBin]->GetParameter(4)*myfit[ptBin][centBin]->GetParameter(5)*nBarPrime[centBin])/2.0;
                jetVolumePlot[ptBin]->SetBinContent(centBin+1, jetVolumeYield);
                
                jetVolumeYieldError = getJetYieldError(myfit[ptBin][centBin]->GetParameter(3), myfit[ptBin][centBin]->GetParameter(4), myfit[ptBin][centBin]->GetParameter(5),
                                                       myfit[ptBin][centBin]->GetParError(3), myfit[ptBin][centBin]->GetParError(4), myfit[ptBin][centBin]->GetParError(5));
                
                cout << jetVolumeYieldError << endl;
                
                jetVolumePlot[ptBin]->SetBinError(centBin+1, jetVolumeYieldError);
                
                cout << "# jet tracks per D0 in pt-bin: " << ptBin << " and Centrality Bin: " << "  = " << jetVolumeYield << endl;
            
            }  

            jetVolumePlot[ptBin]->SetMarkerStyle(20); 
            jetVolumePlot[ptBin]->Draw("p");             
        }   
        
        str1 = jetVolumeYieldLabel + pngLabel;
        jetYieldCanvas->SaveAs(str1);
        
    }
    //else if(FULL_FIT_WITH_AWAY_SIDE_GAUSS) { cout << "Using option: FULL FIT USING AS GAUSSIAN (28 - 7) DF" << endl; }
    //else if(FULL_FIT_WITH_DIPOLE_AND_ETA_GAUSSIAN) { cout << "Using option: FULL FIT USING DIPOLE & ETA_DEPENDENT GAUSSIAN (28 - 8) DF" << endl; }
    //else if(DOUBLE_2D_GAUSS) { cout << "Using option: FULL FIT COMPLICATED JET SHAPE (28 - 10) DF" << endl; }
  
    inputFile->Close();
    outputFile->Close();
    
}

int getPhiBin(int bin){

    /* if(bin == 2 || bin == 3) {return 4;}
    if(bin == 4 || bin == 5) {return 5;}
    if(bin == 6 || bin == 7) {return 6;}
    if(bin == 8 || bin == 9) {return 7;}
    if(bin == 10 || bin == 11) {return 8;} */
    
    if(bin == 1) {return 3;}
    if(bin == 2) {return 4;}
    if(bin == 3) {return 5;}
    if(bin == 4) {return 6;}
    if(bin == 5) {return 7;}
    if(bin == 6) {return 8;}
    if(bin == 7) {return 9;}
	
    
    else return 0;
}

int getEtaBin(int bin){

    /* if(bin == 2 || bin == 3) {return 6;}
    if(bin == 4 || bin == 5) {return 7;}
    if(bin == 6 || bin == 7) {return 8;} */
    
    if(bin == 1) {return 7;}
    if(bin == 2) {return 8;}
    if(bin == 3) {return 9;}
    if(bin == 4) {return 10;}
	if(bin == 5) {return 11;}
    if(bin == 6) {return 12;}
    
   
	
    
    
    else return 0;
}

double getJetYieldError(double A, double sigEta, double sigPhi, double eA, double eSigEta, double eSigPhi){

    double arg = ((eA/A)*(eA/A)) + ((eSigEta/sigEta)*(eSigEta/sigEta)) + ((eSigPhi/sigPhi)*(eSigPhi/sigPhi));
    double product = TMath::Abs((A*sigEta*sigPhi));
    
    cout << "prduct value: " << product << endl;
    cout << "arg value: " << arg << endl;
    cout << "sqrt of arg: " << TMath::Sqrt(arg) << endl;
    
    return TMath::Abs((A*sigEta*sigPhi))*TMath::Sqrt(arg);
    
}
    
 
    
    /* myfit[0][0]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[0][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[0][0]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[0][0]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[0][0]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[0][0]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[0][1]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[0][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[0][1]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[0][1]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[0][1]->SetParLimits(4, 0.1, 2);     //phi width
    myfit[0][1]->SetParLimits(5, 0.1, 2);    //eta width
	
	myfit[0][2]->SetParLimits(0, -.1, .01);    //offset A0
    myfit[0][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[0][2]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[0][2]->SetParLimits(3, 0.0, .04);    //A1 jet amp
    myfit[0][2]->SetParLimits(4, 0.1, 1.2);     //phi width
    myfit[0][2]->SetParLimits(5, 0.1, 1.2);    //eta width
	
	myfit[1][0]->SetParLimits(0, -.08, .01);    //offset A0
    myfit[1][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[1][0]->SetParLimits(2, 0.0,  .06);    //v2     AQ
    myfit[1][0]->SetParLimits(3, 0.0, .1);    //A1 jet amp
    myfit[1][0]->SetParLimits(4, 0.6, 1.3);     //phi width
    myfit[1][0]->SetParLimits(5, 0.5, 1.2);    //eta width
	
	myfit[1][1]->SetParLimits(0, -.07, .01);    //offset A0
    myfit[1][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[1][1]->SetParLimits(2, 0.0,  .06);    //v2     AQ
    myfit[1][1]->SetParLimits(3, 0.0, .04);    //A1 jet amp
    myfit[1][1]->SetParLimits(4, 0.8, 1.2);     //phi width
    myfit[1][1]->SetParLimits(5, 1.1, 1.6);    //eta width
	
	myfit[1][2]->SetParLimits(0, -.04, .01);    //offset A0
    myfit[1][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[1][2]->SetParLimits(2, 0.0,  .04);    //v2     AQ
    myfit[1][2]->SetParLimits(3, 0.0, .01);    //A1 jet amp
    myfit[1][2]->SetParLimits(4, 0.1, 1.2);     //phi width
    myfit[1][2]->SetParLimits(5, 0.1, 1.0);    //eta width */
    ////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////
    /* myfit[3][0]->SetParLimits(0, -.07, .01);    //offset A0
    myfit[3][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[3][0]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[3][0]->SetParLimits(3, 0.0, .25);    //A1 jet amp
    myfit[3][0]->SetParLimits(4, 0.1, .8);     //phi width
    myfit[3][0]->SetParLimits(5, 0.4, .7);    //eta width
	
	myfit[3][1]->SetParLimits(0, -.02, .01);    //offset A0
    myfit[3][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[3][1]->SetParLimits(2, 0.0,  .06);    //v2     AQ
    myfit[3][1]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[3][1]->SetParLimits(4, 0.1, .8);     //phi width
    myfit[3][1]->SetParLimits(5, 0.1, 1.0);    //eta width
	
	myfit[3][2]->SetParLimits(0, -.03, .01);    //offset A0
    myfit[3][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[3][2]->SetParLimits(2, 0.0,  .03);    //v2     AQ
    myfit[3][2]->SetParLimits(3, 0.0, .07);    //A1 jet amp
    myfit[3][2]->SetParLimits(4, 0.6, 1.2);     //phi width
    myfit[3][2]->SetParLimits(5, 0.1, .9);    //eta width
    
    myfit[4][0]->SetParLimits(0, -.15, .01);    //offset A0
    myfit[4][0]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[4][0]->SetParLimits(2, 0.0,  .1);    //v2     AQ
    myfit[4][0]->SetParLimits(3, 0.0, .2);    //A1 jet amp
    myfit[4][0]->SetParLimits(4, 0.4, .8);     //phi width
    myfit[4][0]->SetParLimits(5, 0.4, .8);    //eta width
	
	myfit[4][1]->SetParLimits(0, -.02, .01);    //offset A0
    myfit[4][1]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[4][1]->SetParLimits(2, 0.0,  .05);    //v2     AQ
    myfit[4][1]->SetParLimits(3, 0.0, .08);    //A1 jet amp
    myfit[4][1]->SetParLimits(4, 0.1, .7);     //phi width
    myfit[4][1]->SetParLimits(5, 0.1, .7);    //eta width
	
	myfit[4][2]->SetParLimits(0, -.04, .01);    //offset A0
    myfit[4][2]->SetParLimits(1, -.05,  0.0);    //v1     AD
    myfit[4][2]->SetParLimits(2, 0.0,  .02);    //v2     AQ
    myfit[4][2]->SetParLimits(3, 0.0, .04);    //A1 jet amp
    myfit[4][2]->SetParLimits(4, 0.1, .7);     //phi width
    myfit[4][2]->SetParLimits(5, 0.1, .8);    //eta width 
     */








    