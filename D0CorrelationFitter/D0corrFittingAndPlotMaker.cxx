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


using namespace std;


double calcScaleFactorforFit(double sigma)
{
    double preFactor = TMath::Sqrt(TMath::PiOver2());
    
    return preFactor*sigma*.5;
    
}
 
int getPtBinRange(int k); 

int D0corrFittingAndPlotMaker(){


    bool USE_SIX_PARAMETER_STANDARD     = false;                  //currently using for mid-central
    bool USE_TWO_2D_GAUSSIANS           = false;
    bool USE_TWO_2D_GAUSSIANS_PLUS_QUAD = false;
    bool USE_2D_SS_GAUSS_AS_GAUSS_WITH_DIPOLE_PLUS_QUAD = true;   //CURRENTLY USING FOR CENTRAL BIN
    bool OFFSET_QUAD_JETPEAK_AS_1D_GAUSS = false;                 //CURRENTLY USING FOR PERIPHERAL

	bool AS_ONLY_FIT_PERIPHERAL = false;
    bool AS_ONLY_FIT_MID_CENTRAL = false;
    bool AS_ONLY_FIT_CENTRAL = false;
    
    if((USE_SIX_PARAMETER_STANDARD+USE_TWO_2D_GAUSSIANS+USE_TWO_2D_GAUSSIANS_PLUS_QUAD+USE_2D_SS_GAUSS_AS_GAUSS_WITH_DIPOLE_PLUS_QUAD+OFFSET_QUAD_JETPEAK_AS_1D_GAUSS+
       AS_ONLY_FIT_MID_CENTRAL+AS_ONLY_FIT_CENTRAL) > 1) { cout << endl <<"Multiple flags true!!!!!!!" << endl; return 0;}
    
    double nBarPrime[3] = {74.4, 351, 1027};
    
    double numFitParameters = 0;
    
    int NUM_PT_BINS = 3;
    int FIRST_PT_BIN = 2;
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 13;
    double ETA_RANGE = 1.69; // 1.69;  //1.38
    double phiBinShift = (TMath::Pi()/12.0);
    
    double v2Array[5][3];
    double v2ArrayStatError[5][3];

    TString inputFileName  = "d0HadronCorrMaker_OUTPUT_12_phi_13_eta_new_dStar_range.root";
	TString outputFileName = "fittingOutput_OUTPUT_12_phi_13_eta_new_dStar_range.root";
    TString outputFolder   = "FitParameterPlots/";
    TString path = "C:/Users/ajentsch/desktop/D0_Analysis_Code_Local_Copies/D0CorrelationFitter/";
	
	TString fullySubtractedLabelCustomCentrality = "Full_DStar_Correction_FullSubtractedCorr_SideBand__PtBin_";
	TString ptBinLabel[5] = {"0", "1", "2", "3", "4"};
    TString centralityBin[3] = {"Peripheral", "MidCentral", "Central"};
    TString fileType = ".png";
    
    TString ptLabel     = "PtBin_";
    TString dipoleLabel = "Dipole_AD_";
    
    TString offsetLabel = "Offset_A0_";
    TString quadLabel   = "Quadrupole_AQ_";
    
    TString jetAmpLabel = "Same_Side_Amplitude_A_{SS}_";
    TString sigEtaLabel = "Same_Side_Eta_Width_";
    TString sigPhiLabel = "Same_Side_Phi_Width";
    
    TString awaySideAmpLabel = "AwaySideAmp_";
    TString awaySidePhiWidthLabel = "AwaySide_Phi_Width_";
    TString awaySideEtaWidthLabel = "AwaySide_Eta_Width_";
	
	TString phiProj = "Phi_projection_";
	TString str1;
   
	TH2D* fullySubtractedCorrCustomCentralitySideBand[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProj[5][3];
	TH2D* fullySubtractedCorrCustomCentralitySideBandFit[5][3];
    TH2D* fullySubtractedCorrCustomCentralitySideBandRes[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[5][3];
    
    TH2D* fullySubtractedCorrCustomCentralitySideBandPileupCorrected[5][3];
    
    TH1D* fitParameterOffsetPlot[3];
    TH1D* fitParameterDipolePlot[3];
    TH1D* fitParameterQuadrupolePlot[3];
    TH1D* fitParameterJetAmpPlot[3];
    TH1D* fitParameterSigmaEtaPlot[3];
    TH1D* fitParameterSigmePhiPlot[3];
    
    //TH1D* fitParameterOffsetPlotVsCent = new TH1D(offsetLabel, offsetLabel, 3, 0, 3);
    //TH1D* fitParameterDipolePlotVsCent = new TH1D(dipoleLabel, dipoleLabel, 3, 0, 3);
    //TH1D* fitParameterQuadrupolePlotVsCent = new TH1D(quadLabel, quadLabel, 3, 0, 3);
    //TH1D* fitParameterJetAmpPlotVsCent = new TH1D(jetAmpLabel, jetAmpLabel, 3, 0, 3);
    //TH1D* fitParameterSigmaEtaPlotVsCent = new TH1D(sigEtaLabel, sigEtaLabel, 3, 0, 3);
    //TH1D* fitParameterSigmePhiPlotVsCent = new TH1D(sigPhiLabel, sigPhiLabel, 3, 0, 3);
    
    TH1D* fitParameterOffsetPlotVsCent[5];
    TH1D* fitParameterDipolePlotVsCent[5];
    TH1D* fitParameterQuadrupolePlotVsCent[5];
    TH1D* fitParameterJetAmpPlotVsCent[5];
    TH1D* fitParameterSigmaEtaPlotVsCent[5];
    TH1D* fitParameterSigmePhiPlotVsCent[5];
    TH1D* fitParameterASAmpPlotVsCent[5];
    TH1D* fitParameterASSigmaEtaPlotVsCent[5];
    TH1D* fitParameterASSigmaPhiPlotVsCent[5];
    
    
    TH1D * jetVolumePlot[5];
   // TH2D* errorOnBin[3];
    TH2D* sigmaPerBin[5][3];
    
    //TString offsetLabel = "Offset_A0_";
    //TString quadLabel   = "Quadrupole_AQ_";
    
    //TString jetAmpLabel = "JetAmp_A1_";
    //TString sigEtaLabel = "SigmaEta_";
    //TString sigPhiLabel = "SigmaPhi_";
    
    //TString awaySideAmpLabel = "AwaySideAmp_";
    //TString awaySidePhiWidthLabel = "AwaySide_Phi_Width_";
    //TString awaySideEtaWidthLabel = "AwaySide_Eta_Width_";
    
    for(int k = 2; k < 5; k++){
    
        str1 = offsetLabel + ptLabel + ptBinLabel[k];
        fitParameterOffsetPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = dipoleLabel + ptLabel + ptBinLabel[k];
        fitParameterDipolePlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = quadLabel + ptLabel + ptBinLabel[k];
        fitParameterQuadrupolePlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = jetAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterJetAmpPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = sigEtaLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaEtaPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = sigPhiLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmePhiPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = awaySideAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterASAmpPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = awaySideEtaWidthLabel + ptLabel + ptBinLabel[k];
        fitParameterASSigmaEtaPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = awaySidePhiWidthLabel + ptLabel + ptBinLabel[k];
        fitParameterASSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterQuadrupolePlotVsCent[k]->SetLabelSize(.07);
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterJetAmpPlotVsCent[k]->SetLabelSize(.07);
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterSigmaEtaPlotVsCent[k]->SetLabelSize(.07);
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterSigmePhiPlotVsCent[k]->SetLabelSize(.07);
        fitParameterSigmePhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterSigmePhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterSigmePhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterASSigmaEtaPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterASSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterASSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASSigmaEtaPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetLabelSize(.07);
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
    }
    
    float Lower[4];
    Lower[0] = 0;
    Lower[1] = 1;
    Lower[2] = 2;
    //Lower[3] = 3;
    //Lower[4] = 5;
    Lower[3] = 10;
    
    for(int i = 0; i < 3; i++){
    
        str1 = offsetLabel + centralityBin[i];
        fitParameterOffsetPlot[i]     = new TH1D(str1, str1, 3, Lower);
        str1 = dipoleLabel + centralityBin[i];
        fitParameterDipolePlot[i]     = new TH1D(str1, str1, 3, Lower);
        str1 = quadLabel + centralityBin[i];
        fitParameterQuadrupolePlot[i] = new TH1D(str1, str1, 3, Lower);
        str1 = jetAmpLabel + centralityBin[i];
        fitParameterJetAmpPlot[i]     = new TH1D(str1, str1, 3, Lower);
        str1 = sigEtaLabel + centralityBin[i];
        fitParameterSigmaEtaPlot[i]   = new TH1D(str1, str1, 3, Lower);
        str1 = sigPhiLabel + centralityBin[i];
        fitParameterSigmePhiPlot[i]   = new TH1D(str1, str1, 3, Lower);
     
    }

    TFile * inputFile  = new TFile(inputFileName);
    TFile * outputFile = new TFile(outputFileName, "RECREATE");
 
 
    TCanvas *c = new TCanvas("c3", "Histograms", 1100, 850);
 
    TString fitFunctionLabel = "parameterFitFunction";
 
    //TF2 *myfit   = new TF2("parameterFitFunction","[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
    
    TF2 *myfit[5][3];   
    
    TF1 *myfit1D = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x)" , -1.57, 4.71);
    
    TF1 *fitProjection = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x) + [3]*[4]*(exp(-0.5*(x*x)/([5]*[5]))+exp(-0.5*((x-6.28)*(x-6.28))/([5]*[5])))" , -1.57, 4.71);
    
    TF1 *quadrupole1    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole2    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole3    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    
    TString fitLabel = "_Fit";
    TString resLabel = "_Residual";
    TString pileupLabel = "PileupCorrected";
    
    ofstream fitParameterFile2D;
    TString fitParameterFileName2D = "Extracted_Parameters_from_fitting_2D.txt";
    TString fitParamterOutputFile2D = path + fitParameterFileName2D;
    fitParameterFile2D.open (fitParamterOutputFile2D);
    
    ofstream fitParameterFile1D;
    TString fitParameterFileName1D = "Extracted_Parameters_from_fitting_1D.txt";
    TString fitParamterOutputFile1D = path + fitParameterFileName1D;
    fitParameterFile1D.open (fitParamterOutputFile1D);
    
    double v2_2D;
    double v2_2D_Error;
    double v2_1D;
    double v2_1D_Error;
    
    int bin1;
    int bin2;
    int bin3;
    int bin4;
    
    double GaussIntegralParameter;
    
    int ptBinRange;
    
    for(int k = 2; k < 5; k++){
        for(int i = 0; i < 3; i++){ 
 
            str1 = fitFunctionLabel + ptBinLabel[k] + centralityBin[i];
            if(USE_SIX_PARAMETER_STANDARD) { myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
            if(USE_TWO_2D_GAUSSIANS) { myfit[k][i] = new TF2(str1,"[0]*(exp(-0.5*(y*y)/([1]*[1]))+exp(-0.5*((y-6.28)*(y-6.28))/([1]*[1])))*exp(-0.5*(x*x)/([2]*[2])) + [3]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([4]*[4]))+exp(-0.5*((y-9.42478)*(y-9.42478))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
            if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD) { myfit[k][i] = new TF2(str1,"[0] + 2*[1]*cos(2*y) + [2]*(exp(-0.5*(y*y)/([3]*[3])) + exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4])) + [5]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([6]*[6]))+exp(-0.5*((y+3.14159)*(y+3.14159))/([6]*[6])))*exp(-0.5*(x*x)/([7]*[7]))",-ETA_RANGE, ETA_RANGE, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
            if(USE_2D_SS_GAUSS_AS_GAUSS_WITH_DIPOLE_PLUS_QUAD) { myfit[k][i] = new TF2(str1,"[0] + 2*[1]*cos(2*y) + [2]*(exp(-0.5*(y*y)/([3]*[3])) + exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4])) + [5]*cos(y)*exp(-0.5*(x*x)/([6]*[6]))",-ETA_RANGE, ETA_RANGE, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
            if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { myfit[k][i] = new TF2(str1,"[0] + 2*[1]*cos(2*y) + [2]*(exp(-0.5*(y*y)/([3]*[3])) + exp(-0.5*((y-6.283185)*(y-6.283185))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4])) + [5]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([6]*[6]))+exp(-0.5*((y+3.14159)*(y+3.14159))/([6]*[6])))", -ETA_RANGE, ETA_RANGE, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift); }

			if(AS_ONLY_FIT_PERIPHERAL) { myfit[k][i] = new TF2(str1,"[0] + 2*[1]*cos(2*y) + [2]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([3]*[3]))+exp(-0.5*((y+3.14159)*(y+3.14159))/([3]*[3])))", -ETA_RANGE, ETA_RANGE, TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift); }
            if(AS_ONLY_FIT_MID_CENTRAL) { myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) +2*[2]*cos(2*y)", -ETA_RANGE, ETA_RANGE, TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift); }
            //if(AS_ONLY_FIT_CENTRAL) { myfit[k][i] = new TF2(str1,"[0] + 2*[1]*cos(2*y) +  [2]*cos(y)*exp(-0.5*(x*x)/([3]*[3]))", -ETA_RANGE, ETA_RANGE, TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift); }
			if(AS_ONLY_FIT_CENTRAL) { myfit[k][i] = new TF2(str1,"[0] + 2*[1]*cos(2*y) +  [2]*exp(-0.5*(x*x)/([3]*[3]))*cos(y)", -ETA_RANGE, ETA_RANGE, TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift); }

		}
	}
    
   
    if(USE_SIX_PARAMETER_STANDARD) { //use for mid-central
    
    
    
    
        myfit[2][0]->SetParLimits(0, -.01, 0.0);    //offset A0
        myfit[2][0]->SetParLimits(1, -.01,  0.0);    //v1     AD
        myfit[2][0]->SetParLimits(2, 0.0,  .009);    //v2     AQ
        //myfit[2][0]->FixParameter(2, .0084); //v2 AQ
        myfit[2][0]->SetParLimits(3, 0.06, 0.13);    //A1 jet amp
        myfit[2][0]->SetParLimits(4, 0.2, .5);     //phi width
        myfit[2][0]->SetParLimits(5, 0.2, .6);    //eta width
        
         //                         A0      AD     AQ       A1    sPhi    sEta  
       // myfit[2][1]->SetParameters(-.012, -.0208,   .00478,  .0325,  .67,    4.0);
        
        myfit[2][1]->SetParLimits(0, -.02, 0.0);    //offset A0
        //myfit[2][1]->SetParLimits(1, -.03, -0.003);    //v1     AD
        //myfit[2][1]->SetParLimits(2, 0.003,  .01);    //v2     AQ
        //myfit[2][1]->FixParameter(1, -.0177424); //AD     --smaller eta range
        //myfit[2][1]->FixParameter(2, .00605651); //v2 AQ  --smaller eta range
		
		myfit[2][1]->FixParameter(1, -.0165570); //AD
        myfit[2][1]->FixParameter(2, .00647402); //v2 AQ
        myfit[2][1]->SetParLimits(3, 0.02, .08);    //A1 jet amp
        myfit[2][1]->SetParLimits(4, 0.2, .9);     //phi width
        myfit[2][1]->SetParLimits(5, .5, 4.3);    //eta width
        
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][2]->SetParameters(-.0007, -.002, .003, .02,  .32, .37);
        
        myfit[2][2]->SetParLimits(0, -.01, 0.0);    //offset A0
        myfit[2][2]->SetParLimits(1, -.02,  0.0);    //v1     AD
        myfit[2][2]->SetParLimits(2, 0.00,  .007);    //v2     AQ
        //myfit[2][2]->FixParameter(2, .0036); //v2 AQ
        myfit[2][2]->SetParLimits(3, 0.0, .05);    //A1 jet amp
        myfit[2][2]->SetParLimits(4, .5, 1.0);     //phi width
        myfit[2][2]->SetParLimits(5, .6, 1.9);    //eta width 
    
    
       
    }
    
    if(USE_TWO_2D_GAUSSIANS) {
    //                               A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][0]->SetParameters(-.003,  -.02,   .02,   .08,  .3, .4);
    
        myfit[2][0]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[2][0]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[2][0]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[2][0]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[2][0]->SetParLimits(4, 1.5, 2.0);     //away-side phi width
        myfit[2][0]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        
         //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][1]->SetParameters(-.005, -.01,   .01,  .03,  .7, 1.4);
        
        myfit[2][1]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[2][1]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[2][1]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[2][1]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[2][1]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[2][1]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][2]->SetParameters(-.0007, -.002, .003, .02,  .32, .37);
        
        myfit[2][2]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[2][2]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[2][2]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[2][2]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[2][2]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[2][2]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
    
    
       
 
    }
    
    if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD) {
    //                               A0      AQ    A1    JetPhi  JetEta    ASAmp   ASPhi  ASEta  
        //myfit[2][0]->SetParameters(-.003,  .01,   .8,    .6,      .5,     .02,     .9,    1.2);
    
        myfit[2][0]->SetParLimits(0, -.02, 0.0);    //offset  
        myfit[2][0]->SetParLimits(1, 0.003, .01);    //quadrupole  
        myfit[2][0]->SetParLimits(2, 0.07,  .14);   //jet amp
        myfit[2][0]->SetParLimits(3, 0.1, .4);    //jet phi width
        myfit[2][0]->SetParLimits(4, 0.1, .4);    //jet eta width 
        myfit[2][0]->SetParLimits(5, 0.01, .03);    //AS amp
        myfit[2][0]->SetParLimits(6, 0.4, .8);    //AS phi width
        myfit[2][0]->SetParLimits(7, .3, 4.0);    //AS eta width
        
        
        //                           A0      AQ    A1    JetPhi  JetEta    ASAmp   ASPhi  ASEta  
        //myfit[2][0]->SetParameters(-.003,  .01,   .8,    .6,      .5,     .02,     .9,    .9);
        
        myfit[2][1]->SetParLimits(0, -.02, 0.0);    //offset  
        myfit[2][1]->SetParLimits(1, 0.0, .03);    //quadrupole  
        myfit[2][1]->SetParLimits(2, 0.0,  .03);   //jet amp
        myfit[2][1]->SetParLimits(3, 0.0, 1.0);    //jet phi width
        myfit[2][1]->SetParLimits(4, 0.0, 2.0);    //jet eta width 
        myfit[2][1]->SetParLimits(5, 0.0, .05);    //AS amp
        myfit[2][1]->SetParLimits(6, 0.0, 1.2);    //AS phi width
        myfit[2][1]->SetParLimits(7, 0.0, 3.0);    //AS eta width
        
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][2]->SetParameters(-.0007, -.002, .003, .02,  .32, .37);
        
        //myfit[2][2]->SetParLimits(0, -.004, -0.003);    //offset  
        //myfit[2][2]->SetParLimits(1, 0.002, .0023);    //quadrupole  
        //myfit[2][2]->SetParLimits(2, 0.0095,  .01);   //jet amp
        //myfit[2][2]->SetParLimits(3, 0.48, .52);    //jet phi width
        //myfit[2][2]->SetParLimits(4, 0.68, .72);    //jet eta width 
        //myfit[2][2]->SetParLimits(5, 0.0093, .011);    //AS amp
        //myfit[2][2]->SetParLimits(6, .72, .8);    //AS phi width
        //myfit[2][2]->SetParLimits(7, 0.19, .23);    //AS eta width
    
        myfit[2][2]->SetParLimits(0, -.01, 0.0);    //offset 
        myfit[2][2]->SetParLimits(4, 0.5, 1.2);    //jet eta width         
        
 
    }
 
    if(USE_2D_SS_GAUSS_AS_GAUSS_WITH_DIPOLE_PLUS_QUAD) { //CENTRAL FIT
    
    //                               A0      AQ    A1    JetPhi  JetEta    ASAmp   ASPhi  ASEta  
        //myfit[2][0]->SetParameters(-.003,  .01,   .8,    .6,      .5,     .02,     .9,    1.2);
    
        myfit[2][0]->SetParLimits(0, -.02, 0.0);    //offset  
        myfit[2][0]->SetParLimits(1, 0.0, .015);    //quadrupole  
        myfit[2][0]->SetParLimits(2, 0.0,  .1);   //jet amp
        myfit[2][0]->SetParLimits(3, 0.48, .52);    //jet phi width
        myfit[2][0]->SetParLimits(4, 0.68, .72);    //jet eta width 
        myfit[2][0]->SetParLimits(5, 0.0, .03);    //AS dipole
        myfit[2][0]->SetParLimits(6, 1.0, 2.2);    //AS eta width
        
        
        //                           A0      AQ    A1    JetPhi  JetEta    ASAmp   ASPhi  ASEta  
        myfit[2][1]->SetParameters(-.003,  .01,   .8,    .6,      .5,     .02,     .9,    .9);
        
        myfit[2][1]->SetParLimits(0, -.02, 0.0);    //offset  
        myfit[2][1]->SetParLimits(1, 0.0, .03);    //quadrupole  
        myfit[2][1]->SetParLimits(2, 0.0,  .03);   //jet amp
        myfit[2][1]->SetParLimits(3, 0.0, 1.0);    //jet phi width
        myfit[2][1]->SetParLimits(4, 0.0, 2.0);    //jet eta width 
        myfit[2][1]->SetParLimits(5, 0.0, .05);    //AS dipole
        myfit[2][1]->SetParLimits(6, 0.0, 2.0);    //AS eta width
        
        //                         A0       AQ     A1    sPhi  sEta  AS dipole   AS eta width
        myfit[2][2]->SetParameters(-.001, 0.0,  .04,     .75,    1.32,     -.0145       ,   1.06      );
        
        myfit[2][2]->SetParLimits(0, -.015, 0.0);    //offset  
        myfit[2][2]->SetParLimits(1, 0.0, .003);    //quadrupole  
        myfit[2][2]->SetParLimits(2, 0.008,  .19);   //jet amp
        myfit[2][2]->SetParLimits(3, 0.3, .9);    //jet phi width
        myfit[2][2]->SetParLimits(4, .5, 1.5);    //jet eta width 
        myfit[2][2]->SetParLimits(5, -0.03, 0.0);    //AS dipole
        myfit[2][2]->SetParLimits(6, 0.05, 3.5);    //AS eta width
    
    
    }
    
    if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS){ //use for peripheral
    
	    //AQ = 0.005  (assuming the form 2*A_Q*cos(2*phi_Delta) )
		//AS 1DG amp = 0.028, sigma_phi = 0.52  and flat on eta
		//SS 2DG amp = 0.11, sigma_eta = 0.27, sigma_phi = 0.3 
		//This should be at the true, chi-square minimum.
	
        myfit[2][0]->SetParameters(-0.011, .004, .11, .3, .3, .028, .52);
        
        myfit[2][0]->SetParLimits(0, -.015, -0.002);    //offset A0
        myfit[2][0]->SetParLimits(1, 0.002,  0.005);    //v2     AQ
        myfit[2][0]->SetParLimits(2, 0.01, .14);    //A1 jet amp
        myfit[2][0]->SetParLimits(3, 0.2, .4);     //phi width
        myfit[2][0]->SetParLimits(4, 0.24, .4);    //eta width 
        myfit[2][0]->SetParLimits(5, .01,  0.032);    //AS amp
        myfit[2][0]->SetParLimits(6, 0.48, .7);     //AS Phi Width
        
        
        myfit[2][1]->SetParameters(-0.011, .005, .11, .3, .27, .028, .52);
       
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
   
   
    if(AS_ONLY_FIT_CENTRAL) {
    
		myfit[2][2]->SetParameters(-.00762196, 0.001, -.0145896, 1.06427);
	
        myfit[2][2]->SetParLimits(0, -.02, 0.0);    //offset A0
        myfit[2][2]->SetParLimits(1, 0.0,  0.002);    //v2     AQ
        myfit[2][2]->SetParLimits(2, -0.015, -0.013);    //AS dipole amp
        myfit[2][2]->SetParLimits(3, 0.9, 1.1);     //AS eta width
		
    
    }
	
	if(AS_ONLY_FIT_PERIPHERAL) {
	
		myfit[2][2]->SetParLimits(0, -.02, 0.0);    //offset A0
        myfit[2][2]->SetParLimits(1, 0.001,  0.01);    //v2     AQ
        myfit[2][2]->SetParLimits(2, 0.0, 0.04);    //AS amp
        myfit[2][2]->SetParLimits(3, 0.2, 1.2);     //AS phi
		
	
	
	}
 
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){ 
 
        
            str1 = fullySubtractedLabelCustomCentrality  + ptBinLabel[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBand[k][i] = (TH2D*)inputFile->Get(str1);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetDirectory(0);
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality  + ptBinLabel[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
            //SET THE ETA RANGE HERE
            
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY("",4,10,"");
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetNameTitle(str1, str1);
            //
            
            str1 = fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 = fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 = phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
 
            str1 = fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
 
            fitParameterFile2D << "_________________________________________________________________________________________________________" << endl;
            fitParameterFile1D << "_________________________________________________________________________________________________________" << endl;
 
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone();

            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            
             if(AS_ONLY_FIT_PERIPHERAL){
                fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            }
			
            if(AS_ONLY_FIT_MID_CENTRAL){
                fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            }
            
            if(AS_ONLY_FIT_CENTRAL) {
            
                fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            }
            
            myfit1D->SetParLimits(0, -.5, .01); //offset A0
            myfit1D->SetParLimits(1, -.1, 0.0); //v1     AD
            myfit1D->SetParLimits(2, 0.0, .1); //v2     AQ
            
            //////////////////////2d fitting///////////////////////////////////////
            
            for(int eta = 1; eta < 14; eta++){
                for(int phi = 1; phi < 13; phi++){
            
                    fullySubtractedCorrCustomCentralitySideBand[k][i]->SetBinError(eta, phi, 2.0*fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinError(eta, phi));
                    
                }
            }
            
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Fit(myfit[k][i], "R0");
            
            
            
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone("hfit");
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Eval(myfit[k][i]);
            
            GaussIntegralParameter = calcScaleFactorforFit(myfit[k][i]->GetParameter(5));
            
            cout << "Pt bin: " << k << "   Cent Bin: " << i << "   ProjectionNumber: " << GaussIntegralParameter << endl;
            
            fitParameterFile2D << endl << "pt bin (3 is integrated): " << k << "    Centrality Bin (0 - per, 1 - mid cent, 2 - cent): " << i << endl;
            fitParameterFile2D << endl << "A0" << "\t\t" << "AD" << "\t\t" <<"AQ" << "\t\t" <<"Jet Amp A1" << "\t\t" << "SigEta" << "\t\t" <<"SigPhi" << "\t\t" << endl;
            
            fitParameterFile2D << endl << myfit[k][i]->GetParameter(0) << "\t" << myfit[k][i]->GetParameter(1) << "\t" 
                                       << myfit[k][i]->GetParameter(2) << "\t" << myfit[k][i]->GetParameter(3) << "\t" 
                                       << myfit[k][i]->GetParameter(5) << "\t" << myfit[k][i]->GetParameter(4) << endl; 
            
            fitParameterFile2D << myfit[k][i]->GetParError(0) << "\t" << myfit[k][i]->GetParError(1) << "\t" 
                               << myfit[k][i]->GetParError(2) << "\t" << myfit[k][i]->GetParError(3) << "\t" 
                               << myfit[k][i]->GetParError(5) << "\t" << myfit[k][i]->GetParError(4) << endl;            
            
            //fullySubtractedCorrCustomCentralitySideBandRes[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Add(fullySubtractedCorrCustomCentralitySideBand[k][i], fullySubtractedCorrCustomCentralitySideBandFit[k][i], 1, -1);
            
            str1 = path + outputFolder + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fileType;
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
            c->SaveAs(str1);
 
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Write();
            
            str1 = path + outputFolder + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Draw("SURF1");
            c->SaveAs(str1);
 
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Write();
            
            str1 = path + outputFolder + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Draw("SURF1");
            c->SaveAs(str1); 
            
            /////////////////////////////1d fitting///////////////////////////////////////////
    
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Fit(myfit1D, "qR0");
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hfit");
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Eval(myfit1D);
 
            fitParameterFile1D << endl << "pt bin (3 is integrated): " << k << "    Centrality Bin (0 - per, 1 - mid cent, 2 - cent): " << i << endl;
            fitParameterFile1D << endl << "A0" << "\t\t" << "AD" << "\t\t" <<"AQ" << endl; // "\t\t" << "Jet Amp A1" << "\t\t" << "SigPhi" << endl;
            
            fitParameterFile1D << endl << myfit1D->GetParameter(0) << "\t" << myfit1D->GetParameter(1) << "\t" << myfit1D->GetParameter(2) << endl; //<< "\t" << myfit1D->GetParameter(3) << "\t" << myfit1D->GetParameter(5) << endl; 
            fitParameterFile1D << myfit1D->GetParError(0) << "\t" << myfit1D->GetParError(1) << "\t" << myfit1D->GetParError(2) << endl;//<< "\t" << myfit1D->GetParError(3) << "\t" << myfit1D->GetParError(5) << endl; 
            
            v2_2D = TMath::Sqrt(myfit[k][i]->GetParameter(2));
            v2_2D_Error = v2_2D*(.5)*(myfit[k][i]->GetParError(2)/myfit[k][i]->GetParameter(2));
            
            v2Array[k][i] = v2_2D;
            v2ArrayStatError[k][i] = v2_2D_Error;
            
            
            v2_1D = TMath::Sqrt(myfit1D->GetParameter(2));
            v2_1D_Error = v2_2D*(.5)*(myfit1D->GetParError(2)/myfit1D->GetParameter(2));
            
            //fitParameterFile2D << "V2 = " << v2_2D << " +/- " << v2_2D_Error << endl;
            //fitParameterFile1D << "V2 = " << v2_1D << " +/- " << v2_1D_Error << endl;
            fitParameterFile2D << endl;
            fitParameterFile1D << endl;
            
            fitProjection->SetParameter(0, myfit[k][i]->GetParameter(0));
            fitProjection->SetParameter(1, myfit[k][i]->GetParameter(1));
            fitProjection->SetParameter(2, myfit[k][i]->GetParameter(2));
            fitProjection->SetParameter(3, myfit[k][i]->GetParameter(3));
            fitProjection->SetParameter(5, myfit[k][i]->GetParameter(4));
            fitProjection->SetParameter(4, GaussIntegralParameter);
            
            quadrupole1->SetParameter(0, myfit[k][i]->GetParameter(2));
            quadrupole2->SetParameter(0, myfit1D->GetParameter(2));
            quadrupole3->SetParameter(0, fitProjection->GetParameter(2));
            
            fitProjection->SetLineColor(2); //red
            quadrupole1->SetLineColor(3); // Green, from 2D fit
            quadrupole2->SetLineColor(4); //blue, from 1D fit
            myfit1D->SetLineColor(7); // light blue
            quadrupole3->SetLineColor(41);
            
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Add(myfit1D, -1);
            
            
            str1 = path + outputFolder + phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel + fileType;
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Draw();
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw("EX0");
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetMarkerStyle(20);
            
            
            
            fitProjection->Draw("SAME");
            quadrupole1->Draw("SAME");
            quadrupole2->Draw("SAME");
            //quadrupole3->Draw("SAME");
            myfit1D->Draw("SAME");
            c->SaveAs(str1);
 
            str1 = path + outputFolder + phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Draw();
            c->SaveAs(str1); 
 
            ptBinRange = k+1;
 
            if(k == 0) {continue;}
 
            fitParameterOffsetPlot[i]->Fill(ptBinRange, myfit[k][i]->GetParameter(0));
            fitParameterDipolePlot[i]->Fill(ptBinRange, myfit[k][i]->GetParameter(1));
            fitParameterQuadrupolePlot[i]->Fill(ptBinRange, myfit[k][i]->GetParameter(2));
            fitParameterJetAmpPlot[i]->Fill(ptBinRange, myfit[k][i]->GetParameter(3));
            fitParameterSigmaEtaPlot[i]->Fill(ptBinRange, myfit[k][i]->GetParameter(5));
            fitParameterSigmePhiPlot[i]->Fill(ptBinRange, myfit[k][i]->GetParameter(4));
            
            fitParameterOffsetPlot[i]->SetBinError(ptBinRange, myfit[k][i]->GetParError(0));
            fitParameterDipolePlot[i]->SetBinError(ptBinRange, myfit[k][i]->GetParError(1));
            fitParameterQuadrupolePlot[i]->SetBinError(ptBinRange, myfit[k][i]->GetParError(2));
            fitParameterJetAmpPlot[i]->SetBinError(ptBinRange, myfit[k][i]->GetParError(3));
            fitParameterSigmaEtaPlot[i]->SetBinError(ptBinRange, myfit[k][i]->GetParError(5));
            fitParameterSigmePhiPlot[i]->SetBinError(ptBinRange, myfit[k][i]->GetParError(4));
 
           
 
 
        }
    }    
    
    for(int i = 0; i < 3; i++){
    
        fitParameterOffsetPlot[i]->SetMarkerColor(2);
        fitParameterDipolePlot[i]->SetMarkerColor(2);
        fitParameterQuadrupolePlot[i]->SetMarkerColor(2);
        fitParameterJetAmpPlot[i]->SetMarkerColor(2);
        fitParameterSigmaEtaPlot[i]->SetMarkerColor(2);
        fitParameterSigmePhiPlot[i]->SetMarkerColor(2);
        
        fitParameterOffsetPlot[i]->SetMarkerStyle(8);
        fitParameterDipolePlot[i]->SetMarkerStyle(8);
        fitParameterQuadrupolePlot[i]->SetMarkerStyle(8);
        fitParameterJetAmpPlot[i]->SetMarkerStyle(8);
        fitParameterSigmaEtaPlot[i]->SetMarkerStyle(8);
        fitParameterSigmePhiPlot[i]->SetMarkerStyle(8);
    
        fitParameterOffsetPlot[i]->Write();
        fitParameterDipolePlot[i]->Write();
        fitParameterQuadrupolePlot[i]->Write();
        fitParameterJetAmpPlot[i]->Write();
        fitParameterSigmaEtaPlot[i]->Write();
        fitParameterSigmePhiPlot[i]->Write();
    }
    
    int bin = 0;
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){
    
            bin = i + 1;
        
            if(USE_SIX_PARAMETER_STANDARD){
            
                fitParameterOffsetPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(0));
                fitParameterDipolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(1));
                fitParameterQuadrupolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(2));
                fitParameterJetAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(3));
                fitParameterSigmaEtaPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(5));
                fitParameterSigmePhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(4));
                    
                fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterDipolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
                fitParameterQuadrupolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(2));
                fitParameterJetAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(3));
                fitParameterSigmaEtaPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(5));
                fitParameterSigmePhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(4));
            }
            
            if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD){
            
                fitParameterOffsetPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(0));      //offset
                fitParameterQuadrupolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(1));  //quad
                fitParameterJetAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(2));      //jet amp
                fitParameterSigmePhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(3));    // jet phi
                fitParameterSigmaEtaPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(4));    //jet eta
                fitParameterASAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(5));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(6));  //AS phi
                fitParameterASSigmaEtaPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(7));  // AS eta
        
                    
                fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterQuadrupolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
                fitParameterJetAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(2));
                fitParameterSigmePhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(3));
                fitParameterSigmaEtaPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(4));
                fitParameterASAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParameter(5));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParameter(6));  //AS phi
                fitParameterASSigmaEtaPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParameter(7));  // AS eta
            
            }
            
            if(USE_2D_SS_GAUSS_AS_GAUSS_WITH_DIPOLE_PLUS_QUAD) { 
            
                fitParameterOffsetPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(0));      //offset
                fitParameterQuadrupolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(1));  //quad
                fitParameterJetAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(2));      //jet amp
                fitParameterSigmePhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(3));    // jet phi
                fitParameterSigmaEtaPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(4));    //jet eta
                fitParameterDipolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(5));      //dipole
                fitParameterASSigmaEtaPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(6));  // AS eta
        
                    
                fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterQuadrupolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
                fitParameterJetAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(2));
                fitParameterSigmePhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(3));
                fitParameterSigmaEtaPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(4));
                fitParameterDipolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(5));
                fitParameterASSigmaEtaPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParameter(6));  // AS eta
                
            }
            
            if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { 
            
                fitParameterOffsetPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(0));      //offset
                fitParameterQuadrupolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(1));  //quad
                fitParameterJetAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(2));      //jet amp
                fitParameterSigmePhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(3));    // jet phi
                fitParameterSigmaEtaPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(4));    //jet eta
                fitParameterASAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(5));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(6));  //AS phi
                    
                
                fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterQuadrupolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
                fitParameterJetAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(2));
                fitParameterSigmePhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(3));
                fitParameterSigmaEtaPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(4));
                fitParameterASAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(5));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(6));  //AS phi
                    
                
            }
            
        }
    }        
 
    fitParameterFile2D.close();
    fitParameterFile1D.close();
 
    TCanvas * finalCorrGridCanvas = new TCanvas("", "", 3600, 3400);
    finalCorrGridCanvas->Divide(3,3);
    
    TString ptBinTitle[5] = {"0 GeV/c < p_{t, D^{0}} < 1 GeV/c", "1 GeV/c < p_{t, D^{0}} < 2 GeV/c", "2 GeV/c < p_{t, D^{0}} < 3 GeV/c",   //change made from 2-3 to 2-10
                             "3 GeV/c < p_{t, D^{0}} < 4 GeV/c", "4 GeV/c < p_{t, D^{0}} < 10 GeV/c"}; 
    int padCount = 1;
    for(int i = 0; i < 3; i++){
        for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        
        
            finalCorrGridCanvas->cd(padCount);
            
            //if(i == 0) { fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle(ptBinTitle[k]); }
            //else fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle("");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle("");
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitle("#Delta#eta");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleOffset(1.5);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleSize(.065);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitle("#Delta#phi");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleOffset(1.5);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleSize(.065);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
            
            padCount++;
   
        }
    }
   
    finalCorrGridCanvas->SaveAs("nominal_12_phi_bins_13_eta_3_pt_bins.png");
    
    TCanvas * finalFitCorrGridCanvas = new TCanvas("", "", 4800, 1000);
    finalFitCorrGridCanvas->Divide(4,1);
    
    TH2D* errorOnBin[3];
    
    TString outputFitString = "Fitting_nominal_cuts_PtBin_";
    
    double sigmaValueSquare, sigmaValue;
    double data, model, error;
    double delEtaBinCenter, delPhiBinCenter;
    double maxZ, minZ;
    TString chiPlotLabel = "nSigma Residuals ";
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        
        for(int i = 0; i < 3; i++){
        
            padCount = 1;
           //cout << "CHECK IT!" << endl;
        
            finalFitCorrGridCanvas->cd(padCount);
            
            //if(i == 0) { fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle(ptBinTitle[k]); }
            //else fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle("");
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitle("#Delta#eta");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitle("#Delta#phi");
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitle("#Delta#rho/#rho_{mix}");
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleSize(.09);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleSize(.09); 
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitleSize(.07);
    
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetTitleOffset(1.15);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleOffset(1.15); 
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetTitleOffset(.9);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetLabelSize(.05);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetLabelSize(.05);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetLabelSize(.04);
			fullySubtractedCorrCustomCentralitySideBand[k][i]->GetZaxis()->SetNdivisions(5,3,0, 1);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
            
            //maxZ = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMaximum();
            //minZ = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMinimum();
			fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetTitle("#Delta#eta");
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->SetTitle("#Delta#phi");
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetZaxis()->SetTitle("#Delta#rho/#rho_{mix}");
            
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetTitleSize(.09);
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->SetTitleSize(.09); 
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetZaxis()->SetTitleSize(.07);
    
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetTitleOffset(1.15);
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->SetTitleOffset(1.15); 
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetZaxis()->SetTitleOffset(.9);
			fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetLabelSize(.05);
			fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetYaxis()->SetLabelSize(.05);
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetZaxis()->SetLabelSize(.04);
			fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetZaxis()->SetNdivisions(5,3,0, 1);
            
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetMaximum(fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMaximum());
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetMinimum(fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMinimum());
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetStats(0);
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetLineColor(2);
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetLineWidth(2);
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Draw("SURF SAME");
            
           
            
            padCount++;
            finalFitCorrGridCanvas->cd(padCount);
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetMaximum(fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMaximum());
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetMinimum(fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMinimum());
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetStats(0);
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetLineColor(1);
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetLineWidth(1);
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Draw("SURF1");
            padCount++;
            finalFitCorrGridCanvas->cd(padCount);
			fullySubtractedCorrCustomCentralitySideBandRes[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetTitle("#Delta#eta");
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->SetTitle("#Delta#phi");
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->CenterTitle();
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetZaxis()->SetTitle("#Delta#rho/#rho_{mix}");
            
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetTitleSize(.09);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->SetTitleSize(.09); 
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetZaxis()->SetTitleSize(.07);
    
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetTitleOffset(1.15);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->SetTitleOffset(1.15); 
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetZaxis()->SetTitleOffset(.9);
			fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetLabelSize(.05);
			fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetYaxis()->SetLabelSize(.05);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetZaxis()->SetLabelSize(.04);
			fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetZaxis()->SetNdivisions(5,3,0, 1);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->SetMaximum(fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMaximum());
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->SetMinimum(fullySubtractedCorrCustomCentralitySideBand[k][i]->GetMinimum());
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Draw("SURF1");
            padCount++;
            finalFitCorrGridCanvas->cd(padCount);
            
            
            str1 = chiPlotLabel + centralityBin[i];
            sigmaPerBin[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
            if(AS_ONLY_FIT_MID_CENTRAL){
               sigmaPerBin[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
               
            }
            
            if(AS_ONLY_FIT_CENTRAL) {
            
                sigmaPerBin[k][i]->GetYaxis()->SetRangeUser(TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
                
            }
			
			fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetStats(0);
            sigmaPerBin[k][i]->GetXaxis()->SetTitle("#Delta#eta");
            sigmaPerBin[k][i]->GetXaxis()->CenterTitle();
            sigmaPerBin[k][i]->GetYaxis()->SetTitle("#Delta#phi");
            sigmaPerBin[k][i]->GetYaxis()->CenterTitle();
            //sigmaPerBin[k][i]->GetZaxis()->SetTitle("nSigma");
            
            sigmaPerBin[k][i]->GetXaxis()->SetTitleSize(.09);
            sigmaPerBin[k][i]->GetYaxis()->SetTitleSize(.09); 
            sigmaPerBin[k][i]->GetZaxis()->SetTitleSize(.04);
    
            sigmaPerBin[k][i]->GetXaxis()->SetTitleOffset(.55);
            sigmaPerBin[k][i]->GetYaxis()->SetTitleOffset(.55); 
            sigmaPerBin[k][i]->GetZaxis()->SetTitleOffset(.65);
			sigmaPerBin[k][i]->GetXaxis()->SetLabelSize(.04);
			sigmaPerBin[k][i]->GetYaxis()->SetLabelSize(.04);
            sigmaPerBin[k][i]->GetZaxis()->SetLabelSize(.05);
            
            for(int eta = 1; eta < NUM_ETA_BINS+1; eta++){
                for(int phi = 1; phi < NUM_PHI_BINS+1; phi++){
                
                   data = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinContent(eta, phi);
                   //delEtaBinCenter = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->GetBinCenter(eta);
                   //delPhiBinCenter = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->GetBinCenter(phi);
                   model = fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetBinContent(eta, phi);
                   error = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinError(eta, phi);
                   sigmaValueSquare = ((data-model)*(data-model))/(error*error); //*((data-model)/error)
                   sigmaValue = TMath::Abs(sigmaValueSquare);
                   //sigmaValue = TMath::Abs(sigmaValue);
                   //sigmaPerBin
                
                   sigmaPerBin[k][i]->SetBinContent(eta, phi, sigmaValue);
               }
            }
            
            sigmaPerBin[k][i]->SetStats(0);
            sigmaPerBin[k][i]->SetMinimum(0);
            sigmaPerBin[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            sigmaPerBin[k][i]->Draw("COLZ");  

            sigmaValueSquare = 0;
            
            for(int eta = 2; eta < 12; eta++){
                for(int phi = 1; phi < 13; phi++){
                
                   data = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinContent(eta, phi);
                   model = fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetBinContent(eta, phi);
                   error = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetBinError(eta, phi);
                   sigmaValueSquare += ((data-model)*(data-model))/(error*error); //*((data-model)/error)
                }
            }
            
         
            if(USE_SIX_PARAMETER_STANDARD) { numFitParameters = 6; }
            if(USE_TWO_2D_GAUSSIANS) { numFitParameters = 6; }
            if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD) { numFitParameters = 8; }
            if(USE_2D_SS_GAUSS_AS_GAUSS_WITH_DIPOLE_PLUS_QUAD) { numFitParameters = 7; }
            if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { numFitParameters = 7; }
           
            cout << "chi-Square for cent bin " << i << " : " << sigmaValueSquare/(35-numFitParameters) << endl;
            
            padCount++;
   
            str1 = outputFitString + ptBinLabel[k] + centralityBin[i] + fileType; 
            finalFitCorrGridCanvas->SaveAs(str1);
   
        }
        //str1 = outputFitString + ptBinLabel[k] + centralityBin[i] + fileType; 
        //finalFitCorrGridCanvas->SaveAs(str1);
    }
   
    TCanvas * finalFitCorrGridCanvasPhiProj = new TCanvas("", "", 3600, 900);
    finalFitCorrGridCanvasPhiProj->Divide(3,1);
   
    TString outputPhiProjString = "phiProjection_3_pt_bins_nominal_cuts_";
   
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        padCount = 1;
        for(int i = 0; i < 3; i++){
   
            finalFitCorrGridCanvasPhiProj->cd(padCount);
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw("EX0");
            padCount++;
        }
        str1 = outputPhiProjString + ptBinLabel[k] + fileType;
        finalFitCorrGridCanvasPhiProj->SaveAs(str1);
    }
   
   
   // finalFitCorrGridCanvas->SaveAs("Fitting_testing_for_final_results_symmeterized_April_2018.png");
    
    TCanvas * finalFitParametersCorrGridCanvas = new TCanvas("", "", 4200, 3600);
    finalFitParametersCorrGridCanvas->Divide(3,3);
    
    TString fitParameterLabel = "FIT_PARAMETERS_3_pt_bins_nominal_cuts_PtBin_";
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        
        
        if(USE_SIX_PARAMETER_STANDARD){
            finalFitParametersCorrGridCanvas->cd(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerColor(2);
            fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterOffsetPlotVsCent[k]->SetMarkerSize(2);
            fitParameterOffsetPlotVsCent[k]->SetLineWidth(2);
            fitParameterOffsetPlotVsCent[k]->SetLineColor(1);
            fitParameterOffsetPlotVsCent[k]->SetStats(0);
            fitParameterOffsetPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterDipolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterDipolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterDipolePlotVsCent[k]->SetMarkerSize(2);
            fitParameterDipolePlotVsCent[k]->SetLineWidth(2);
            fitParameterDipolePlotVsCent[k]->SetLineColor(1);
            fitParameterDipolePlotVsCent[k]->SetStats(0);
            fitParameterDipolePlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(3);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerSize(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineWidth(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineColor(1);
            fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
            fitParameterQuadrupolePlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(4);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerColor(2);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerSize(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineWidth(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineColor(1);
            fitParameterJetAmpPlotVsCent[k]->SetStats(0);
            fitParameterJetAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(5);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmaEtaPlotVsCent[k]->SetStats(0);
            fitParameterSigmaEtaPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(6);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmePhiPlotVsCent[k]->SetStats(0);
            fitParameterSigmePhiPlotVsCent[k]->Draw("E1");
       
            str1 = fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
        }
        
        if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD){
        
            finalFitParametersCorrGridCanvas->cd(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerColor(2);
            fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterOffsetPlotVsCent[k]->SetMarkerSize(2);
            fitParameterOffsetPlotVsCent[k]->SetLineWidth(2);
            fitParameterOffsetPlotVsCent[k]->SetLineColor(1);
            fitParameterOffsetPlotVsCent[k]->SetStats(0);
            fitParameterOffsetPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerSize(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineWidth(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineColor(1);
            fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
            fitParameterQuadrupolePlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(3);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerColor(2);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerSize(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineWidth(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineColor(1);
            fitParameterJetAmpPlotVsCent[k]->SetStats(0);
            fitParameterJetAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(4);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmaEtaPlotVsCent[k]->SetStats(0);
            fitParameterSigmaEtaPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(5);
             fitParameterSigmePhiPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmePhiPlotVsCent[k]->SetStats(0);
            fitParameterSigmePhiPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(6);
            fitParameterASAmpPlotVsCent[k]->SetMarkerColor(2);
            fitParameterASAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASAmpPlotVsCent[k]->SetStats(0);
            fitParameterASAmpPlotVsCent[k]->Draw();
            finalFitParametersCorrGridCanvas->cd(7);
            fitParameterASSigmaEtaPlotVsCent[k]->SetMarkerColor(2);
            fitParameterASSigmaEtaPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASSigmaEtaPlotVsCent[k]->SetStats(0);
            fitParameterASSigmaEtaPlotVsCent[k]->Draw();
            finalFitParametersCorrGridCanvas->cd(8);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerColor(2);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASSigmaPhiPlotVsCent[k]->SetStats(0);
            fitParameterASSigmaPhiPlotVsCent[k]->Draw();
            
             str1 = fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
            
        }
        
        if(USE_2D_SS_GAUSS_AS_GAUSS_WITH_DIPOLE_PLUS_QUAD) { 
        
            finalFitParametersCorrGridCanvas->cd(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerColor(2);
            fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterOffsetPlotVsCent[k]->SetMarkerSize(2);
            fitParameterOffsetPlotVsCent[k]->SetLineWidth(2);
            fitParameterOffsetPlotVsCent[k]->SetLineColor(1);
            fitParameterOffsetPlotVsCent[k]->SetStats(0);
            fitParameterOffsetPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerSize(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineWidth(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineColor(1);
            fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
            fitParameterQuadrupolePlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(3);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerColor(2);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerSize(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineWidth(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineColor(1);
            fitParameterJetAmpPlotVsCent[k]->SetStats(0);
            fitParameterJetAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(4);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmaEtaPlotVsCent[k]->SetStats(0);
            fitParameterSigmaEtaPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(5);
             fitParameterSigmePhiPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmePhiPlotVsCent[k]->SetStats(0);
            fitParameterSigmePhiPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(6);
            fitParameterDipolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterDipolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterDipolePlotVsCent[k]->SetMarkerSize(2);
            fitParameterDipolePlotVsCent[k]->SetLineWidth(2);
            fitParameterDipolePlotVsCent[k]->SetLineColor(1);
            fitParameterDipolePlotVsCent[k]->SetStats(0);
            fitParameterDipolePlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(7);
            fitParameterASSigmaEtaPlotVsCent[k]->SetMarkerColor(2);
            fitParameterASSigmaEtaPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASSigmaEtaPlotVsCent[k]->SetMarkerSize(2);
            fitParameterASSigmaEtaPlotVsCent[k]->SetLineWidth(2);
            fitParameterASSigmaEtaPlotVsCent[k]->SetLineColor(1);
            fitParameterASSigmaEtaPlotVsCent[k]->SetStats(0);
            fitParameterASSigmaEtaPlotVsCent[k]->Draw("E1");
            
            
             str1 = fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
        
        
        }
        
        if(OFFSET_QUAD_JETPEAK_AS_1D_GAUSS) { 
        
            finalFitParametersCorrGridCanvas->cd(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerColor(2);
            fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterOffsetPlotVsCent[k]->SetMarkerSize(2);
            fitParameterOffsetPlotVsCent[k]->SetLineWidth(2);
            fitParameterOffsetPlotVsCent[k]->SetLineColor(1);
            fitParameterOffsetPlotVsCent[k]->SetStats(0);
            fitParameterOffsetPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerSize(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineWidth(2);
            fitParameterQuadrupolePlotVsCent[k]->SetLineColor(1);
            fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
            fitParameterQuadrupolePlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(3);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerColor(2);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerSize(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineWidth(2);
            fitParameterJetAmpPlotVsCent[k]->SetLineColor(1);
            fitParameterJetAmpPlotVsCent[k]->SetStats(0);
            fitParameterJetAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(4);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmaEtaPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmaEtaPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmaEtaPlotVsCent[k]->SetStats(0);
            fitParameterSigmaEtaPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(5);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmePhiPlotVsCent[k]->SetMarkerSize(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineWidth(2);
            fitParameterSigmePhiPlotVsCent[k]->SetLineColor(1);
            fitParameterSigmePhiPlotVsCent[k]->SetStats(0);
            fitParameterSigmePhiPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(6);
            
            fitParameterASAmpPlotVsCent[k]->SetMarkerColor(2);
            fitParameterASAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASAmpPlotVsCent[k]->SetMarkerSize(2);
            fitParameterASAmpPlotVsCent[k]->SetLineWidth(2);
            fitParameterASAmpPlotVsCent[k]->SetLineColor(1);
            fitParameterASAmpPlotVsCent[k]->SetStats(0);
            fitParameterASAmpPlotVsCent[k]->Draw("E1");
            
            finalFitParametersCorrGridCanvas->cd(7);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerColor(2);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerSize(2);
            fitParameterASSigmaPhiPlotVsCent[k]->SetLineWidth(2);
            fitParameterASSigmaPhiPlotVsCent[k]->SetLineColor(1);
            fitParameterASSigmaPhiPlotVsCent[k]->SetStats(0);
            fitParameterASSigmaPhiPlotVsCent[k]->Draw("E1");
            
            
             str1 = fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
        
        
        }
    }
   
    double jetVolumeYield;
    double jetVolumeYieldError;
   
    TString jetVolumeYieldLabel = "Jet per-trigger yield ";
    TCanvas * jetYieldCanvas = new TCanvas("", "", 3300, 1100);
    jetYieldCanvas->Divide(3,1);
   
    if(USE_SIX_PARAMETER_STANDARD) { 
        for(int ptBin = 2; ptBin < NUM_PT_BINS; ptBin++){
        
            str1 = jetVolumeYieldLabel + ptBinLabel[ptBin];
            jetVolumePlot[ptBin] = new TH1D(str1, str1, 3, 0, 3);
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(1, "50-100%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(2, "20-50%");
            jetVolumePlot[ptBin]->GetXaxis()->SetBinLabel(3, "0-20%");
        
            
        
            for(int centBin = 0; centBin < 3; centBin++){
             
                jetYieldCanvas->cd(centBin+1);
                
                jetVolumeYield = (myfit[ptBin][centBin]->GetParameter(3)*myfit[ptBin][centBin]->GetParameter(4)*myfit[ptBin][centBin]->GetParameter(5)*nBarPrime[centBin])/2.0;
                jetVolumePlot[ptBin]->SetBinContent(centBin+1, jetVolumeYield);
                
                jetVolumeYieldError = getJetYieldError(myfit[ptBin][centBin]->GetParameter(3), myfit[ptBin][centBin]->GetParameter(4), myfit[ptBin][centBin]->GetParameter(5),
                                                       myfit[ptBin][centBin]->GetParError(3), myfit[ptBin][centBin]->GetParError(4), myfit[ptBin][centBin]->GetParError(5));
                
                //cout << jetVolumeYieldError << endl;
                
                jetVolumePlot[ptBin]->SetBinError(centBin+1, jetVolumeYieldError);
                
                cout << "# jet tracks per D0 in pt-bin: " << ptBin << " and Centrality Bin: " << centBin <<"  = " << jetVolumeYield << " +/- " << jetVolumeYieldError << endl;
            
            }  

            jetVolumePlot[ptBin]->SetMarkerStyle(20); 
            jetVolumePlot[ptBin]->Draw("p");             
        }   
        
        str1 = jetVolumeYieldLabel + fileType;
        jetYieldCanvas->SaveAs(str1);
        
    }
    
    
   
    inputFile->Close();
    outputFile->Close();
   


}


int getPtBinRange(int k){

    if(k < 3)  { return k+1; }
    if(k == 3) { return 5;   }
    if(k == 4) { return 8;   }
    
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

/*  myfit[0][0]->SetParLimits(0, -.1, .01);    //offset A0
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
    myfit[1][2]->SetParLimits(5, 0.1, 1.0);    //eta width
    ////////////////////////////////////////////////////////////////////////////////////
    //myfit[2][0]->SetParLimits(0, -.04, 0.0);    //offset A0
    //myfit[2][0]->SetParLimits(1, -.06,  0.0);    //v1     AD
    //myfit[2][0]->SetParLimits(2, 0.001,  .02);    //v2     AQ
    //myfit[2][0]->SetParLimits(3, 0.08, .16);    //A1 jet amp
    //myfit[2][0]->SetParLimits(4, 0.2, .9);     //phi width
    //myfit[2][0]->SetParLimits(5, 0.05, .5);    //eta width
	
	//myfit[2][1]->SetParLimits(0, -.04, 0.0);    //offset A0
    //myfit[2][1]->SetParLimits(1, -.03,  0.0);    //v1     AD
    //myfit[2][1]->SetParLimits(2, 0.0,  .03);    //v2     AQ
    //myfit[2][1]->SetParLimits(3, 0.02, .09);    //A1 jet amp
    //myfit[2][1]->SetParLimits(4, 0.4, .8);     //phi width
    //myfit[2][1]->SetParLimits(5, 0.6, 1.0);    //eta width
	
	//myfit[2][2]->SetParLimits(0, -.02, .01);    //offset A0
    //myfit[2][2]->SetParLimits(1, -.02,  0.0);    //v1     AD
    //myfit[2][2]->SetParLimits(2, 0.0,  .02);    //v2     AQ
    //myfit[2][2]->SetParLimits(3, 0.0, .05);    //A1 jet amp
    //myfit[2][2]->SetParLimits(4, 0.3, .8);     //phi width
    //myfit[2][2]->SetParLimits(5, 0.3, .8);    //eta width  */
    
    
    
    
    /*2D Gauss plus Quad
    //////////////////////////////////////////////////////////////////////////////////////
        myfit[3][0]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[3][0]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[3][0]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[3][0]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[3][0]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[3][0]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        
         //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][1]->SetParameters(-.005, -.01,   .01,  .03,  .7, 1.4);
        
        myfit[3][1]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[3][1]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[3][1]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[3][1]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[3][1]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[3][1]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][2]->SetParameters(-.0007, -.002, .003, .02,  .32, .37);
        
        myfit[3][2]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[3][2]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[3][2]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[3][2]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[3][2]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[3][2]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        ///////////////////////////////////////////////////////////////////////////////////////
        myfit[4][0]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[4][0]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[4][0]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[4][0]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[4][0]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[4][0]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        
         //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][1]->SetParameters(-.005, -.01,   .01,  .03,  .7, 1.4);
        
        myfit[4][1]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[4][1]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[4][1]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[4][1]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[4][1]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[4][1]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][2]->SetParameters(-.0007, -.002, .003, .02,  .32, .37);
        
        myfit[4][2]->SetParLimits(0, 0.0, .15);    //same-side amp
        myfit[4][2]->SetParLimits(1, 0.0, 2.0);    //same-side phi width
        myfit[4][2]->SetParLimits(2, 0.0,  2.0);    //same-side eta width
        myfit[4][2]->SetParLimits(3, 0.0, .04);    //away-side amp
        myfit[4][2]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        myfit[4][2]->SetParLimits(5, 0.0, 2.0);    //away-side eta width
        */