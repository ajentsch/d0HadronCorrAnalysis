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

int D0PhiProjectionFitter(){


    bool USE_SIX_PARAMETER_STANDARD     = false;
    bool USE_TWO_2D_GAUSSIANS           = false;
    bool USE_TWO_2D_GAUSSIANS_PLUS_QUAD = false;
    
    bool FIT_NON_JET_COMPONENT_WITH_ETA_GAP = true;
    
    
    //////////////////////////////////////////////
    
    bool USE_ETA_GAP = true;
    int ETA_GAP_BIN_LEFT_LOW = 3;
    int ETA_GAP_BIN_LEFT_HIGH = 4;
    int ETA_GAP_BIN_RIGHT_LOW = 9;
    int ETA_GAP_BIN_RIGHT_HIGH = 11;
    int ETA_RANGE_LOW  = 2;
    int ETA_RANGE_HIGH = 11;


    //output strings
    TString outputStringRawCorr = "12phi_13eta_no_DStar_pion_PID.png";
    TString outputFitString = "Fitting_12phi_13eta__no_DStar_pion_PID_PtBin_";
    outputFitString = "12phi_13eta_no_DStar_pion_PID_Fits_Superposed_onto_data_";
    TString fitParameterLabel = "FIT_PARAMETERS_12phi_13eta_no_DStar_pion_PID_PtBin_";
    

    double NDF = 1;
    
    double projectionScaleFactor;
    
    if(USE_ETA_GAP) { projectionScaleFactor = 1.0/((ETA_GAP_BIN_LEFT_HIGH-ETA_GAP_BIN_LEFT_LOW)+1); }
    else projectionScaleFactor = 1.0/((ETA_RANGE_HIGH-ETA_RANGE_LOW)+1);
    
    cout << "Scale Factor: " << projectionScaleFactor << endl;

    int NUM_PT_BINS = 3;
    int FIRST_PT_BIN = 2;
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 13;
    double ETA_RANGE = 1.69;
    double phiBinShift = (TMath::Pi()/12.0);
    
    if(USE_SIX_PARAMETER_STANDARD) { NDF = NUM_PHI_BINS - 5.0; }
    if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD) { NDF = NUM_PHI_BINS - 6.0; }
    
    double v2Array[5][3];
    double v2ArrayStatError[5][3];

    TString inputFileName  = "d0HadronCorrMaker_OUTPUT_dStar_correction_with_no_pion_PID_12_phi_13_eta_2_10_pt.root";
	TString outputFileName = "fittingOutput_OUTPUT_dStar_correction_with_no_pion_PID_12_phi_13_eta_2_10_pt.root";
    TString outputFolder   = "deltaPhiProjectionFits/";
    TString path = "C:/Users/ajentsch/desktop/D0_Analysis_Code_Local_Copies/D0CorrelationFitter/";
	
	TString fullySubtractedLabelCustomCentrality = "Full_DStar_Correction_FullSubtractedCorr_SideBand__PtBin_";
	TString ptBinLabel[5] = {"0", "1", "2", "3", "4"};
    TString centralityBin[3] = {"Peripheral", "MidCentral", "Central"};
    TString fileType = ".png";
    
    TString ptLabel     = "PtBin_";
    TString dipoleLabel = "Dipole_AD_";
    
    TString offsetLabel = "Offset_A0_";
    TString quadLabel   = "Quadrupole_AQ_";
    
    TString jetAmpLabel = "JetAmp_A1_";
    TString sigEtaLabel = "SigmaEta_";
    TString sigPhiLabel = "SigmaPhi_";
    
    TString awaySideAmpLabel = "AwaySideAmp_";
    TString awaySidePhiWidthLabel = "AwaySide_Phi_Width_";
    
    TString leftLabel = "_left_";
    TString rightLabel = "_right_";
	
	TString phiProj = "Phi_projection_";
	TString str1;
   
	TH2D* fullySubtractedCorrCustomCentralitySideBand[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProj[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProjLeft[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProjRight[5][3];
	TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[5][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[5][3];
   
    TH1D* fitParameterOffsetPlotVsCent[5];
    TH1D* fitParameterDipolePlotVsCent[5];
    TH1D* fitParameterQuadrupolePlotVsCent[5];
    TH1D* fitParameterJetAmpPlotVsCent[5];
    TH1D* fitParameterSigmaPhiPlotVsCent[5];
    TH1D* fitParameterASAmpPlotVsCent[5];
    TH1D* fitParameterASSigmaPhiPlotVsCent[5];
    
   // TH2D* errorOnBin[3];
    TH1D* sigmaPerBin[5][3];
    
   
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
    
        str1 = offsetLabel + ptLabel + ptBinLabel[k];
        fitParameterOffsetPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = dipoleLabel + ptLabel + ptBinLabel[k];
        fitParameterDipolePlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = quadLabel + ptLabel + ptBinLabel[k];
        fitParameterQuadrupolePlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = jetAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterJetAmpPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = sigPhiLabel + ptLabel + ptBinLabel[k];
        fitParameterSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        
        str1 = awaySideAmpLabel + ptLabel + ptBinLabel[k];
        fitParameterASAmpPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        str1 = awaySidePhiWidthLabel + ptLabel + ptBinLabel[k];
        fitParameterASSigmaPhiPlotVsCent[k] = new TH1D(str1, str1, 3, 0, 3);
        
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterOffsetPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterDipolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterQuadrupolePlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterJetAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterSigmaPhiPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
        
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(1, "50-100%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(2, "20-50%");
        fitParameterASAmpPlotVsCent[k]->GetXaxis()->SetBinLabel(3, "0-20%");
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
    
    TFile * inputFile  = new TFile(inputFileName);
    TFile * outputFile = new TFile(outputFileName, "RECREATE");
 
 
    TCanvas *c = new TCanvas("c3", "Histograms", 1100, 850);
 
    TString fitFunctionLabel = "parameterFitFunction";
 
    //TF2 *myfit   = new TF2("parameterFitFunction","[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
    
    TF1 *myfit[5][3];   
    
    TF1 *myfit1D = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x)" , -1.57, 4.71);
    
    TF1 *fitProjection = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x) + [3]*[4]*(exp(-0.5*(x*x)/([5]*[5]))+exp(-0.5*((x-6.28)*(x-6.28))/([5]*[5])))" , -1.57, 4.71);
    
    TF1 *quadrupole    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TF1 *dipole        = new TF1("dipole", "[0]*cos(x)", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TF1 *offset        = new TF1("Offset", "[0]", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TF1 *SSGauss       = new TF1("SSGauss", "[0]*((exp(-0.5*((x-6.28)*(x-6.28))/([1]*[1])))+exp(-0.5*(x*x)/([1]*[1])))", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TF1 *ASGauss       = new TF1("ASGauss", "[0]*(exp(-0.5*((x-3.14159)*(x-3.14159))/([1]*[1]))+exp(-0.5*((x+3.14159)*(x+3.14159))/([1]*[1])))", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TF1 *zeroLine      = new TF1("zeroLine", "0", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    
    TString fitLabel = "_Fit";
    TString resLabel = "_Residual";
    
    /* ofstream fitParameterFile1D;
    TString fitParameterFileName1D = "Extracted_Parameters_from_fitting_1D.txt";
    TString fitParamterOutputFile1D = path + fitParameterFileName1D;
    fitParameterFile1D.open (fitParamterOutputFile1D);*/
    
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
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){ 
 
            str1 = fitFunctionLabel + ptBinLabel[k] + centralityBin[i];
            if(USE_SIX_PARAMETER_STANDARD) { myfit[k][i] = new TF1(str1,"[0] + [1]*cos(x) + [2]*2*cos(2*x) + [3]*((exp(-0.5*((x-6.28)*(x-6.28))/([4]*[4])))+exp(-0.5*(x*x)/([4]*[4])))", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
            if(USE_TWO_2D_GAUSSIANS) { myfit[k][i] = new TF1(str1,"[0] + [1]*(exp(-0.5*(x*x)/([2]*[2]))+exp(-0.5*((x-6.28)*(x-6.28))/([2]*[2]))) + [3]*(exp(-0.5*((x-3.14159)*(x-3.14159))/([4]*[4]))+exp(-0.5*((x+3.14159)*(x+3.14159))/([4]*[4])))", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
            if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD) { myfit[k][i] = new TF1(str1,"[0] + 2*[1]*cos(2*x) + [2]*(exp(-0.5*(x*x)/([3]*[3])) + exp(-0.5*((x-6.28)*(x-6.28))/([3]*[3]))) + [4]*(exp(-0.5*((x-3.14159)*(x-3.14159))/([5]*[5]))+exp(-0.5*((x+3.14159)*(x+3.14159))/([5]*[5])))", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
            if(FIT_NON_JET_COMPONENT_WITH_ETA_GAP) { myfit[k][i] = new TF1(str1, "[0] + [1]*cos(x) + 2*[2]*cos(2*x)", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);}
        
        }
	}
    
   
    if(USE_SIX_PARAMETER_STANDARD) {
    //                               A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][0]->SetParameters(-.003,  -.02,   .02,   .08,  .3, .4);
    
        myfit[2][0]->SetParLimits(0, -.01, 0.0);    //offset A0
        myfit[2][0]->SetParLimits(1, -.01,  0.0);    //v1     AD
        //myfit[2][0]->SetParLimits(2, 0.0,  .011);    //v2     AQ
        myfit[2][0]->FixParameter(2, .00938519);    //AQ
        myfit[2][0]->SetParLimits(3, 0.0, 0.02);    //A1 jet amp
        myfit[2][0]->SetParLimits(4, 0.75, 1.1);     //phi width
        
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][1]->SetParameters(-.005, -.01,   .01,  .03,  .7, 1.4);
        
        myfit[2][1]->SetParLimits(0, -.1, 0.0);    //offset A0
        myfit[2][1]->SetParLimits(1, -.025,  0.0);    //v1     AD
        myfit[2][1]->SetParLimits(2, 0.0,  .01);    //v2     AQ
        //myfit[2][1]->FixParameter(2, .00639478);    //AQ
        myfit[2][1]->SetParLimits(3, 0.0, .2);    //A1 jet amp
        myfit[2][1]->SetParLimits(4, 0.4, 1.2);     //phi width
       
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][2]->SetParameters(-.0007, -.002, .003, .02,  .32, .37);
        
        myfit[2][2]->SetParLimits(0, -.02, 0.0);    //offset A0
        myfit[2][2]->SetParLimits(1, -.03,  0.0);    //v1     AD
        myfit[2][2]->SetParLimits(2, 0.0,  .05);    //v2     AQ
        //myfit[2][2]->FixParameter(2, .00211007);    //AQ
        myfit[2][2]->SetParLimits(3, 0.0, .07);    //A1 jet amp
        myfit[2][2]->SetParLimits(4, 0.0, 1.2);     //phi width
    }
    
    if(USE_TWO_2D_GAUSSIANS) {
    //                               A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][0]->SetParameters(-.003,  -.02,   .02,   .08,  .3, .4);
        myfit[2][0]->SetParLimits(0, -0.05, -.04);
        myfit[2][0]->SetParLimits(1, 0.04, .09);    //same-side amp
        myfit[2][0]->SetParLimits(2, 0.2, .8);    //same-side phi width
        myfit[2][0]->SetParLimits(3, 0.02, .09);    //away-side amp
        myfit[2][0]->SetParLimits(4, 0.2, .8);     //away-side phi width
       
        
         //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][1]->SetParameters(-.005, -.01,   .01,  .03,  .7, 1.4);
        myfit[2][1]->SetParLimits(0, -0.03, 0.0);
        myfit[2][1]->SetParLimits(1, 0.0, .15);    //same-side amp
        myfit[2][1]->SetParLimits(2, 0.0, 2.0);    //same-side phi width
        myfit[2][1]->SetParLimits(3, 0.0, .1);    //away-side amp
        myfit[2][1]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
        
        //                         A0      AD     AQ    A1   sPhi sEta  
        //myfit[2][2]->SetParameters(-.0007, -.002, .003, .02,  .32, .37);
        myfit[2][2]->SetParLimits(0, -0.03, 0.0);
        myfit[2][2]->SetParLimits(1, 0.0, .15);    //same-side amp
        myfit[2][2]->SetParLimits(2, 0.0, 2.0);    //same-side phi width
        myfit[2][2]->SetParLimits(3, 0.0, .1);    //away-side amp
        myfit[2][2]->SetParLimits(4, 0.0, 2.0);     //away-side phi width
   
     
    }
    
    if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD) {
    //                              A0      AQ     A1    JetPhi    ASAmp   ASPhi    
        //myfit[2][0]->SetParameters(-.000187,    .008427,   .3,    .9,       .4,     .6);
    
        myfit[2][0]->SetParLimits(0, -.035, -.01);    //offset 
        myfit[2][0]->FixParameter(1, .00938519);  
        //myfit[2][0]->SetParLimits(1, .0045, .0055);    //quadrupole  
        myfit[2][0]->SetParLimits(2, 0.0,  .1);   //jet amp
        myfit[2][0]->SetParLimits(3, 0.4, 1.5);    //jet phi width
        myfit[2][0]->SetParLimits(4, 0.0, .09);    //AS amp
        myfit[2][0]->SetParLimits(5, .7, 1.2);    //AS phi width
        
        //Quadrupole from 2D fit - .00938519
        
        //                          A0      AQ    A1    JetPhi     ASAmp   ASPhi   
       //myfit[2][1]->SetParameters(-.003,  .01,   .8,    .6,         .02,     .9);
        
        myfit[2][1]->SetParLimits(0, -.03, -.01);    //offset  
        myfit[2][1]->FixParameter(1, .00639478);
        //myfit[2][1]->SetParLimits(1, 0.005, .015);    //quadrupole  
        myfit[2][1]->SetParLimits(2, 0.012,  .05);   //jet amp
        myfit[2][1]->SetParLimits(3, 0.5, 1.2);    //jet phi width
        myfit[2][1]->SetParLimits(4, 0.006, .05);    //AS amp
        myfit[2][1]->SetParLimits(5, .8, 1.7);    //AS phi width
        
        //quadrupole from 2D fit - .00639478
        
        //                         A0      AQ     JetAmp JetWidth    ASAMP   ASWIDTH
        //myfit[2][2]->SetParameters(-.001, .0035,  .0023,    .4,      .0012,    2.0);
        
        myfit[2][2]->SetParLimits(0, -.01, 0.0);    //offset  
        myfit[2][2]->FixParameter(1, .00211007);
        //myfit[2][2]->SetParLimits(1, 0.003, .0036);    //quadrupole  
        myfit[2][2]->SetParLimits(2, 0.001,  .02);   //jet amp
        myfit[2][2]->SetParLimits(3, 0.5, 1.0);    //jet phi width
        myfit[2][2]->SetParLimits(4, 0.001, .02);    //AS amp
        myfit[2][2]->SetParLimits(5, 1.0, 1.6);    //AS phi width
        
        //quadrupole from 2D fit - .00211007
    }
 
    if(FIT_NON_JET_COMPONENT_WITH_ETA_GAP){ 
 
        //myfit[2][0]->SetParLimits(0, -.4, 0.0);    //offset
        //myfit[2][0]->SetParLimits(0, 0.0045, .0055);    //quadrupole  
       
        //myfit[2][1]->SetParLimits(0, -.2, 0.0);    //offset  
        //myfit[2][1]->SetParLimits(0, 0.0065, .008);    //quadrupole  
        
        //myfit[2][2]->SetParLimits(0, -.07, 0.0);    //offset  
        //myfit[2][2]->SetParLimits(0, 0.0034, .0038);    //quadrupole  
    }
 
    double projScale;
 
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){ 
 
        
            str1 = fullySubtractedLabelCustomCentrality  + ptBinLabel[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBand[k][i] = (TH2D*)inputFile->Get(str1);
            fullySubtractedCorrCustomCentralitySideBand[k][i]->SetDirectory(0);
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality  + ptBinLabel[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality + leftLabel + ptBinLabel[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBandPhiProjLeft[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 = phiProj + fullySubtractedLabelCustomCentrality + rightLabel + ptBinLabel[k] + centralityBin[i];
            fullySubtractedCorrCustomCentralitySideBandPhiProjRight[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
            //SET THE ETA RANGE HERE
           
            if(USE_ETA_GAP) { 
            
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY("",ETA_GAP_BIN_LEFT_LOW, ETA_GAP_BIN_LEFT_HIGH,"");
                //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY("",ETA_GAP_BIN_RIGHT_LOW, ETA_GAP_BIN_RIGHT_HIGH,"");
                //fullySubtractedCorrCustomCentralitySideBandPhiProjLeft[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY("",ETA_GAP_BIN_LEFT_LOW, ETA_GAP_BIN_LEFT_HIGH,"");
                //fullySubtractedCorrCustomCentralitySideBandPhiProjRight[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY("",ETA_GAP_BIN_RIGHT_LOW, ETA_GAP_BIN_RIGHT_HIGH,"");
                //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Add(fullySubtractedCorrCustomCentralitySideBandPhiProjLeft[k][i], fullySubtractedCorrCustomCentralitySideBandPhiProjRight[k][i], .5, .5);
            
            }
            else fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY("",3,11,"");
            
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetNameTitle(str1, str1);
           
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Scale(projectionScaleFactor);
            
            //cout << "code gets here: " << endl;
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            str1 = phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
 
            
            /////////////////////////////1d fitting///////////////////////////////////////////
    
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Fit(myfit[k][i], "R0");
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hfit");
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Eval(myfit[k][i]);
 
           /* fitParameterFile1D << endl << "pt bin (3 is integrated): " << k << "    Centrality Bin (0 - per, 1 - mid cent, 2 - cent): " << i << endl;
            fitParameterFile1D << endl << "A0" << "\t\t" << "AD" << "\t\t" <<"AQ" << endl; // "\t\t" << "Jet Amp A1" << "\t\t" << "SigPhi" << endl;
            
            fitParameterFile1D << endl << myfit1D->GetParameter(0) << "\t" << myfit1D->GetParameter(1) << "\t" << myfit1D->GetParameter(2) << endl; //<< "\t" << myfit1D->GetParameter(3) << "\t" << myfit1D->GetParameter(5) << endl; 
            fitParameterFile1D << myfit1D->GetParError(0) << "\t" << myfit1D->GetParError(1) << "\t" << myfit1D->GetParError(2) << endl;//<< "\t" << myfit1D->GetParError(3) << "\t" << myfit1D->GetParError(5) << endl; 
            
            v2_2D = TMath::Sqrt(myfit[k][i]->GetParameter(2));
            v2_2D_Error = v2_2D*(.5)*(myfit[k][i]->GetParError(2)/myfit[k][i]->GetParameter(2));
            
            v2Array[k][i] = v2_2D;
            v2ArrayStatError[k][i] = v2_2D_Error;
            
            
            v2_1D = TMath::Sqrt(myfit1D->GetParameter(2));
            v2_1D_Error = v2_2D*(.5)*(myfit1D->GetParError(2)/myfit1D->GetParameter(2));*/
            
            TF1 *quadrupole    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            TF1 *dipole        = new TF1("dipole", "[0]*cos(x)", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            TF1 *offset        = new TF1("Offset", "[0]", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            TF1 *SSGauss        = new TF1("SSGauss", "[0]*((exp(-0.5*((x-6.28)*(x-6.28))/([1]*[1])))+exp(-0.5*(x*x)/([1]*[1])))", -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
            
            
            //quadrupole1->SetParameter(0, myfit[k][i]->GetParameter(2));
            //quadrupole2->SetParameter(0, myfit1D->GetParameter(2));
            //quadrupole3->SetParameter(0, fitProjection->GetParameter(2));
            
            //fitProjection->SetLineColor(2); //red
            //quadrupole1->SetLineColor(3); // Green, from 2D fit
            //quadrupole2->SetLineColor(4); //blue, from 1D fit
            //myfit1D->SetLineColor(7); // light blue
            //quadrupole3->SetLineColor(41);
            
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Add(myfit[k][i], -1);
            
            //str1 = path + outputFolder + phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + fitLabel + fileType;
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Draw();
            //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw("EX0");
            //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetMarkerStyle(20);
            
            //fitProjection->Draw("SAME");
            //quadrupole1->Draw("SAME");
            //quadrupole2->Draw("SAME");
            //quadrupole3->Draw("SAME");
            //myfit1D->Draw("SAME");
            //c->SaveAs(str1);
 
            //str1 = path + outputFolder + phiProj + fullySubtractedLabelCustomCentrality + ptBinLabel[k] + centralityBin[i] + resLabel + fileType;
            //fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Draw();
            //c->SaveAs(str1); 
 
            ptBinRange = k+1;
 
            if(k == 0) {continue;}
 
        }
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
                fitParameterSigmaPhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(4));
                    
                fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterDipolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
                fitParameterQuadrupolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(2));
                fitParameterJetAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(3));
                fitParameterSigmaPhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(4));
            }
            
            if(USE_TWO_2D_GAUSSIANS){
            
                fitParameterOffsetPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(0));
                fitParameterJetAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(1));      //jet amp
                fitParameterSigmaPhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(2));    // jet phi
                fitParameterASAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(3));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(4));  //AS phi
                
                fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterJetAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
                fitParameterSigmaPhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(2));
                fitParameterASAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(3));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(4));  //AS phi
                
            
            }
            
            if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD){
            
                fitParameterOffsetPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(0));      //offset
                fitParameterQuadrupolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(1));  //quad
                fitParameterJetAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(2));      //jet amp
                fitParameterSigmaPhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(3));    // jet phi
                fitParameterASAmpPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(4));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(5));  //AS phi
                
                fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterQuadrupolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
                fitParameterJetAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(2));
                fitParameterSigmaPhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(3));
                fitParameterASAmpPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(4));       //As amp
                fitParameterASSigmaPhiPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(5));  //AS phi
                
            
            }
            
            if(FIT_NON_JET_COMPONENT_WITH_ETA_GAP){ 
 
                //fitParameterOffsetPlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(0));      //offset
                fitParameterQuadrupolePlotVsCent[k]->SetBinContent(i+1, myfit[k][i]->GetParameter(1));  //quad
               
                //fitParameterOffsetPlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(0));
                fitParameterQuadrupolePlotVsCent[k]->SetBinError(i+1, myfit[k][i]->GetParError(1));
               
            }
        }
    }        
 
    //fitParameterFile2D.close();
    //fitParameterFile1D.close();
 
    TCanvas * finalCorrGridCanvas = new TCanvas("", "", 3600, 3400);
    finalCorrGridCanvas->Divide(3,3);
    
    TString ptBinTitle[5] = {"0 GeV/c < p_{t, D^{0}} < 1 GeV/c", "1 GeV/c < p_{t, D^{0}} < 2 GeV/c", "2 GeV/c < p_{t, D^{0}} < 10 GeV/c",   //change made from 2-3 to 2-10
                             "3 GeV/c < p_{t, D^{0}} < 4 GeV/c", "4 GeV/c < p_{t, D^{0}} < 10 GeV/c"}; 
    int padCount = 1;
    for(int i = 0; i < 3; i++){
        for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        
        
            finalCorrGridCanvas->cd(padCount);
            
            if(i == 0) { fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetTitle(ptBinTitle[k]); }
            else fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetTitle("");
            
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetMarkerStyle(20);
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetMarkerColor(2);
            
            double scaleFactor = 1.0;
            
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Scale(scaleFactor);
            
            //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetXaxis()->SetTitle("#Delta#Phi");
            //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetXaxis()->CenterTitle();
            //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetXaxis()->SetTitleOffset(1.5);
            //fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetXaxis()->SetTitleSize(.065);
            
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw();
            
            padCount++;
   
        }
    }
    
    
    str1 = path + outputFolder + outputStringRawCorr;
   
    finalCorrGridCanvas->SaveAs(str1);
    
    TCanvas * finalFitCorrGridCanvas = new TCanvas("", "", 4800, 2700);
    finalFitCorrGridCanvas->Divide(4,3);
    
    
    
    
    
    double sigmaValueSquare, sigmaValue;
    double data, model, error;
    double delEtaBinCenter, delPhiBinCenter;
    double maxZ, minZ;
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        padCount = 1;
        for(int i = 0; i < 3; i++){
        
        
           //cout << "CHECK IT!" << endl;
        
            finalFitCorrGridCanvas->cd(padCount);
            
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw();
            
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->SetMaximum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMaximum());
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->SetMinimum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMinimum());
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->SetStats(0);
            
            padCount++;
            finalFitCorrGridCanvas->cd(padCount);
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->SetMaximum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMaximum());
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->SetMinimum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMinimum());
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->SetStats(0);
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetLineColor(1);
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i]->SetLineWidth(1);
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Draw();
            padCount++;
            finalFitCorrGridCanvas->cd(padCount);
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->SetMaximum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMaximum());
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->SetMinimum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMinimum());
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Draw();
            padCount++;
            finalFitCorrGridCanvas->cd(padCount);
            
            sigmaPerBin[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Clone();
            
            for(int phi = 1; phi < NUM_PHI_BINS+1; phi++){
                
               data = fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetBinContent(phi);
               //delEtaBinCenter = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->GetBinCenter(eta);
               //delPhiBinCenter = fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->GetBinCenter(phi);
               model = fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->GetBinContent(phi);
               error = fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetBinError(phi);
               sigmaValueSquare = ((data-model)*(data-model))/(error*error); //*((data-model)/error)
               sigmaValue = TMath::Abs(sigmaValueSquare);
               //sigmaValue = TMath::Abs(sigmaValue);
               //sigmaPerBin
            
               sigmaPerBin[k][i]->SetBinContent(phi, sigmaValueSquare);
               
            }
            
            sigmaPerBin[k][i]->SetStats(0);
            sigmaPerBin[k][i]->SetMinimum(0);
            //errorOnBin[i]->SetMaximum(.03);
            sigmaPerBin[k][i]->Draw();            
            
            padCount++;
   
        }
        str1 = path + outputFolder + outputFitString + ptBinLabel[k] + fileType; 
        finalFitCorrGridCanvas->SaveAs(str1);
    }
   
    
    TCanvas * finalFitCorrGridCanvasFitOnData = new TCanvas("", "", 3600, 1100);
    finalFitCorrGridCanvasFitOnData->Divide(3,1);
   
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        
        for(int i = 0; i < 3; i++){
        
            padCount = 1;
        
           //cout << "CHECK IT!" << endl;
        
            finalFitCorrGridCanvasFitOnData->cd(padCount);
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw();
            padCount++;
            finalFitCorrGridCanvasFitOnData->cd(padCount);
            
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw();
            myfit[k][i]->SetLineColor(4);
            myfit[k][i]->SetLineWidth(3);
            myfit[k][i]->Draw("SAME");
            
            if(USE_SIX_PARAMETER_STANDARD){
                
                offset->SetParameter(0, myfit[k][i]->GetParameter(0));
                dipole->SetParameter(0, myfit[k][i]->GetParameter(1));
                quadrupole->SetParameter(0, myfit[k][i]->GetParameter(2));
                SSGauss->SetParameter(0, myfit[k][i]->GetParameter(3));
                SSGauss->SetParameter(1, myfit[k][i]->GetParameter(4));
                offset->SetLineColor(2); //red
                dipole->SetLineColor(3); //green
                quadrupole->SetLineColor(6); //magenta
                SSGauss->SetLineColor(7);   //light blue
                //offset->Draw("SAME");//red
                //dipole->Draw("SAME"); //green
                quadrupole->Draw("SAME"); //??????
                //SSGauss->Draw("SAME");   //light blue
            }
            
            if(USE_TWO_2D_GAUSSIANS){
            
                SSGauss->SetParameter(0, myfit[k][i]->GetParameter(1));
                SSGauss->SetParameter(1, myfit[k][i]->GetParameter(2));
                ASGauss->SetParameter(0, myfit[k][i]->GetParameter(3));
                ASGauss->SetParameter(1, myfit[k][i]->GetParameter(4));
                //offset->SetLineColor(2); //red
                quadrupole->SetLineColor(6); //magenta
                //SSGauss->SetLineColor(7);   //light blue
                //ASGauss->SetLineColor(3);   //green
                //offset->Draw("SAME"); //red
                quadrupole->Draw("SAME"); //??????
                //SSGauss->Draw("SAME");   //light blue
                //ASGauss->Draw("SAME");   //light blue
            }
            
            if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD){
                
                offset->SetParameter(0, myfit[k][i]->GetParameter(0));
                quadrupole->SetParameter(0, myfit[k][i]->GetParameter(1));
                SSGauss->SetParameter(0, myfit[k][i]->GetParameter(2));
                SSGauss->SetParameter(1, myfit[k][i]->GetParameter(3));
                ASGauss->SetParameter(0, myfit[k][i]->GetParameter(4));
                ASGauss->SetParameter(1, myfit[k][i]->GetParameter(5));
                offset->SetLineColor(2); //red
                quadrupole->SetLineColor(6); //magenta
                SSGauss->SetLineColor(7);   //light blue
                ASGauss->SetLineColor(3);   //green
                //offset->Draw("SAME"); //red
                quadrupole->Draw("SAME"); //??????
                //SSGauss->Draw("SAME");   //light blue
                //ASGauss->Draw("SAME");   //light blue
                
                //dipole->SetParameter(0, -.02);
                //dipole->SetLineColor(3); //green
                //dipole->Draw("SAME"); //green
            }
            padCount++;
            finalFitCorrGridCanvasFitOnData->cd(padCount);
            
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->SetMaximum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMaximum());
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->SetMinimum(fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->GetMinimum());
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->SetStats(0);
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Draw();
            zeroLine->Draw("SAME");
            
            padCount++;
   
            str1 = path + outputFolder + outputFitString + centralityBin[i] + ptBinLabel[k] + fileType; 
            finalFitCorrGridCanvasFitOnData->SaveAs(str1);
   
        }
        
    }
   
   
    // finalFitCorrGridCanvas->SaveAs("Fitting_testing_for_final_results_symmeterized_April_2018.png");
    
    TCanvas * finalFitParametersCorrGridCanvas = new TCanvas("", "", 3600, 2000);
    finalFitParametersCorrGridCanvas->Divide(3,2);
    
    
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){
        
        
        if(USE_SIX_PARAMETER_STANDARD){
            finalFitParametersCorrGridCanvas->cd(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerColor(2);
            fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterOffsetPlotVsCent[k]->SetStats(0);
            fitParameterOffsetPlotVsCent[k]->Draw();
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterDipolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterDipolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterDipolePlotVsCent[k]->SetStats(0);
            fitParameterDipolePlotVsCent[k]->Draw();
            finalFitParametersCorrGridCanvas->cd(3);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
            fitParameterQuadrupolePlotVsCent[k]->Draw();
            finalFitParametersCorrGridCanvas->cd(4);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerColor(2);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterJetAmpPlotVsCent[k]->SetStats(0);
            fitParameterJetAmpPlotVsCent[k]->Draw();
            finalFitParametersCorrGridCanvas->cd(5);
            fitParameterSigmaPhiPlotVsCent[k]->SetMarkerColor(2);
            fitParameterSigmaPhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmaPhiPlotVsCent[k]->SetStats(0);
            fitParameterSigmaPhiPlotVsCent[k]->Draw();
       
            str1 = path + outputFolder +  fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
        }
        
        if(USE_TWO_2D_GAUSSIANS){
        
            finalFitParametersCorrGridCanvas->cd(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerColor(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterOffsetPlotVsCent[k]->SetStats(0);
            fitParameterOffsetPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerColor(1);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterJetAmpPlotVsCent[k]->SetStats(0);
            fitParameterJetAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(3);
            fitParameterSigmaPhiPlotVsCent[k]->SetMarkerColor(1);
            fitParameterSigmaPhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmaPhiPlotVsCent[k]->SetStats(0);
            fitParameterSigmaPhiPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(4);
            fitParameterASAmpPlotVsCent[k]->SetMarkerColor(1);
            fitParameterASAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASAmpPlotVsCent[k]->SetStats(0);
            fitParameterASAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(5);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerColor(1);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASSigmaPhiPlotVsCent[k]->SetStats(0);
            fitParameterASSigmaPhiPlotVsCent[k]->Draw("E1");
            
            str1 = path + outputFolder + fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
        }
        
        if(USE_TWO_2D_GAUSSIANS_PLUS_QUAD){
        
            finalFitParametersCorrGridCanvas->cd(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerColor(1);
            fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterOffsetPlotVsCent[k]->SetStats(0);
            fitParameterOffsetPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerColor(1);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
            fitParameterQuadrupolePlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(3);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerColor(1);
            fitParameterJetAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterJetAmpPlotVsCent[k]->SetStats(0);
            fitParameterJetAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(4);
            fitParameterSigmaPhiPlotVsCent[k]->SetMarkerColor(1);
            fitParameterSigmaPhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterSigmaPhiPlotVsCent[k]->SetStats(0);
            fitParameterSigmaPhiPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(5);
            fitParameterASAmpPlotVsCent[k]->SetMarkerColor(1);
            fitParameterASAmpPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASAmpPlotVsCent[k]->SetStats(0);
            fitParameterASAmpPlotVsCent[k]->Draw("E1");
            finalFitParametersCorrGridCanvas->cd(6);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerColor(1);
            fitParameterASSigmaPhiPlotVsCent[k]->SetMarkerStyle(8);
            fitParameterASSigmaPhiPlotVsCent[k]->SetStats(0);
            fitParameterASSigmaPhiPlotVsCent[k]->Draw("E1");
            
            str1 = path + outputFolder + fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
            
        }
        
         if(FIT_NON_JET_COMPONENT_WITH_ETA_GAP){ 
         
            //finalFitParametersCorrGridCanvas->cd(1);
            //fitParameterOffsetPlotVsCent[k]->SetMarkerColor(2);
            //fitParameterOffsetPlotVsCent[k]->SetMarkerStyle(8);
            //fitParameterOffsetPlotVsCent[k]->SetStats(0);
            //fitParameterOffsetPlotVsCent[k]->Draw();
            finalFitParametersCorrGridCanvas->cd(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerColor(2);
            fitParameterQuadrupolePlotVsCent[k]->SetMarkerStyle(8);
            fitParameterQuadrupolePlotVsCent[k]->SetStats(0);
            fitParameterQuadrupolePlotVsCent[k]->Draw();
           
            
            str1 = path + outputFolder + fitParameterLabel + ptBinLabel[k] + fileType;
       
            finalFitParametersCorrGridCanvas->SaveAs(str1);
        }
    }

    double XSquarePDF;
    
    for(int k = FIRST_PT_BIN; k < NUM_PT_BINS; k++){    
        for(int i = 0; i < 3; i++) { 
        
            cout << "Chi-Square = " << myfit[k][i]->GetChisquare() << endl;
            XSquarePDF = myfit[k][i]->GetChisquare()/NDF;
            cout << "chi-square p.d.f. for cent bin " << i << "  = " << XSquarePDF << endl;
        }
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


