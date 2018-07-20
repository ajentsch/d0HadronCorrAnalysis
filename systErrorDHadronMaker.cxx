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

#include "C:/Users/ajentsch/desktop/systErrorDHadronMaker.h"

int systErrorDHadronMaker(){

    
    int NUM_ETA_BINS = 9;
    int NUM_PHI_BINS = 12;
    
    TH2D *secondariesErrors            = new TH2D("secondariesErrors","secondariesErrors", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *bFeedDownErrors              = new TH2D("bFeedDownErrors","bFeedDownErrors", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    //TH2D *pileupErrors               = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *fBarFactorErrors             = new TH2D("fBarFactorErrors","fBarFactorErrors", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *topologyCutErrors            = new TH2D("topologyCutErrors","topologyCutErrors", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *backgroundSubtractionErrors  = new TH2D("backgroundSubtractionErrors","backgroundSubtractionErrors", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    
    TH2D *finalErrorResults            = new TH2D("systematicErrors","systematicErrors", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *subtract                     = new TH2D("subtract","subtract", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    
    
    TH2D *tmp             = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *tmp1            = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *tmp2            = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *tmp3            = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *tmp4            = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *tmp5            = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *nominalRootFile = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *nominalRootFile1 = new TH2D("","", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    
    TString mainPath       = "C:/Users/ajentsch/Desktop/";    
    TString mainFolder     = "RootFiles_for_Systematics/";
    TString nominalFile    = "Main.root";
    TString LSFile         = "MainLS.root";
    TString topCut1        = "Main180.root";
    TString topCut2        = "Main260.root";
    TString fBarOption1    = "fBarOption1.root";
    TString fBarOption2    = "fBarOption2.root";
    
    TString finalErrorFile = "finalSystematicErrors.root";
    
    TString inputFileName  = mainPath + mainFolder + nominalFile;
    TString inputFileName1  = mainPath + mainFolder + fBarOption1;
    TString inputFileName2  = mainPath + mainFolder + fBarOption2;
    TString inputFileName3  = mainPath + mainFolder + topCut1;
    TString inputFileName4  = mainPath + mainFolder + topCut2;
    TString inputFileName5  = mainPath + mainFolder + LSFile;
    
    TString outputFileName = mainPath + mainFolder + finalErrorFile;
    
    TFile *file = new TFile(inputFileName); //Root file with nominal histograms 
    
    TFile *file1 = new TFile(inputFileName1); //Root file with  histograms
    TFile *file2 = new TFile(inputFileName2); //Root file with histograms
    TFile *file3 = new TFile(inputFileName3); //Root file with histograms
    TFile *file4 = new TFile(inputFileName4); //Root file with  histograms
    TFile *file5 = new TFile(inputFileName5); //Root file with  histograms
    
    TString fullySubtractedCorrLabel = "FullSubtractedCorr";
    TString sideBandLabel = "_SideBand_";
    TString LSLabel = "_LS_";
    
    TString centBins[3] = {"Peripheral", "MidCentral", "Central"};
    TString ptBins[4] = {"_PtBin_0", "_PtBin_1", "_PtBin_2", ""};
    
    TString mainHistogram = fullySubtractedCorrLabel + sideBandLabel + ptBins[1] + centBins[1]; // only ptBin 1, centBin 1 for now
    TString mainHistogramLS = fullySubtractedCorrLabel + LSLabel + ptBins[1] + centBins[1];     // "                              "
    
    cout << mainHistogram << endl;
    
    //Extract and save nominal file here
    tmp = (TH2D*) file->Get(mainHistogram);
    nominalRootFile = (TH2D*) tmp->Clone();
    nominalRootFile1 = (TH2D*) tmp->Clone();
    
    
    
    secondariesErrors = (TH2D*) errorsFromSecondaries(nominalRootFile);  //systematic error 1
    bFeedDownErrors   = (TH2D*) errorsFromBFeedDown(nominalRootFile1);    //systematic error 2
    
    //close file here to move onto next one.
    
    
    //Extract, save hists, and calculate fbarErrors
    
    tmp1 = (TH2D*) file1->Get(mainHistogram);
    tmp2 = (TH2D*) file2->Get(mainHistogram);
    
    subtract->Add(tmp, tmp1, 1, -1);
    
    fBarFactorErrors = (TH2D*) errorsFromFBarFactor(tmp, tmp1, tmp2, NUM_ETA_BINS, NUM_PHI_BINS);
    
    
    
    //Close file and move on to next one
    
    
    //Extract, save hists, and calculate topCutErrors
    
    tmp3 = (TH2D*) file3->Get(mainHistogram);
    tmp4 = (TH2D*) file4->Get(mainHistogram);
    
    topologyCutErrors = (TH2D*) errorsFromTopologyCut(tmp, tmp3, tmp4, NUM_ETA_BINS, NUM_PHI_BINS);
   
    
    //Close file and move on to next one
    
    
    //Extract, save hists, and calculate LS/SB errors
    
    tmp5 = (TH2D*) file5->Get(mainHistogramLS);
    
    backgroundSubtractionErrors = (TH2D*) errorsFromBackgroundSubtraction(tmp, tmp5, NUM_ETA_BINS, NUM_PHI_BINS);
    
    
    
    
    //Close file and move on to next one
    
    
    finalErrorResults = (TH2D*) calculateFinalSystematicErrors(secondariesErrors, bFeedDownErrors, fBarFactorErrors, topologyCutErrors, backgroundSubtractionErrors, NUM_ETA_BINS, NUM_PHI_BINS);
    
    TFile *output = new TFile(outputFileName, "RECREATE"); //root file to store calculated and scaled histograms
    
    double errorProjection;
    double binError;
    
    for(int i = 1; i < NUM_PHI_BINS+1; i++){
    
        errorProjection = 0;
        binError = 0;
        for(int j = 1; j < NUM_ETA_BINS+1; j++){
        
           binError = finalErrorResults->GetBinContent(j, i);
           errorProjection = errorProjection + binError;
        
    
        }
     
        binError = binError/9.0;
        cout << "Bin " << i << " " << "Error: " << binError << endl;     
        
    }
    
    
    secondariesErrors->SetNameTitle("secondariesErrors", "secondariesErrors");
    bFeedDownErrors->SetNameTitle("bFeedDownErrors", "bFeedDownErrors");
    fBarFactorErrors->SetNameTitle("fBarFactorErrors", "fBarFactorErrors");
    topologyCutErrors->SetNameTitle("topologyCutErrors", "topologyCutErrors");
    backgroundSubtractionErrors->SetNameTitle("backgroundSubtractionErrors", "backgroundSubtractionErrors");
    
    subtract->Write();
    secondariesErrors->Write();
    bFeedDownErrors->Write();
    fBarFactorErrors->Write();
    topologyCutErrors->Write();
    backgroundSubtractionErrors->Write();
    finalErrorResults->Write();
    
    output->Close();

    file->Close();
    file1->Close();
    file2->Close();
    file3->Close();
    file4->Close();
    file5->Close();
    
    TCanvas *c = new TCanvas("canvas","canvas", 800, 800);
    
    TF1 *effWeightPions   = new TF1("effWeightPions", "[0]*TMath::Exp(-([1]/x)**[2])", 0, 1.4); 
    effWeightPions->SetParameter(0, .809);
    effWeightPions->SetParameter(1, .109);
    effWeightPions->SetParameter(2, 3.224);
    
    //effWeightPions->Draw();
    
    //cout << effWeightPions->Eval(.27, 0,0,0) << endl;
    
    TF1 *effWeightD0   = new TF1("effWeightPD0", "TMath::Exp([0]*([1]-[2]*TMath::Exp(-x/[3])))", 2, 10); 
    TH1D *drawPlot = new TH1D("","", 20, 2, 10);
    drawPlot->SetAxisRange(0,.14,"Y");
    effWeightD0->SetParameter(0, 2.30259);
    effWeightD0->SetParameter(1, -1.306);
    effWeightD0->SetParameter(2, 1.9760001);
    effWeightD0->SetParameter(3, 3.3699999);
    
    drawPlot->SetStats(0);
    drawPlot->Draw();
    effWeightD0->Draw("SAME");
    
    TString tmpStr = "C:/Users/ajentsch/Desktop/func.png";
    
    c->SaveAs(tmpStr);
    
    cout << effWeightD0->Eval(.27, 0,0,0) << endl;
    
    



}