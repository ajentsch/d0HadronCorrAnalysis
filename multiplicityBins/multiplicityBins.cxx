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


double correctedFromRaw(double raw, double alpha, double beta){

    return ((raw + alpha*raw*raw)/(beta));
    
}    

double getCorrectedPowerLawValue(double nChargedMinusThreeFourth, double NchCorr, double alpha, double beta){

    double final = 0;
    double denominator = TMath::Sqrt(1+(4*alpha*beta*NchCorr));
    double nChargedCorrThreeFourth = TMath::Power(NchCorr, .75);
    
    final = (beta*nChargedMinusThreeFourth*nChargedCorrThreeFourth)/(denominator);
    
    return final;
}    

int multiplicityBins(){

    TString inputFileName  = "best_nominal_new_eff_correction.root";
    TString outputFileName = "multiplicityOutput.root"; 
    TString run4FileName = "Run4PowerLaw.txt";

    TFile *file = new TFile(inputFileName); //Root file with raw histograms to use for calculations
    
    TFile *output = new TFile(outputFileName, "RECREATE"); //root file to store calculated and scaled histograms
    
    TString multiplicityDistLabel = "number_of_tracks_per_event";
    
    TH1I* h = new TH1I("number_of_tracks_per_event", "number_of_tracks_per_event", 1500, 0, 1499);
    TH1I* multiplicityDist = new TH1I(, , 1500, 0, 1499);
    TH1D* powerLawPlot = new TH1D("Power Law", "Power Law", 100, 0, 7);
    TH1D* normPowerLawPlot = new TH1D("Power Law -- Normalized", "Power Law -- Normalized", 100, 0, 7);
    TH1D* normCorrectedPowerLawPlot = new TH1D("Power Law -- Normalized & Corrected", "Power Law -- Normalized & Corrected", 100, 0, 7);
    TH1D* normRun4PowerLawPlot = new TH1D("run 4 Power Law -- Normalized", "run 4 Power Law -- Normalized", 100, 0, 7);
    
    //BEGIN -- Read in text file to generate run4 power law histogram 
    //**********************
    ifstream run4Power(run4FileName);
    string line;
    
    while(getline(run4Power, line)){
    
        stringstream linestream(line);
        string data;
        double val1;
        double val2;
        int bin;
        
        linestream >> val1 >> val2;
        bin = normRun4PowerLawPlot->FindBin(val1);
        normRun4PowerLawPlot->SetBinContent(bin, val2);
    
    }
    
    run4Power.close();
    //**********************
    //END
    
    double numEntriesPowerLaw = 0;
    
    double numEntriesPowerLawRun4Range = 0;
    double numEntriesPowerLawRange = 0;
    
    powerLawPlot->GetXaxis()->SetTitle("N_{ch}^{1/4}");
    powerLawPlot->GetXaxis()->SetTitleOffset(1.3);
    
    powerLawPlot->GetYaxis()->SetTitle("#frac{dN}{N_{ch}^{1/4}}");
    powerLawPlot->GetYaxis()->SetTitleOffset(1.2);
    
    normPowerLawPlot->GetXaxis()->SetTitle("N_{ch}^{1/4}");
    normPowerLawPlot->GetXaxis()->SetTitleOffset(1.3);
    
    normPowerLawPlot->GetYaxis()->SetTitle("#frac{dN}{N_{ch}^{1/4}}");
    normPowerLawPlot->GetYaxis()->SetTitleOffset(1.3);
    
    normCorrectedPowerLawPlot->GetXaxis()->SetTitle("N_{ch}^{1/4}");
    normCorrectedPowerLawPlot->GetXaxis()->SetTitleOffset(1.3);
    
    normCorrectedPowerLawPlot->GetYaxis()->SetTitle("#frac{dN}{N_{ch}^{1/4}}");
    normCorrectedPowerLawPlot->GetYaxis()->SetTitleOffset(1.3);
    
    normRun4PowerLawPlot->GetXaxis()->SetTitle("N_{ch}^{1/4}");
    normRun4PowerLawPlot->GetXaxis()->SetTitleOffset(1.3);
    
    normRun4PowerLawPlot->GetYaxis()->SetTitle("#frac{dN}{N_{ch}^{1/4}}");
    normRun4PowerLawPlot->GetYaxis()->SetTitleOffset(1.3);
    
    TH1D* h1 = new TH1D("unlikeSign", "unlikeSign", 50, 1.6, 2.1);
    TH1D* invMass = new TH1D(, , 50, 1.6, 2.1);
    
    
    h = (TH1I*) file->Get(multiplicityDistLabel);
    multiplicityDist = (TH1I*)h->Clone(multiplicityDistLabel);
    
    multiplicityDist->SetTitle("Track Multiplicity per Event");
    
    multiplicityDist->GetXaxis()->SetTitle("Log(N_{ch})");
    multiplicityDist->GetXaxis()->SetTitleOffset(1.1);
    
    multiplicityDist->GetYaxis()->SetTitle("Log(#frac{dN}{N_{ch}})");
    multiplicityDist->GetYaxis()->SetTitleOffset(1.1);
    
    h1 = (TH1D*) file->Get("unlikeSign");
    invMass = (TH1D*)h1->Clone("unlikeSign");
    
    int numBins = 15;  //  0   1   2   3   4    5    6    7    8    9    10   11   12   13   14   15   16   17
    int multBins[16]  = {0, 16, 31, 51, 79, 116, 164, 224, 295, 335, 376, 420, 470, 550, 650, 775};  
                       //  0   1   2   3   4    5    6    7      8     9          10
    int integralCheck = 0;                                    
    
    for( int i = 0; i < numBins; i++){
        
        cout << "Bin: " << i << endl;
        cout << "Bin Range: " << multBins[i] << " - " << multBins[i+1] << "." << endl;
        double integralValue = multiplicityDist->Integral(multBins[i]+1, multBins[i+1]);
        cout << "Integral of this range: " << integralValue << endl << endl;
        integralCheck = integralCheck + integralValue;
        
    }
    
    cout << "integral Check: " << integralCheck << endl << endl;
    
    double integralValue = multiplicityDist->Integral(1, 775);
    double binSize = integralValue/numBins;
    
    cout << "Integral of multiplicity distribution: " << integralValue << endl;
    cout << "Equal space for " << numBins << " bins: " << binSize << endl;
    
    
    //Generate power law plot here
    double binContent = 0;
    double nChargedOneFourth = 0;
    double nChargedThreeFourth = 0;
    double nChargedMinusThreeFourth = 0;
    double powerLawValue = 0;
    double powerLawBin = 0;
    double alpha = 0.0000014;
    double beta = .97;
    double dNevOverdNChRawOneFourth = 0;
    double NchCorr = 0;
    double powerLawBinCorr = 0;
    double nChargedCorrOneFourth = 0;
    double numEntriesPowerLawCorr = 0;
    
    //double nChargedCorrThreeFourth = 0;
    double correctedPowerLawValue = 0;
    
    for(int i = 1; i < 1000; i++){//i is N_Ch,Raw
    
        binContent = multiplicityDist->GetBinContent(i);
        nChargedOneFourth   = TMath::Power(i, .25);        //N_Ch,Raw^1/4
        nChargedThreeFourth = TMath::Power(i, .75);        //N_Ch,Raw^3/4
        nChargedMinusThreeFourth = TMath::Power(i, -.75); //N_Ch,Raw^-3/4
        
        
        
        //This section generates the simple Jacobian transform from dN_ev/dN_Ch,Raw -> dN_ev/dN_Ch,Raw^1/4
        powerLawValue = 4*nChargedThreeFourth*binContent;       
        powerLawBin = powerLawPlot->FindBin(nChargedOneFourth);
        powerLawPlot->SetBinContent(powerLawBin, powerLawValue);
        normPowerLawPlot->SetBinContent(powerLawBin, powerLawValue); //just a clone of the above plot to later have a normalized and raw version
        
        //This section generates the CORRECTED Jacobian transform from dN_ev/dN_Ch,Raw^1/4 -> dN_ev/dN_Ch,Corr^1/4
        dNevOverdNChRawOneFourth = powerLawValue;
        
        NchCorr = correctedFromRaw(i, alpha, beta);
        nChargedCorrOneFourth = TMath::Power(NchCorr, .25);
        
        correctedPowerLawValue = getCorrectedPowerLawValue(nChargedMinusThreeFourth, NchCorr, alpha, beta)*dNevOverdNChRawOneFourth;
        powerLawBinCorr = powerLawPlot->FindBin(nChargedCorrOneFourth);
        
        normCorrectedPowerLawPlot->SetBinContent(powerLawBinCorr, correctedPowerLawValue);
        
        
    }

    for(int i = 1; i < 16; i++){
    
        cout << "Bin " << i << ":" << TMath::Power(multBins[i], .25) << endl;
    }
    
    numEntriesPowerLaw = powerLawPlot->Integral(1,100);
    numEntriesPowerLawCorr = normCorrectedPowerLawPlot->Integral(1,100);
    numEntriesPowerLawRun4Range = normRun4PowerLawPlot->Integral(62,71);
    numEntriesPowerLawRange = powerLawPlot->Integral(62,71);
   
    //normPowerLawPlot->Scale(numEntriesPowerLawRun4Range/numEntriesPowerLawRange);
    normPowerLawPlot->Scale(1/numEntriesPowerLaw);
    normPowerLawPlot->SetLineColor(kRed);
    normCorrectedPowerLawPlot->Scale(1/numEntriesPowerLawCorr);
    
    normCorrectedPowerLawPlot->SetLineColor(kRed);
    
    TCanvas *c = new TCanvas("c2", "Histograms", 1100, 850);
    
    //Draw and write plots to file
    powerLawPlot->Draw();
    powerLawPlot->Write();
    c->SaveAs("powerLaw.png");    
    normPowerLawPlot->Draw();
    normPowerLawPlot->Write();
    c->SaveAs("powerLaw_normalized.png"); 
    normRun4PowerLawPlot->Draw("C");
    normRun4PowerLawPlot->Write();
    normCorrectedPowerLawPlot->Draw("C");
    normCorrectedPowerLawPlot->Write();
    c->SaveAs("powerLaw_normalized_corrected.png");
    multiplicityDist->Draw();
    multiplicityDist->Write();    
    
    normCorrectedPowerLawPlot->Draw("C");
    normRun4PowerLawPlot->Draw("C SAME");
    c->SaveAs("mult.png");
    
    normPowerLawPlot->Draw("C");
    normRun4PowerLawPlot->Draw("C SAME");
    c->SaveAs("mult2.png");
    
    int binOne = 0;
    int binTwo = 0;
    
    binOne = invMass->FindBin(1.62);
    binTwo = invMass->FindBin(1.7);
    
    cout << "Integral of sideBandLeft : " << invMass->Integral(binOne, binTwo) << endl;
    binOne = invMass->FindBin(1.82);
    binTwo = invMass->FindBin(1.9);
    cout << "Integral of signalRegion : " << invMass->Integral(binOne, binTwo) << endl;
    binOne = invMass->FindBin(2.0);
    binTwo = invMass->FindBin(2.1);
    cout << "Integral of sideBandRight : " << invMass->Integral(binOne, binTwo) << endl;
    
    ofstream powerLawPlotText;
    powerLawPlotText.open("powerLawPlotText.txt");
    
    for(int i = 1; i < 101; i++){
    
        powerLawPlotText << powerLawPlot->GetBinCenter(i) <<"\t" << powerLawPlot->GetBinContent(i) << endl;
    
    }
    
    powerLawPlotText.close();
    
    //////////////////////////////////
    
    ofstream normPowerLawPlotText;
    normPowerLawPlotText.open("normPowerLawPlotText.txt");
    
    for(int i = 1; i < 101; i++){
    
        normPowerLawPlotText << normPowerLawPlot->GetBinCenter(i) <<"\t" << normPowerLawPlot->GetBinContent(i) << endl;
    
    }
    
    normPowerLawPlotText.close();
    
    ///////////////////////////////////
    
    ofstream normCorrectedPowerLawPlotText;
    normCorrectedPowerLawPlotText.open("normCorrectedPowerLawPlotText.txt");
    
    for(int i = 1; i < 101; i++){
    
        normCorrectedPowerLawPlotText << normCorrectedPowerLawPlot->GetBinCenter(i) <<"\t" << normCorrectedPowerLawPlot->GetBinContent(i) << endl;
    
       
    
    }
    
    normCorrectedPowerLawPlotText.close();
    
    
    
    file->Close();
    output->Close();

}



