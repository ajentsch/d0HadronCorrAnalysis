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


int pythiaD0CorrelationFitter(){

    TFile * inputFile = new TFile("5M_events_unsymmeterized_7_11_2018.root");
    TFile * outputFile = new TFile("pythia_D0_Hadron_output_5M_events_unsymmeterized_7_11_2018.root", "RECREATE");
    
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 13;
    double phiBinShift = (TMath::Pi()/12.0); 

    int NUM_PT_BINS = 3;
    
    double ETA_RANGE = 2;
    double ETA_RANGE_USER = 1.69;
    
    TString sibLabel = "sibling_D0_Hadron_Corr_";
    TString mixLabel = "mixed_D0_Hadron_Corr_";
    TString ptRangeLabel[3] = {"2_3GeV", "3_4GeV", "4_10GeV"};
    TString softPionLabel = "_NoSoftPion";
    TString softPionPeakLabel = "sibling_D0_Soft_Pion_Peak_";
    TString trackMultLabel = "track_multiplicity";
    
    
    TH2D* pythiaCorrSib[3];
    TH2D* pythiaCorrMix[3];
    
    TH2D* pythiaDStarPeakSib[3];
    TH2D* dStarSibTotal = new TH2D("total_dStar_Sib", "total_dStar_Sib", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH2D* nSigmaPlot = new TH2D("nSigmaPerBin", "nSigmaPerBin", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH1D* trackMult;
    
    TH2D* sibTotal = new TH2D("total_d0_hadron_sibling", "total_d0_hadron_sibling", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* mixTotal = new TH2D("total_d0_hadron_mix", "total_d0_hadron_mix", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* delRhoTotal = new TH2D("total_d0_hadron_delRho", "total_d0_hadron_delRho", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* ratioTotal = new TH2D("total_d0_hadron_ratio", "total_d0_hadron_ratio", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH2D* fit = new TH2D("fit", "fit", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* dipoleOffsetFit = new TH2D("dipole_offset_fit", "dipole_offset_fit", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* ASGaussFit = new TH2D("AS_Gauss_fit", "AS_Gauss_fit", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH2D* residual = new TH2D("residual", "residual", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH2D* JetPeakResidual = new TH2D("jet_peak_residual", "jet_peak_residual", NUM_ETA_BINS, -ETA_RANGE, ETA_RANGE, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH1D* fitParameters = new TH1D("fit_parameters", "fit_parameters", 7, 0, 7);
    
    fitParameters->GetXaxis()->SetBinLabel(1, "offset");
    //fitParameters->GetXaxis()->SetBinLabel(2, "dipole");
	fitParameters->GetXaxis()->SetBinLabel(2, "AS Amp");
	fitParameters->GetXaxis()->SetBinLabel(3, "AS Width");
    fitParameters->GetXaxis()->SetBinLabel(4, "jet_amplitude");
    fitParameters->GetXaxis()->SetBinLabel(5, "phi_width");
    fitParameters->GetXaxis()->SetBinLabel(6, "eta_width");
    fitParameters->GetXaxis()->SetBinLabel(7, "Gaussian_mod");
    
    
    
    //calculate n-bar here
    
    trackMult = (TH1D*) inputFile->Get(trackMultLabel);
    
    double nBar;
    double nTotalTracks = 0;
    
    for(int i = 1; i < 50; i++){
    
       nTotalTracks += i*trackMult->GetBinContent(i);
       
    }
    
    nBar = nTotalTracks/trackMult->GetEntries();
    
    cout << "Avg. tracks per event: " << nBar << endl;
    
    TString str1;
    double NSib1, NSib2, NMix1, NMix2, NSib, NMix, ratio;
    
    for(int i = 0; i < 3; i++){
    
        str1 = softPionPeakLabel + ptRangeLabel[i];
		pythiaDStarPeakSib[i] = (TH2D*) inputFile->Get(str1);
		pythiaDStarPeakSib[i]->SetDirectory(0);
        
        dStarSibTotal->Add(pythiaDStarPeakSib[i]);
       
    }    
    
    
    for(int i = 0; i < 3; i++){
    
        str1 = sibLabel + ptRangeLabel[i] + softPionLabel;
        pythiaCorrSib[i] = (TH2D*) inputFile->Get(str1);
        pythiaCorrSib[i]->SetDirectory(0);
        
        str1 = mixLabel + ptRangeLabel[i];
        pythiaCorrMix[i] = (TH2D*) inputFile->Get(str1);
        pythiaCorrMix[i]->SetDirectory(0);
     
        sibTotal->Add(pythiaCorrSib[i]);
        mixTotal->Add(pythiaCorrMix[i]);
     
    }    
    
    double eta, phi, currentBinValue;
    ofstream outputTextFilesForLanny;
    
    TString textFileForLannySib = "D0_Hadron_Pythia_2_10GeV_Sib.txt";
    outputTextFilesForLanny.open(textFileForLannySib);
        
    for(int bin1 = 1; bin1 < NUM_ETA_BINS+1; bin1++){
        for(int bin2 = 1; bin2 < NUM_PHI_BINS+1; bin2++){
        
            eta = sibTotal->GetXaxis()->GetBinCenter(bin1);
            phi = sibTotal->GetYaxis()->GetBinCenter(bin2);
            currentBinValue = sibTotal->GetBinContent(bin1, bin2);
            
            if(bin2 == 6) {phi = 0.0;} 
            
            outputTextFilesForLanny << eta << "\t" << phi << "\t" << currentBinValue << endl;
        }
    }
        
    outputTextFilesForLanny.close();
    
    
    TString textFileForLannyMix = "D0_Hadron_Pythia_2_10GeV_Mix.txt";
    outputTextFilesForLanny.open(textFileForLannyMix);
        
    for(int bin1 = 1; bin1 < NUM_ETA_BINS+1; bin1++){
        for(int bin2 = 1; bin2 < NUM_PHI_BINS+1; bin2++){
        
            eta = mixTotal->GetXaxis()->GetBinCenter(bin1);
            phi = mixTotal->GetYaxis()->GetBinCenter(bin2);
            currentBinValue = mixTotal->GetBinContent(bin1, bin2);
            
            if(bin2 == 6) {phi = 0.0;} 
            
            outputTextFilesForLanny << eta << "\t" << phi << "\t" << currentBinValue << endl;
        }
    }
        
    outputTextFilesForLanny.close();
    
   
    
    
        
    NSib = sibTotal->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    NMix = mixTotal->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    
    //NSib1 = sibTotal->Integral(1, 11, 1, NUM_PHI_BINS);
    //NMix1 = mixTotal->Integral(1, 11, 1, NUM_PHI_BINS);
    
    //NSib2 = sibTotal->Integral(17, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    //NMix2 = mixTotal->Integral(17, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    
    //NSib = NSib1 + NSib2;
    //NMix = NMix1 + NMix2;
    
    ratio = NSib/NMix;    
    
    //ratio = 1.0/mixTotal->GetBinContent(13,6);
        
    cout << ratio << endl;    
        
    delRhoTotal->Add(sibTotal, mixTotal, 1, -ratio);
    ratioTotal->Divide(delRhoTotal, mixTotal, 1, ratio);
    //ratioTotal->Divide(sibTotal, mixTotal, 1, ratio);
    //double mixRatio = 1/ratio;
    sibTotal->Write();
    
    TH1D* sibProjection = (TH1D*) sibTotal->ProjectionY();
    
    mixTotal->Scale(ratio);
    
    mixTotal->Write();
    
    double oneOverD0Yield = 1.0/94131; //10M 200 GeV -- no tune, .15 associated
    //double oneOverD0Yield = 1.0/70009; 10M 7 TeV   --  no tune, .15 associated
    //double oneOverD0Yield = 1.0/3685;

    sibProjection->Scale(oneOverD0Yield);
    
    sibProjection->Write();
    
    //ratioTotal->Scale(oneOverD0Yield);
    
    
	ratioTotal->GetXaxis()->SetTitle("#Delta#eta");
	ratioTotal->GetXaxis()->CenterTitle();
	ratioTotal->GetYaxis()->SetTitle("#Delta#phi");
	ratioTotal->GetYaxis()->CenterTitle();
	ratioTotal->GetZaxis()->SetTitle("#Delta#rho/#rho_{mix}");
    ratioTotal->GetXaxis()->SetTitleSize(.09);
	ratioTotal->GetYaxis()->SetTitleSize(.09); 
	ratioTotal->GetZaxis()->SetTitleSize(.07);
	ratioTotal->GetXaxis()->SetTitleOffset(1.05);
	ratioTotal->GetYaxis()->SetTitleOffset(1.05); 
	ratioTotal->GetZaxis()->SetTitleOffset(.6);
	ratioTotal->GetXaxis()->SetLabelSize(.05);
	ratioTotal->GetYaxis()->SetLabelSize(.05);
	ratioTotal->GetZaxis()->SetLabelSize(.04);
	ratioTotal->GetZaxis()->SetNdivisions(5,3,0, 1);
    ratioTotal->SetStats(0);
    ratioTotal->SetTitle("");
    ratioTotal->Draw("SURF1");
    c1->SaveAs("d0_hadron_corr_pythia_2_10_GeV.png");

    dipoleOffsetFunc = new TF2("dipole_offset", "[0] + [1]*cos(y)" ,-ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ASGaussOffsetFunc = new TF2("offset_ASGauss", "[0] + [1]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([2]*[2]))+exp(-0.5*((y+3.14159)*(y+3.14159))/([2]*[2])))" ,-ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

	//myfit = new TF2(str1,"[0] + [1]*cos(y) + [2]*(exp(-0.5*(y*y)/([3]*[3]))+exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4]))", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    //0 - offset, 1 - dipole, 2 - amp, 3 - phi width, 4 - eta width 
    
    //myfit = new TF2(str1,"[0] + [1]*cos(y) + [2]*(exp(-0.5*(y*y)/([3]*[3]))+exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4])) + [5]*exp(-sqrt(((y*y)/([6]*[6]))+((x*x)/([7]*[7]))))", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    //0 - offset, 1 - dipole, 2 - amp, 3 - phi width, 4 - eta width, 5- exponetnial amplitude, 6 - exponential phi width, 7 - exponential eta width
    
    //myfit = new TF2(str1,"[0] + [1]*cos(y) + [2]*exp(-sqrt(((y*y)/([3]*[3]))+((x*x)/([4]*[4]))))", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    //0 - offset, 1 - dipole, 2 - amp, 3 - phi width, 4 - eta width 
    
    //myfit = new TF2(str1,"[0] + [1]*cos(y) + [2]*(exp(-0.5*(y*y)/([3]*[3]))+exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*exp(-0.5*(x*x)/([4]*[4])) + [5]*(exp(-0.5*(y*y)/([6]*[6]))+exp(-0.5*((y-6.28)*(y-6.28))/([6]*[6])))*exp(-0.5*(x*x)/([7]*[7]))", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    //0 - offset, 1 - dipole, 2 - amp, 3 - phi width, 4 - eta width 
    
    //myfit = new TF2(str1,"[0] + [1]*cos(y) + [2]*(exp(-0.5*(y*y)/([3]*[3]))+exp(-0.5*((y-6.28)*(y-6.28))/([3]*[3])))*(([4]*[4])/((x*x)+([4]*[4])))", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
	//myfit = new TF2(str1,"[0] + [1]*cos(y) + [2]*(([3]*[3])/((y*y)+([3]*[3])))*(([4]*[4])/((x*x)+([4]*[4])))", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

	//AS Guassian fit
	//myfit = new TF2(str1,"[0] + [1]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([2]*[2]))+exp(-0.5*((y+3.14159)*(y+3.14159))/([2]*[2]))) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    //Generalized Gaussian
    myfit = new TF2(str1,"[0] + [1]*(exp(-0.5*((y-3.14159)*(y-3.14159))/([2]*[2]))+exp(-0.5*((y+3.14159)*(y+3.14159))/([2]*[2]))) + [3]*(exp(-(0.5*(y*y)/([4]*[4]))**[6])+exp(-(0.5*((y-6.28)*(y-6.28))/([4]*[4]))**[6]))*exp(-(0.5*(x*x)/([5]*[5]))**[6])", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
	
	//Best fit Lanny method 200 GeV data -- with dipole
    //myfit->SetParameters(-0.268, -.2, 1.5, .31, .31);
    
    //myfit->SetParLimits(0, -.28, -0.025);
    //myfit->SetParLimits(1, -.25, 0.0);
    //myfit->SetParLimits(2, 1.45, 1.69);
    //myfit->SetParLimits(3, 0.29, .5);
    //myfit->SetParLimits(4, 0.29, .5);
	
	//Best fit Lanny method 200 GeV data -- with dipole -- tunes from Mustafa and Xin
    //myfit->SetParameters(-0.3, -.4, 2.0, .31, .31);
    
    //myfit->SetParLimits(0, -.5, 0.0);
    //myfit->SetParLimits(1, -.5, 0.0);
    //myfit->SetParLimits(2, 1.4, 2.6);
    //myfit->SetParLimits(3, 0.29, .6);
    //myfit->SetParLimits(4, 0.29, .6);
	
	
	//Best fit Lanny method 200 GeV data -- with AS gaussian
    //myfit->SetParameters(-0.22, .286, 1.11, 1.35, .31, .31);
    
    //myfit->SetParLimits(0, -.28, -0.025);
	//myfit->SetParLimits(1, 0.0, .5);
    //myfit->SetParLimits(2, 0.2, 1.3);
	//myfit->SetParLimits(3, 1.25, 2.5);
    //myfit->SetParLimits(4, 0.29, .5);
    //myfit->SetParLimits(5, 0.29, .5);
	
	//Best fit Lanny method 200 GeV data -- with AS gaussian -- tunes from Mustafa and Xin
    //myfit->SetParameters(-0.1, .4, .8, 2.1, .31, .35);
    
    //myfit->SetParLimits(0, -.5, 0.0);
	//myfit->SetParLimits(1, 0.0, .5);
    //myfit->SetParLimits(2, 0.2, 1.3);
	//myfit->SetParLimits(3, 1.25, 2.5);
    //myfit->SetParLimits(4, 0.25, .5);
    //myfit->SetParLimits(5, 0.25, .5);
    
    
    //Best fit Lanny method 200 GeV data -- generalized Gaussian near-side and away-side -- tunes from Mustafa and Xin
    //myfit->SetParameters(-0.1, .4, .8, 2.1, .31, .35);
    
    //myfit = new TF2(str1,"[0] + [1]*(exp((-0.5*((y-3.14159)/[2])^[3]))+exp((-0.5*((y+3.14159)/[2])^[3]))) + [4]*(exp(-0.5*(y/[5])^[6])+exp(-0.5*((y-6.28)/[5])^[6]))*exp(-0.5*(x/[7])^[8])", -ETA_RANGE_USER, ETA_RANGE_USER, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    myfit->SetParLimits(0, -.5, 0.0);  //offset
	myfit->SetParLimits(1, 0.0, .4);   //AS Amp
    myfit->SetParLimits(2, 0.2, 1.2);  //AS Width
	myfit->SetParLimits(3, .5, 1.5);   //NS Amp
    myfit->SetParLimits(4, 0.2, .8);  //NS phi width
    myfit->SetParLimits(5, 0.2, .8);  //NS eta width
    myfit->SetParLimits(6, 0.3, 1.1);  //NS exponent for Gaussian
   
	
	//Best fit Lanny method 7 TeV data, .3 associated -- with AS gaussian
    //myfit->SetParameters(-0.22, .286, 1.11, 1.35, .31, .31);
    
    //myfit->SetParLimits(0, -.28, -0.025);
	//myfit->SetParLimits(1, 0.0, .3);
    //myfit->SetParLimits(2, 0.2, 1.3);
	//myfit->SetParLimits(3, 2.2, 2.5);
    //myfit->SetParLimits(4, 0.29, .5);
    //myfit->SetParLimits(5, 0.25, .5);
	
	//Best fit ALICE method 200 GeV data
    //myfit->SetParameters(0.0088, -.002, .022, .34, .34);
    
    //myfit->SetParLimits(0, 0.0088, .0088);
    //myfit->SetParLimits(1, -.0022, -0.0018);
    //myfit->SetParLimits(2, 0.02, .026);
    //myfit->SetParLimits(3, .3, .4);
    //myfit->SetParLimits(4, 0.3, .4);
	
	//Best fit ALICE method 200 GeV data -- with AS gaussian
    //myfit->SetParameters(0.25, .0095, .9, .049, .38, .38);
    
    //myfit->SetParLimits(0, 0.023, .027);
	//myfit->SetParLimits(1, 0.009, .01);
    //myfit->SetParLimits(2, 0.87, .95);
	//myfit->SetParLimits(3, .046, .055);
    //myfit->SetParLimits(4, 0.29, .45);
    //myfit->SetParLimits(5, 0.29, .45);
    
   
    
    ratioTotal->Fit(myfit, "R0");
    
    fit->Eval(myfit);
    //cout << "DO IT!!!!!!!!!!" << endl;
    residual->Add(ratioTotal, fit, 1, -1);
    
   // dipoleOffsetFunc->SetParameter(0, myfit->GetParameter(0));
    //dipoleOffsetFunc->SetParameter(1, myfit->GetParameter(1));
    
    //dipoleOffsetFit->Eval(dipoleOffsetFunc);
	
	//ASGaussOffsetFunc->SetParameter(0, myfit->GetParameter(0));
	//ASGaussOffsetFunc->SetParameter(1, myfit->GetParameter(1));
	//ASGaussOffsetFunc->SetParameter(2, myfit->GetParameter(2));
	
	//ASGaussFit->Eval(ASGaussOffsetFunc);
    
    ratioTotal->GetXaxis()->SetRangeUser(-ETA_RANGE_USER, ETA_RANGE_USER);
    fit->GetXaxis()->SetRangeUser(-ETA_RANGE_USER, ETA_RANGE_USER);
    residual->GetXaxis()->SetRangeUser(-ETA_RANGE_USER, ETA_RANGE_USER);
    residual->SetMaximum(ratioTotal->GetMaximum());
    
    double data, model, error, rawnSigma, nSigma;
    
    for(int etabin = 1; etabin < NUM_ETA_BINS+1; etabin++){
        for(int phibin = 1; phibin < NUM_PHI_BINS+1; phibin++){
        
            data = ratioTotal->GetBinContent(etabin, phibin);
            model = fit->GetBinContent(etabin, phibin);
            error = ratioTotal->GetBinError(etabin, phibin);
            
            rawnSigma = (data-model)/error;
            nSigma = TMath::Abs(rawnSigma);
            
            nSigmaPlot->SetBinContent(etabin, phibin, nSigma);
        }
    }
    
    //JetPeakResidual->Add(ratioTotal, dipoleOffsetFit, 1, -1);
	//JetPeakResidual->Add(ratioTotal, ASGaussFit, 1, -1);
    
    cout << "Jet yield (if using ALICE method): " << JetPeakResidual->Integral(9,17, 3, 9) << endl;
	
	double pythiaNBar = 8.49823;
    //double pythiaNBar = 5;
    //double pythiaIntEta = .520653;
	
	double jetAmp = myfit->GetParameter(3);
	double phiWidth = myfit->GetParameter(4);
	double etaWidth = myfit->GetParameter(5);
	double jetAmpError = myfit->GetParError(3);
	double phiWidthError = myfit->GetParError(4);
	double etaWidthError = myfit->GetParError(5);
    
    double jetYieldLanny = ((pythiaNBar*jetAmp*phiWidth*etaWidth)/2.0);
	
	double jetYieldError = (pythiaNBar/2.0)*getJetYieldError(jetAmp, etaWidth , phiWidth, jetAmpError, etaWidthError , phiWidthError);

   
	
	cout << "Jet yield (if using Lanny method): " << jetYieldLanny << " +/- " << jetYieldError << endl;
    
    fitParameters->SetBinContent(1, myfit->GetParameter(0));
    fitParameters->SetBinContent(2, myfit->GetParameter(1));
    fitParameters->SetBinContent(3, myfit->GetParameter(2));
    fitParameters->SetBinContent(4, myfit->GetParameter(3));
    fitParameters->SetBinContent(5, myfit->GetParameter(4));
	fitParameters->SetBinContent(6, myfit->GetParameter(5));
    fitParameters->SetBinContent(7, myfit->GetParameter(6));

	
    
    fitParameters->SetBinError(1, myfit->GetParError(0));
    fitParameters->SetBinError(2, myfit->GetParError(1));
    fitParameters->SetBinError(3, myfit->GetParError(2));
    fitParameters->SetBinError(4, myfit->GetParError(3));
    fitParameters->SetBinError(5, myfit->GetParError(4));
	fitParameters->SetBinError(6, myfit->GetParError(5));
    fitParameters->SetBinError(7, myfit->GetParError(6));

    
    fitParameters->Write();
    ratioTotal->Write();
    fit->Write();
    residual->Write();
    
    dStarSibTotal->Write();
    
    TCanvas * finalCan = new TCanvas("canvas", "canvas", 4400, 1000);
    finalCan->Divide(4,1);
    
    finalCan->cd(1);
    ratioTotal->Draw("SURF1");
    finalCan->cd(2);
    fit->SetStats(0);
    fit->Draw("SURF1");
    finalCan->cd(3);
    residual->SetStats(0);
    residual->Draw("SURF1");
    finalCan->cd(4);
    fitParameters->SetMarkerStyle(20);
    fitParameters->SetMarkerColor(1);
    fitParameters->SetMarkerSize(2);
    fitParameters->SetLineColor(1);
    //fitParameters->Draw("P");
    nSigmaPlot->GetXaxis()->SetRangeUser(-ETA_RANGE_USER, ETA_RANGE_USER);
    nSigmaPlot->SetStats(0);
    nSigmaPlot->Draw("COLZ");
    
    //finalCan->cd(5);
    //dipoleOffsetFit->Draw("SURF1");
    //ASGaussFit->Draw("SURF1");
	//finalCan->cd(6);
    //JetPeakResidual->Draw("SURF1");
    //finalCan->cd(7);
    //TH1D* jetPeakProj = (TH1D*)JetPeakResidual->ProjectionY();
    //jetPeakProj->Draw();
    //jetPeakProj->Write();
    finalCan->SaveAs("final_pythia_results.png");
    
    TCanvas * dStarCan = new TCanvas("canvas", "canvas", 1800, 1300);
    dStarSibTotal->GetXaxis()->SetTitle("#Delta#eta");
	dStarSibTotal->GetXaxis()->CenterTitle();
	dStarSibTotal->GetYaxis()->SetTitle("#Delta#phi");
	dStarSibTotal->GetYaxis()->CenterTitle();
	//dStarSibTotal->GetZaxis()->SetTitle("#Delta#rho/#rho_{ref}");
    dStarSibTotal->GetXaxis()->SetTitleSize(.06);
	dStarSibTotal->GetYaxis()->SetTitleSize(.06); 
	//dStarSibTotal->GetZaxis()->SetTitleSize(.07);
	dStarSibTotal->GetXaxis()->SetTitleOffset(.8);
	dStarSibTotal->GetYaxis()->SetTitleOffset(.7); 
	//dStarSibTotal->GetZaxis()->SetTitleOffset(.6);
	dStarSibTotal->GetXaxis()->SetLabelSize(.05);
	dStarSibTotal->GetYaxis()->SetLabelSize(.05);
    
    dStarSibTotal->GetXaxis()->SetRangeUser(-.8, .8);
    dStarSibTotal->GetYaxis()->SetRangeUser(-.8, .8);
	//dStarSibTotal->GetZaxis()->SetLabelSize(.04);
	//dStarSibTotal->GetZaxis()->SetNdivisions(5,3,0, 1);
	dStarSibTotal->SetTitle("");
    dStarSibTotal->Draw("COLZ2");
    dStarSibTotal->SetStats(0);
	dStarCan->cd(1)->SetLeftMargin(0.1);
	dStarCan->cd(1)->SetRightMargin(0.13);
	dStarCan->cd(1)->SetBottomMargin(0.11);
    dStarCan->SaveAs("dStar_peak_corr_plot.png");
    
    
    fitParameters->Draw("P");
    dStarCan->SaveAs("fit_parameters.png");
    
    inputFile->Close();
    outputFile->Close();
    
    
}
        
        
        
        
double getJetYieldError(double A, double sigEta, double sigPhi, double eA, double eSigEta, double eSigPhi){

    cout << endl <<"error on A, sigPhi, sigEta - " << eA << " , " << eSigPhi << " , " << eSigEta << endl;

    double arg = ((eA/A)*(eA/A)) + ((eSigEta/sigEta)*(eSigEta/sigEta)) + ((eSigPhi/sigPhi)*(eSigPhi/sigPhi));
    double product = TMath::Abs((A*sigEta*sigPhi));
    
    //cout << "prduct value: " << product << endl;
    //cout << "arg value: " << arg << endl;
    //cout << "sqrt of arg: " << TMath::Sqrt(arg) << endl;
    
    return TMath::Abs((A*sigEta*sigPhi))*TMath::Sqrt(arg);
    
}     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    