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

#include "C:/Users/ajentsch/desktop/corrHistogramMaker.h"
//void formatCorrHist(TH2D* hist);
//void formatCorrHist(TH2D* hist, TString title);







Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if ( (reject && (x[0] > 1.81 && x[0] < 1.91)) || (reject && (x[0] > 1.6 && x[0] < 1.7)) ){     //|| (reject && (x[0] > 2.0 && x[0] < 2.1))
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

double calcScaleFactorforFit(double sigma)
{
    double preFactor = TMath::Sqrt(TMath::PiOver2());
    
    return preFactor*sigma*.5;
    
}
        

int D0corrFittingAndPlotMaker(){

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------FITTING AND PARAMETER EXTRACTION DONE HERE---------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    TCanvas *c = new TCanvas("c3", "Histograms", 1100, 850);
 
    TString fitFunctionLabel = "parameterFitFunction";
 
    //TF2 *myfit   = new TF2("parameterFitFunction","[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
    
    TF2 *myfit[4][3];   
    
    TF1 *myfit1D = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x)" , -1.57, 4.71);
    
    TF1 *fitProjection = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x) + [3]*[4]*(exp(-0.5*(x*x)/([5]*[5]))+exp(-0.5*((x-6.28)*(x-6.28))/([5]*[5])))" , -1.57, 4.71);
    
    TF1 *quadrupole1    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole2    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole3    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    
    TF1 *pileupLine1    = new TF1("PileupLine1", "[0] + [1]*x", -2.0, -1.0);
    TF1 *pileupLine2    = new TF1("PileupLine2", "[0] + [1]*x", -1.0, 0);
    TF1 *pileupLine3    = new TF1("PileupLine3", "[0] + [1]*x", 0, 1.0);
    TF1 *pileupLine4    = new TF1("PileupLine4", "[0] + [1]*x", 1.0, 2.0);
    
    TH1D *testPileupFunc = new TH1D("test", "test", 100, -2, 2);
    
    //testPileupFunc->SetMaximum(1);
    //testPileupFunc->SetMinimum(-1);
    
    pileupLine1->FixParameter(0, -2);
    pileupLine1->FixParameter(1, -1);
    pileupLine2->FixParameter(0, 1);
    pileupLine2->FixParameter(1, 2);
    pileupLine3->FixParameter(0, 1);
    pileupLine3->FixParameter(1, -2);
    pileupLine4->FixParameter(0, -2);
    pileupLine4->FixParameter(1, 1);
    
    str1 = path + outputFolders[5] + outputFolders[8] + fileType;
    testPileupFunc->Draw();
    pileupLine1->Draw("SAME");
    pileupLine2->Draw("SAME");
    pileupLine3->Draw("SAME");
    pileupLine4->Draw("SAME");
    
    c->SaveAs(str1); 
 
    TH2D* fullySubtractedCorrCustomCentralitySideBandFit[4][3];
    TH2D* fullySubtractedCorrCustomCentralitySideBandRes[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[4][3];
    
    TH2D* fullySubtractedCorrCustomCentralitySideBandPileupCorrected[4][3];
    //TH2D* fullySubtractedCorrCustomCentralitySideBandRes[4][3];
    //TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[4][3];
    //TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[4][3];
    
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
 
    for(int k = 0; k < 4; k++){
        for(int i = 0; i < 3; i++){ 
 
            str1 = fitFunctionLabel + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
            myfit[k][i] = new TF2(str1,"[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
           
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            str1 = phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
 
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + pileupLabel;
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
 
            fitParameterFile2D << "_________________________________________________________________________________________________________" << endl;
            fitParameterFile1D << "_________________________________________________________________________________________________________" << endl;
 
           
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone();

              
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);

            myfit1D->SetParLimits(0, -.5, .01); //offset A0
            myfit1D->SetParLimits(1, -.1, 0.0); //v1     AD
            myfit1D->SetParLimits(2, 0.0, .1); //v2     AQ
            
            myfit[k][i]->SetParLimits(0, -.1, .01);    //offset A0
            myfit[k][i]->SetParLimits(1, -.05,  0.0);    //v1     AD
            myfit[k][i]->SetParLimits(2, 0.0,  .1);    //v2     AQ
            myfit[k][i]->SetParLimits(3,  0.0,  .19);    //Jetamp A1
            myfit[k][i]->SetParLimits(4, 0.1, 2);     //sigPhi
            myfit[k][i]->SetParLimits(5, 0.1, 2);    //sigEta
            
            //////////////////////2d fitting///////////////////////////////////////
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Fit(myfit[k][i], "R0E");
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
            
            fullySubtractedCorrCustomCentralitySideBandRes[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Add(fullySubtractedCorrCustomCentralitySideBandFit[k][i], -1);
            
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Draw("SURF1");
            c->SaveAs(str1);
 
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Write();
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Draw("SURF1");
            c->SaveAs(str1); 
            
            /////////////////////////////1d fitting///////////////////////////////////////////
    
            
    
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Fit(myfit1D, "qR0E");
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
            
            fitParameterFile2D << "V2 = " << v2_2D << " +/- " << v2_2D_Error << endl;
            fitParameterFile1D << "V2 = " << v2_1D << " +/- " << v2_1D_Error << endl;
            
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
            
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel + fileType;
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Draw();
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw("EX0");
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetMarkerStyle(20);
            
            fitProjection->Draw("SAME");
            quadrupole1->Draw("SAME");
            quadrupole2->Draw("SAME");
            //quadrupole3->Draw("SAME");
            myfit1D->Draw("SAME");
            c->SaveAs(str1);
 
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Draw();
            c->SaveAs(str1); 
 
 
        }
    }    
 
    fitParameterFile2D.close();
    fitParameterFile1D.close();
 
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------Histogram formatting and output for talks----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    TCanvas *talkCanvas = new TCanvas("c3", "Histograms", 3700, 1700);
    talkCanvas->Divide(3,2);
    
    int padNumber = 1;
    TString talkOutput = "Histograms_for_talk";
    TString ptRanges[2] = {" 1 < p_{t} < 4 GeV ", " 4 < p_{t} < 20 GeV "};
    TString centRanges[3] = {" 50-100% ", " 20-50% ", " 0-20% "};
    
    double xValueforPave[3] = {1000, 2100, 3200};
    double yValueforPave[2] = {300, 1150};
    
    TString firstLine;
    TString secondLine;
    
    //PaveText *textBox = new TPaveText(
    
    for(int k = 1; k < 3; k++){
            for(int i = 0; i < 3; i++){
 
                
                formatCorrHist(fullySubtractedCorrCustomCentralitySideBand[k][i]);
                
                talkCanvas->cd(padNumber);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle();
                //fullySubtractedCorrCustomCentralitySideBand[k][i]->GetYaxis()->SetTitleOffset(
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
                
                TPaveText *textBox = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
                
                firstLine = ptRanges[k-1];
                secondLine = centRanges[i];
                
                textBox->AddText(firstLine.Data());
                textBox->AddText(secondLine.Data());
                textBox->Draw();
 
 
                padNumber++;
        }    
    }
 
    str1 = path + talkOutput + fileType;
    talkCanvas->SaveAs(str1);

    //////////////////////////////////////////////////////////////////////////////
	
	TCanvas *talkCanvas2 = new TCanvas("c3", "Histograms", 1200, 1700);
    talkCanvas2->Divide(1,2);
    
    int padNumber = 1;
    TString talkOutput2 = "more_histograms_for_talk";
   
                
				firstLine = ptRanges[0];
                secondLine = centRanges[1];
	             
                talkCanvas2->cd(1);
                oneDUSHistos[1][1]->SetStats(false);
				oneDSubtractedInvMassHistosFunction[1][1]->SetStats(false);
				
				oneDUSHistos[1][1]->SetTitle();
				oneDSubtractedInvMassHistosFunction[1][1]->SetTitle();
                oneDSubtractedInvMassHistosFunction[1][1]->GetXaxis()->SetTitleSize(.04);
                oneDSubtractedInvMassHistosFunction[1][1]->GetYaxis()->SetTitleSize(.04);
                oneDSubtractedInvMassHistosFunction[1][1]->GetXaxis()->SetTitleFont(62);
                oneDSubtractedInvMassHistosFunction[1][1]->GetYaxis()->SetTitleFont(62);
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelSize(.04, "X");
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelSize(.04, "Y");
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelFont(62, "X");
                oneDSubtractedInvMassHistosFunction[1][1]->SetLabelFont(62, "Y");
                
                oneDUSHistos[1][1]->GetXaxis()->SetTitleSize(.04);
                oneDUSHistos[1][1]->GetYaxis()->SetTitleSize(.04);
                oneDUSHistos[1][1]->GetXaxis()->SetTitleFont(62);
                oneDUSHistos[1][1]->GetYaxis()->SetTitleFont(62);
                oneDUSHistos[1][1]->SetLabelSize(.04, "X");
                oneDUSHistos[1][1]->SetLabelSize(.04, "Y");
                oneDUSHistos[1][1]->SetLabelFont(62, "X");
                oneDUSHistos[1][1]->SetLabelFont(62, "Y");
                
                oneDUSHistos[1][1]->Draw();
                
                TPaveText *textBox2 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
				textBox2->AddText(firstLine.Data());
                textBox2->AddText(secondLine.Data());
                textBox2->Draw();
                
                talkCanvas2->cd(2);
                
				oneDSubtractedInvMassHistosFunction[1][1]->Draw();
				
                TPaveText *textBox3 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
				textBox3->AddText(firstLine.Data());
                textBox3->AddText(secondLine.Data());
                textBox3->Draw();
 
 
				str1 = path + talkOutput2 + fileType;
				talkCanvas2->SaveAs(str1);

				//////////////////////////////////////////////////////////////////////////////
                
				TString talkOutput3 = "more_more_output_histograms";
				
				TCanvas *talkCanvas3 = new TCanvas("c3", "Histograms", 3700, 850);
                talkCanvas3->Divide(3,1);

                

                formatCorrHist(fullySubtractedCorrCustomCentralitySideBandFit[2][1]);
				formatCorrHist(fullySubtractedCorrCustomCentralitySideBandRes[2][1]);
                
                talkCanvas3->cd(1);
                fullySubtractedCorrCustomCentralitySideBand[2][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBand[2][1]->SetTitle();
                fullySubtractedCorrCustomCentralitySideBand[2][1]->Draw("SURF1");
                
                TPaveText *textBox4 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
                
                firstLine = ptRanges[1];
                secondLine = centRanges[1];
                
                textBox4->AddText(firstLine.Data());
                textBox4->AddText(secondLine.Data());
                textBox4->Draw();

                talkCanvas3->cd(2);
                fullySubtractedCorrCustomCentralitySideBandFit[2][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBandFit[2][1]->SetTitle();
                fullySubtractedCorrCustomCentralitySideBandFit[2][1]->Draw("SURF1");
                
                TPaveText *textBox5 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
                
                firstLine = ptRanges[1];
                secondLine = centRanges[1];
                
                textBox5->AddText(firstLine.Data());
                textBox5->AddText(secondLine.Data());
                textBox5->Draw();

                talkCanvas3->cd(3);
                fullySubtractedCorrCustomCentralitySideBandRes[2][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBandRes[2][1]->SetTitle();
				fullySubtractedCorrCustomCentralitySideBandRes[2][1]->SetAxisRange(-.02,.1,"Z");
				//fullySubtractedCorrCustomCentralitySideBandRes[1][1]->SetMinimum(-.02)
                fullySubtractedCorrCustomCentralitySideBandRes[2][1]->Draw("SURF1");
                
                TPaveText *textBox6 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
                
                firstLine = ptRanges[1];
                secondLine = centRanges[1];
                
                textBox6->AddText(firstLine.Data());
                textBox6->AddText(secondLine.Data());
                textBox6->Draw();

                str1 = path + talkOutput3 + fileType;
				talkCanvas3->SaveAs(str1);

                //////////////////////////////////////////////////////--------------------------------------------------------
				
				TString talkOutput4 = "more_more_more_output_histograms";
				
				TCanvas *talkCanvas4 = new TCanvas("c3", "Histograms", 1300, 850);
                
                TF1 *offset = new TF1("offset" , "[0]", -1.57, 4.71);
                TF1 *offsetPlusDipole = new TF1("offset+dipole", "[0]+[1]*cos(x)", -1.57, 4.71);
                TF1 *offsetPlusDipolePlusQuad = new TF1("offset+dipole+quad", "[0]+[1]*cos(x)+2*[2]*cos(2*x)", -1.57, 4.71);
                TF1 *offsetPlusQuad = new TF1("offsetquad", "[0]+2*[1]*cos(2*x)", -1.57, 4.71);
                TF1 *dipole = new TF1("dipole", "[0]*cos(x)", -1.57, 4.71);
                //TF1 *offsetPlusDipole = new TF1("offset+dipole", "[0]+[1]*cos(x)", -1.57, 4.71);
                
                phiProjSyst = new TH1D("", "", NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                
                
                //USE QUADRUPOLE 3
                //TF1 *gaussian = new TF1("
				
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetStats(false);
				//fullySubtractedCorrCustomCentralitySideBandPhiProjFit[1][1]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetXaxis()->SetTitle("#Delta#phi");
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetXaxis()->CenterTitle();
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetXaxis()->SetTitleOffset(1.3);
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetYaxis()->SetTitle("#frac{#Delta#rho}{#rho_{ref}}");
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetYaxis()->SetTitleOffset(1.5);
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetTitle();
				
                double phiValue[12];
                double corrValue[12];
                double ePhi[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
                double eCorrSyst[12] = {0.000694801, 0.000293818, 0.000769675, 0.000769675, 0.000293818, 0.000694801,
                                        0.000758256, 0.000476921, 0.000569372, 0.000569372, 0.000476921, 0.000758256};
                
                
                
                
                
                for( int i = 1; i < NUM_PHI_BINS+1; i++){
                
                    phiValue[i-1] = fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetBinCenter(i);
                    corrValue[i-1] = fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->GetBinContent(i);
                }
                
                
                
                TGraphErrors *PhiDataSyst = new TGraphErrors(12, phiValue, corrValue, ePhi, eCorrSyst);
  
                //PhiDataSyst->SetMarkerColor(16);
                //PhiDataSyst->SetMarkerStyle(21);
               // PhiDataSyst->SetMarkerSize(1.1);
                PhiDataSyst->SetLineColor(1);
                //PhiDataSyst->SetLineWidth(3);
                
				fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->Draw("EX0");
                fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->SetMarkerStyle(20);
                PhiDataSyst->Draw("[] SAME");
            
			    fullySubtractedCorrCustomCentralitySideBandPhiProj[1][1]->Fit(myfit1D, "qR0E");
				fullySubtractedCorrCustomCentralitySideBand[1][1]->Fit(myfit[1][1], "qR0E");
			
			    GaussIntegralParameter = calcScaleFactorforFit(myfit[k][i]->GetParameter(5));
			
			    fitProjection->SetParameter(0, myfit[1][1]->GetParameter(0));
				fitProjection->SetParameter(1, myfit[1][1]->GetParameter(1));
				fitProjection->SetParameter(2, myfit[1][1]->GetParameter(2));
				fitProjection->SetParameter(3, myfit[1][1]->GetParameter(3));
				fitProjection->SetParameter(5, myfit[1][1]->GetParameter(4));
				fitProjection->SetParameter(4, GaussIntegralParameter);
			
                offset->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusDipole->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusDipole->SetParameter(1, myfit[1][1]->GetParameter(1));
                offsetPlusDipolePlusQuad->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusDipolePlusQuad->SetParameter(1, myfit[1][1]->GetParameter(1));
                offsetPlusDipolePlusQuad->SetParameter(2, myfit[1][1]->GetParameter(2));
                offsetPlusQuad->SetParameter(0, myfit[1][1]->GetParameter(0));
                offsetPlusQuad->SetParameter(1, myfit[1][1]->GetParameter(2));
                dipole->SetParameter(0, myfit[1][1]->GetParameter(1));
                quadrupole3->SetParameter(0, myfit[1][1]->GetParameter(2));
                //quadrupole3->SetParameter(0, myfit[k][i]->GetParameter(2));
            
                dipole->SetLineColor(4); //blue
                //offsetPlusDipole->SetLineColor(3); //green
                //offsetPlusDipolePlusQuad->SetLineColor(2); //red
                quadrupole3->SetLineColor(2);
                quadrupole3->SetLineWidth(2);
                fitProjection->SetLineColor(7); //cyan
            
                dipole->Draw("SAME");
                //offsetPlusDipole->Draw("SAME");
                quadrupole3->Draw("SAME");
                fitProjection->Draw("SAME"); 

				TPaveText *textBox7 = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
                
                firstLine = ptRanges[0];
                secondLine = centRanges[1];
                
                textBox7->AddText(firstLine.Data());
                textBox7->AddText(secondLine.Data());
                textBox7->Draw();

                str1 = path + talkOutput4 + fileType;
				talkCanvas4->SaveAs(str1);
				

    //------------------------------------------------------v2 plots-----------------------------------------------

    TCanvas *v2Canvas = new TCanvas("c3", "Histograms", 1100, 850);
    
    TString v2 = "v2";
    
    TH1D *v2Histogram = new TH1D("V_{2} 20-50%","V_{2} 20-50%", 10, 0, 9);

    v2Histogram->SetBinContent(1, v2Array[0][1]);  
    v2Histogram->SetBinContent(3, v2Array[1][1]);
    v2Histogram->SetBinContent(5, v2Array[2][1]);

    v2Histogram->SetBinError(1, v2ArrayStatError[0][1]);  
    v2Histogram->SetBinError(3, v2ArrayStatError[1][1]);
    v2Histogram->SetBinError(5, v2ArrayStatError[2][1]);    
        
    str1 = path + v2 + fileType;     
    v2Histogram->Draw("E");     
    v2Canvas->SaveAs(str1);    
    
    
    TString christinaCheck = "check";
    TString christinaCheck1 = "check1";
    TString christinaCheck2 = "check2";
    
    
    double mixedIntegralSBL;
    double mixedIntegralSBR;
    double mixedIntegralAVG;
    double siblingIntegral;
    
    double mixedIntegralUS;
    
    TCanvas *mixCheck = new TCanvas("c3", "Histograms", 3700, 850);
    mixCheck->Divide(3,1);
    
    TH2D *mixingCheck = new TH2D("mixing check", "mixing check", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *mixingCheckSBL = new TH2D("SBL", "SBL", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *mixingCheckSBR = new TH2D("SBR", "SBR", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *mixingCheckAVG = new TH2D("avg1", "avg1", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    TH2D *SBAVG = new TH2D("avg", "avg", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
    
    //         band    pt      vz      multiplcity
    //         left    int     5          10
    
    //siblingIntegral = sibCorrBin[0][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    mixedIntegralSBL = mixCorrBin[0][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    mixedIntegralSBR = mixCorrBin[2][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    mixedIntegralUS  = mixCorrBin[1][3][5][10]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
    
    SBAVG->Add(mixCorrBin[0][3][5][10], mixCorrBin[2][3][5][10], .5, .5);
    
    mixedIntegralAVG = SBAVG->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                
    //sibCorrBin[0][3][5][10]->Scale(1/siblingIntegral);     
    mixCorrBin[0][3][5][10]->Scale(1/mixedIntegralSBL);
    mixCorrBin[2][3][5][10]->Scale(1/mixedIntegralSBR);
    SBAVG->Scale(1/mixedIntegralAVG);
    mixCorrBin[1][3][5][10]->Scale(1/mixedIntegralUS);
    
   // mixingCheck->Divide(sibCorrBin[0][3][5][10], mixCorrBin[0][3][5][10], 1, 1);
    
   // mixingCheck->Draw("SURF1");
    
    mixingCheckSBL->Divide(mixCorrBin[0][3][5][10], mixCorrBin[1][3][5][10], 1, 1);
    mixingCheckSBR->Divide(mixCorrBin[2][3][5][10], mixCorrBin[1][3][5][10], 1, 1);
    mixingCheckAVG->Divide(SBAVG, mixCorrBin[1][3][5][10], 1, 1);
    
    mixCheck->cd(1);
    mixingCheckSBL->Draw("SURF1");
    mixCheck->cd(2);
    mixingCheckSBR->Draw("SURF1");
    mixCheck->cd(3);
    mixingCheckAVG->Draw("SURF1");
    
    str1 = path + christinaCheck + fileType;
    mixCheck->SaveAs(str1);        
        
    file->Close();
    output->Close();
   
}

}