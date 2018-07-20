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

plotv2Comparison(){

    //v2 values
    double v2TwoPartAngCorr2D[11] = {.035, .050, .055, .065, .069, .067, .064, .052, .036, .019, .002};
    double v2SquareTwoPartCumulant[9] = { .006, .006, .0062, .0062, .0056, .0043, .0027, .0014, .0006};
    //double v2TwoPartCumulant[9];
    double v2FourthFourPartCumulant[8] = { .000007, .000014, .0000175, .000016, .000009, .000003, .0000005, 0.0};
    //double v2FourPartCumulant[8];
    //for(int i = 0; i < 9; i++){ v2TwoPartCumulant[i] = TMath::Sqrt(v2SquareTwoPartCumulant[i]); cout << v2TwoPartCumulant[i] << endl;}
    //for(int i = 0; i < 8; i++){ v2FourPartCumulant[i] = TMath::Sqrt(TMath::Sqrt(v2FourthFourPartCumulant[i])); cout << v2FourPartCumulant[i] << endl;}
    
    double v2TwoPartCumulant[9] = { .025, .035, .05, .065, .074, .077, .076, .074, .0745};
    double v2FourPartCumulant[8] = {.024, .04, .055, .062, .063, .06, .055, .064};
    
    //centrality bins
    double centLanny[11] = {.89, .79, .69, .59, .49, .42, .33, .23, .14, .07, .025};
    //double cent2PartCumulant[9] = {.8, .76, .55, .45, .33, .15, .07, .025};
    //double cent4PartCumulant[8] = {.65, .51, .47, .31, .21, .09, .025 };
    
    double cent2PartCumulant[9] = {.025, .075, .15, .25, .35, .45, .54, .65, .75 };
    double cent4PartCumulant[8] = {.075, .15, .25, .35, .45, .55, .65, .75 };
    
    double lannyAQValues[11] = {.002, .011, .028, .070, .136, .201, .270, .268, .179, .063, .001};
    double lannyErrPercent[11] = {.01, .01, .01, .01, .02, .04, .04, .05, .06, .09, .08};
    double lannydNdEta[11] = {5.2, 13.9, 28.8, 52.8, 89, 139, 209, 307, 440, 564, 671};
    double lannyErrors[11];
    
    for(int i = 0; i < 11; i++){
    
        //lannyErrors[i] = TMath::Sqrt((TMath::Pi()*lannyAQValues[i]*lannyErrPercent[i])/lannydNdEta[i]);
        
        lannyErrors[i] = TMath::Sqrt(TMath::Pi()/lannydNdEta[i])*.5*(1.0/TMath::Sqrt(lannyAQValues[i]))*(lannyErrPercent[i]*lannyAQValues[i]);
    }
    
    double errorsX1[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double errorsX2[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double errorsX3[8] = { 0, 0, 0, 0, 0, 0, 0, 0};
    
    double twoPartCumulantErrors[9] = { .000096, .000088, .000066, .000075, .000094, .00013, .00018, .00029, .00052};
    double fourPartCumulantErrors[8] = {.00066, .00022, .0002, .00026, .00045, .001, .0035, 0.01 };
    
    TGraphErrors * v2LannyPlot = new TGraphErrors(11, centLanny, v2TwoPartAngCorr2D, errorsX1, lannyErrors);
    TGraphErrors * v2TwoPartCumulantPlot = new TGraphErrors(9, cent2PartCumulant,v2TwoPartCumulant, errorsX2, twoPartCumulantErrors);
    TGraphErrors * v2FourPartCumulantPlot = new TGraphErrors(8, cent4PartCumulant, v2FourPartCumulant, errorsX3, fourPartCumulantErrors );
    
    
    v2LannyPlot->SetMarkerColor(2);
    v2LannyPlot->SetMarkerStyle(21);
    v2TwoPartCumulantPlot->SetMarkerColor(4);
    v2TwoPartCumulantPlot->SetMarkerStyle(19);
    v2FourPartCumulantPlot->SetMarkerColor(6);
    v2FourPartCumulantPlot->SetMarkerStyle(19);
    
    //TLegend * leg = new TLegend( .22, .12, .89 , .35);
    //leg->AddEntry(v2LannyPlot, "v2{2D} (multi-parameter fit -> method/data in https://arxiv.org/pdf/1109.4380.pdf),Table III", "p");
    //leg->AddEntry(v2TwoPartCumulantPlot, "v2{2} (event-plane) (method/data https://arxiv.org/pdf/nucl-ex/0409033.pdf, fig. 30)", "p");
    //leg->AddEntry(v2FourPartCumulantPlot, "v2{4} (method/data https://arxiv.org/pdf/nucl-ex/0409033.pdf, fig. 30)", "p");
    
    TCanvas * can = new TCanvas("","", 1100, 1000);
    
    v2TwoPartCumulantPlot->SetMinimum(0);
    v2TwoPartCumulantPlot->SetTitle("v_{2} comparisons");
    v2TwoPartCumulantPlot->GetXaxis()->SetTitle("Centrality");
    v2TwoPartCumulantPlot->GetYaxis()->SetTitle("v_{2}");
    v2TwoPartCumulantPlot->GetYaxis()->SetTitleOffset(1.4);
    v2TwoPartCumulantPlot->GetYaxis()->CenterTitle();
    //v2TwoPartCumulantPlot->Draw("APC");
    //v2FourPartCumulantPlot->Draw("SAME PC");
    //v2LannyPlot->Draw("SAME PC");
    //leg->Draw("SAME");
    
	
    //can->SaveAs("v_2_comparison_plot.png");


	TH1D* testGraph = new TH1D("testing", "testing", 12, 1, 13);
	//testGraph->GetXaxis()->SetBinLabel(1, "PYTHIA");
	//testGraph->GetXaxis()->SetBinLabel(2, "");
	//testGraph->GetXaxis()->SetBinLabel(3, "");
	//testGraph->GetXaxis()->SetBinLabel(4, "");
	//testGraph->GetXaxis()->SetBinLabel(5, "");
	//testGraph->GetXaxis()->SetBinLabel(6, "");
	//testGraph->GetXaxis()->SetBinLabel(7, "");
	//testGraph->GetXaxis()->SetBinLabel(8, "");
	//testGraph->GetXaxis()->SetBinLabel(9, "");
	//testGraph->GetXaxis()->SetBinLabel(10, "");
	//testGraph->GetXaxis()->SetBinLabel(11, "");
	//testGraph->GetXaxis()->SetBinLabel(12, "");
	//testGraph->GetXaxis()->SetBinLabel(3, "");
	
	
	double centralityWithPythia[8] = {6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.25, 12.75};
	double centWithPythiaLow[8]    = {.5, .5, .5, .5, .5, .5, .25, .25};
	double centWithPythiaHigh[8]   = {.5, .5, .5, .5, .5, .5, .25, .25};
	
	double kettlerPt1EtaData[8]      = {.45, .65, .95, 1.4, 1.8, 2.1, 2.45, 2.1};
	double kettlerPt1EtaErrorLow[8]  = {.05*.45, .05*.65, .09*.95, .09*1.4, .09*1.8, .09*2.1, .09*2.45, .09*2.1};
	double kettlerPt1EtaErrorHigh[8] = {.05*.45, .05*.65, .09*.95, .09*1.4, .09*1.8, .09*2.1, .09*2.45, .09*2.1};
	
	double kettlerPt3EtaData[8]      = {.5, .48, .4, .83, .65, .55, 1.85, 1.35};
	double kettlerPt3EtaErrorLow[8]  = {.1*.5, .1*.48, .15*.4, .15*.83, .15*.65, .17*.55, .17*1.85, .17*1.35};
	double kettlerPt3EtaErrorHigh[8] = {.1*.5, .1*.48, .15*.4, .15*.83, .15*.65, .17*.55, .17*1.85, .17*1.35};
	
	double kettlerPt1PhiData[8]      = {.43, .43, .49, .53, .58, .59, .57, .54};
	double kettlerPt1PhiErrorLow[8]  = {.05*.43, .05*.43, .05*.49, .05*.53, .05*.58, .05*.59, .05*.57, .05*.54};
	double kettlerPt1PhiErrorHigh[8] = {.05*.43, .05*.43, .05*.49, .05*.53, .05*.58, .05*.59, .05*.57, .05*.54};
	
	double kettlerPt3PhiData[8]      = {.45, .52, .37, .39, .35, .38, .59, .42};
	double kettlerPt3PhiErrorLow[8]  = {.12*.45, .12*.52, .11*.37, .11*.39, .11*.35, .1*.38, .1*.59, .1*.42};
	double kettlerPt3PhiErrorHigh[8] = {.12*.45, .12*.52, .11*.37, .11*.39, .11*.35, .1*.38, .1*.59, .1*.42};
	
	TGraphAsymmErrors *kettlerPt1Eta = new TGraphAsymmErrors(8, centralityWithPythia, kettlerPt1EtaData, centWithPythiaLow, centWithPythiaHigh, kettlerPt1EtaErrorLow, kettlerPt1EtaErrorHigh);
	TGraphAsymmErrors *kettlerPt3Eta = new TGraphAsymmErrors(8, centralityWithPythia, kettlerPt3EtaData, centWithPythiaLow, centWithPythiaHigh, kettlerPt3EtaErrorLow, kettlerPt3EtaErrorHigh);
	TGraphAsymmErrors *kettlerPt1Phi = new TGraphAsymmErrors(8, centralityWithPythia, kettlerPt1PhiData, centWithPythiaLow, centWithPythiaHigh, kettlerPt1PhiErrorLow, kettlerPt1PhiErrorHigh);
	TGraphAsymmErrors *kettlerPt3Phi = new TGraphAsymmErrors(8, centralityWithPythia, kettlerPt3PhiData, centWithPythiaLow, centWithPythiaHigh, kettlerPt3PhiErrorLow, kettlerPt3PhiErrorHigh);
	
	double centralityAlexWithPythia[3] = {6.5, 9.5, 12.0};
	double centAlexWithPythiaLow[3]  = {1.5, 1.5, 1.0};
	double centAlexWithPythiaHigh[3] = {1.5, 1.5, 1.0};
	
	double alexEtaData[3]    = {.312136, 1.40393, 1.24592};
	double alexEtaErrorLow[3] = {.112368, .353493, .221779};
	double alexEtaErrorHigh[3] = {.112368, .353493, .221779};
	
	double alexPhiData[3]    = {.355024, .669247, .754439};
	double alexPhiErrorLow[3] = {.0675623, .0618335, .0807210};
	double alexPhiErrorHigh[3] = {.0675623, .0618335, .0807210};
	
	double alexNSYieldData[3]      = {.371115, 4.95248, 17.4999};
	double alexNSYieldErrorLow[3]  = {.138868, 1.8081, 7.16079};
	double alexNSYieldErrorHigh[3] = { .138868, 1.8081, 7.16079};
	
	double alexQuadData[3]      = {.00463036, .00647402, 0.0};
	double alexQuadErrorLow[3]  = {.00268507, .00320154 , .000866073};
	double alexQuadErrorHigh[3] = {.00268507, .00320154 , .000866073};
	
	double centralityPythia[1] = {2.0};
	double centPythiaLow[1]  = {1.0};
	double centPythiaHigh[1] = {1.0};
	
	double pythiaEtaData[1] = {.330681};
	double pythiaEtaErrorLow[1] = {.0310280};
	double pythiaEtaErrorHigh[1] = {.0310280};
	
	double pythiaPhiData[1] = { .310772};
	double pythiaPhiErrorLow[1] = {.0343975};
	double pythiaPhiErrorHigh[1] = {.0343975};
	
	//pythiaResultsJetAmp[2]->SetBinContent(1, 2.21967);
    //pythiaResultsJetAmp[2]->SetBinError(1, .764791);
       
    //pythiaResultsEtaWidth[2]->SetBinContent(1, .330681);
    //pythiaResultsEtaWidth[2]->SetBinError(1,    .0310280);
        
    //pythiaResultsPhiWidth[2]->SetBinContent(1, .310772);
    //pythiaResultsPhiWidth[2]->SetBinError(1,   .0343975);
	
	//fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(2, .355024);
    //fitParameterSigmaPhiPlotVsCent[2]->SetBinError(2, .0675623);
	//fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(3, .669247);
    //fitParameterSigmaPhiPlotVsCent[2]->SetBinError(3, .0618335);
	//fitParameterSigmaPhiPlotVsCent[2]->SetBinContent(4,    .754439);
    //fitParameterSigmaPhiPlotVsCent[2]->SetBinError(4,       .0807210);
	
	//fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(2, .306846);
    //fitParameterSigmaEtaPlotVsCent[2]->SetBinError(2, .112368);
	//fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(3, 1.40393);
    //fitParameterSigmaEtaPlotVsCent[2]->SetBinError(3, .353493);
	//fitParameterSigmaEtaPlotVsCent[2]->SetBinContent(4,     1.24592);
    //fitParameterSigmaEtaPlotVsCent[2]->SetBinError(4,       .221779);
	
	TGraphAsymmErrors *alexEta = new TGraphAsymmErrors(3, centralityAlexWithPythia, alexEtaData, centAlexWithPythiaLow, centAlexWithPythiaHigh, alexEtaErrorLow, alexEtaErrorHigh);
	TGraphAsymmErrors *alexPhi = new TGraphAsymmErrors(3, centralityAlexWithPythia, alexPhiData, centAlexWithPythiaLow, centAlexWithPythiaHigh, alexPhiErrorLow, alexPhiErrorHigh);
	
	TGraphAsymmErrors *pythiaEta = new TGraphAsymmErrors(1, centralityPythia, pythiaEtaData, centPythiaLow, centPythiaHigh, pythiaEtaErrorLow, pythiaEtaErrorHigh);
	TGraphAsymmErrors *pythiaPhi = new TGraphAsymmErrors(1, centralityPythia, pythiaPhiData, centPythiaLow, centPythiaHigh, pythiaPhiErrorLow, pythiaPhiErrorHigh);
	
	
	////////Formatting
	
	pythiaEta->SetMarkerStyle(22);
    pythiaEta->SetMarkerSize(3);
    pythiaEta->SetMarkerColor(8);
    pythiaEta->SetLineWidth(3);
    pythiaEta->SetLineColor(8);
	
	pythiaPhi->SetMarkerStyle(22);
    pythiaPhi->SetMarkerSize(3);
    pythiaPhi->SetMarkerColor(8);
    pythiaPhi->SetLineWidth(3);
    pythiaPhi->SetLineColor(8);
	
	alexEta->SetMarkerStyle(29);
    alexEta->SetMarkerSize(5);
    alexEta->SetLineWidth(5);
	alexEta->SetMarkerColor(2);
	alexEta->SetLineColor(2);
    alexEta->SetMaximum(2.5);
	alexEta->SetTitle("");
	alexEta->GetYaxis()->SetTitle("#sigma_{NS, #Delta#eta}");
	alexEta->GetYaxis()->SetTitleSize(.06);
	alexEta->GetYaxis()->SetTitleOffset(1.0);
	
	alexPhi->SetMarkerStyle(29);
    alexPhi->SetMarkerSize(5);
    alexPhi->SetLineWidth(5);
	alexPhi->SetMarkerColor(2);
	alexPhi->SetLineColor(2);
    alexPhi->SetMaximum(1.2);
	alexPhi->SetTitle("");
	alexPhi->GetYaxis()->SetTitle("#sigma_{NS, #Delta#phi}");
	alexPhi->GetYaxis()->SetTitleSize(.06);
	alexPhi->GetYaxis()->SetTitleOffset(1.0);
    
	kettlerPt1Eta->SetMarkerStyle(9);
    kettlerPt1Eta->SetMarkerSize(2);
    kettlerPt1Eta->SetMarkerColor(38);
    kettlerPt1Eta->SetLineWidth(4);
    kettlerPt1Eta->SetLineColor(38);
	
	kettlerPt3Eta->SetMarkerStyle(9);
    kettlerPt3Eta->SetMarkerSize(2);
    kettlerPt3Eta->SetMarkerColor(4);
    kettlerPt3Eta->SetLineWidth(4);
    kettlerPt3Eta->SetLineColor(4);
	
	kettlerPt1Phi->SetMarkerStyle(9);
    kettlerPt1Phi->SetMarkerSize(2);
    kettlerPt1Phi->SetMarkerColor(38);
    kettlerPt1Phi->SetLineWidth(4);
    kettlerPt1Phi->SetLineColor(38);
	
	kettlerPt3Phi->SetMarkerStyle(9);
    kettlerPt3Phi->SetMarkerSize(2);
    kettlerPt3Phi->SetMarkerColor(4);
    kettlerPt3Phi->SetLineWidth(4);
    kettlerPt3Phi->SetLineColor(4);
	
	/////////////////
	
	
    
    double etaWidthValues[3] = {.312136, 1.40393, 1.24592};
    double etaWidthSystematics[3] = { .032325296, .200152368 , .129029439};
   
    //SYSTEMATICS ERRORS HERE
    
	
	
	testGraph->SetMaximum(2.7);
	testGraph->SetStats(0);
	
	TCanvas * etaCan = new TCanvas("","", 1100, 1000);
	
	//testGraph->Draw();
	alexEta->Draw("AP");
	alexEta->GetXaxis()->SetLimits(1,13);
	for(int i = 0; i < 3; i++) {
        double x1 = centralityAlexWithPythia[i]-.1;
        double x2 = centralityAlexWithPythia[i]+.1;
        double y1 = etaWidthValues[i]-etaWidthSystematics[i];
        double y2 = etaWidthValues[i]+etaWidthSystematics[i];

   
    
        TLine *la = new TLine(x1, y1, x1, y1+0.01);
        la->SetLineColor(2);
        la->SetLineWidth(5);
        la->Draw("same");
        TLine *lb = new TLine(x2, y1, x2, y1+0.01);
        lb->SetLineColor(2);
        lb->SetLineWidth(5);
        lb->Draw("same");
        TLine *lc = new TLine(x1, y2, x1, y2-0.01);
        lc->SetLineColor(2);
        lc->SetLineWidth(5);
        lc->Draw("same");
        TLine *ld = new TLine(x2, y2, x2, y2-0.01);
        ld->SetLineColor(2);
        ld->SetLineWidth(5);
        ld->Draw("same");
        TLine *le = new TLine(x1, y1, x2, y1);
        le->SetLineWidth(4);
        le->SetLineColor(2);
        le->Draw("same");
        TLine *lf = new TLine(x1, y2, x2, y2);
        lf->SetLineWidth(4);
        lf->SetLineColor(2);
        lf->Draw("same");
    }
	kettlerPt1Eta->Draw("SAME P");
	kettlerPt3Eta->Draw("SAME P");
	pythiaEta->Draw("SAME P");
	
	etaCan->SaveAs("test_delEta.png");

	
	double phiWidthValues[3] = {.349718, .669247, .754439}; //these are the nominal values
    double phiWidthSystematics[3] = { .037051958, .075790561 , .077030402};
    
   
    //SYSTEMATICS ERRORS HERE
    
	TCanvas * phiCan = new TCanvas("","", 1100, 1000);
	
	alexPhi->Draw("AP");
	alexPhi->GetXaxis()->SetLimits(1,13);

    for(int i = 0; i < 3; i++) {
        double x1 = centralityAlexWithPythia[i]-.1;
        double x2 = centralityAlexWithPythia[i]+.1;
        double y1 = phiWidthValues[i]-phiWidthSystematics[i];
        double y2 = phiWidthValues[i]+phiWidthSystematics[i];

   
    
        TLine *la = new TLine(x1, y1, x1, y1+0.01);
        la->SetLineColor(2);
        la->SetLineWidth(5);
        la->Draw("same");
        TLine *lb = new TLine(x2, y1, x2, y1+0.01);
        lb->SetLineColor(2);
        lb->SetLineWidth(5);
        lb->Draw("same");
        TLine *lc = new TLine(x1, y2, x1, y2-0.01);
        lc->SetLineColor(2);
        lc->SetLineWidth(5);
        lc->Draw("same");
        TLine *ld = new TLine(x2, y2, x2, y2-0.01);
        ld->SetLineColor(2);
        ld->SetLineWidth(5);
        ld->Draw("same");
        TLine *le = new TLine(x1, y1, x2, y1);
        le->SetLineWidth(4);
        le->SetLineColor(2);
        le->Draw("same");
        TLine *lf = new TLine(x1, y2, x2, y2);
        lf->SetLineWidth(4);
        lf->SetLineColor(2);
        lf->Draw("same");
    }
	
	kettlerPt1Phi->Draw("SAME P");
	kettlerPt3Phi->Draw("SAME P");
	pythiaPhi->Draw("SAME P");
	
	TLegend * leg = new TLegend(.17, .88, .89, .6);
    leg->AddEntry(kettlerPt3Phi, "di-Hadron, Mean Trigger p_{T} =  5.7 GeV/c [3]");
    leg->AddEntry(kettlerPt1Phi, "di-Hadron, Mean Trigger p_{T} =  2.56 GeV/c [3]");
    
    leg->AddEntry(pythiaPhi,         "Pythia D^{0}-Hadron, Mean D^{0} p_{T} = 3 GeV/c");
    leg->AddEntry(alexPhi, "D^{0}-Hadron AuAu 200 GeV, Mean D^{0} p_{T} = 3 GeV/c");
    
	leg->SetTextAlign(12);
	leg->SetFillColor(0);
    leg->SetMargin(.15);

	leg->SetBorderSize(0);
	leg->Draw("SAME");
	
	phiCan->SaveAs("test_delPhi.png");
	
}

