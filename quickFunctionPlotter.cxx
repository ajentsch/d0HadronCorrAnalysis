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



int quickFunctionPlotter(){


    /* TF1 * function = new TF1("func","gaus", -1.5, 1.5);
    
    function->FixParameter(0, .4);
    function->FixParameter(1, 0);
    function->FixParameter(2, .5); */
    
    function  = new TF1("effWeightPions", "[0]*TMath::Exp(-([1]/x)**[2])", 0, 6);
  
    function->SetParameter(0, .809);
    function->SetParameter(1, .109);
    function->SetParameter(2, 3.224);
    
    
    TF1 * function2 = new TF1("withHFT", "[0]*TMath::Exp(-([1]/x)**[2])", 0, 6);
    
    function2->SetParameter(0, .201);
    function2->SetParameter(1, .109);
    function2->SetParameter(2, 3.224);
    
    function3  = new TF1("effWeightKaons", "[0]*TMath::Exp(-([1]/x)**[2])", 0, 6);
  
    function3->SetParameter(0, .809);
    function3->SetParameter(1, .231);
    function3->SetParameter(2, 3.968);
    //function3->SetParameter(3, .152);
    
    //function2->FixParameter(0, .01); //amplitude
    //function2->FixParameter(0,50);   //degrees of freedom
    //function2->FixParameter(1, .01);  //amplitude
    //function2->FixParameter(2, 0.0); //mean
    
    function2->SetLineColor(3);
    function3->SetLineColor(4);
    
    function->Draw();
    function2->Draw("SAME");
    function3->Draw("SAME");





















}