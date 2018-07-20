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

//#include "C:/Users/ajentsch/desktop/corrHistogramMaker.h"


int D0HadronCorrErrors(){



/***********************/
//This code functions by reading in a root file containing single particle 
//information for (kpi) pairs and hadrons. All (kpi) pairs are considered D0  
//candidates. The root file will read in the information, obtain the occupancies 
//that give the values for (m,n) in the error calculations, and then use these values 
//to calculate the bin-by-bin statisical errors in the (delEta, delPhi) bins of the
//final correlation quantities.
//
//
/************************/


//int numEventsWithD0 = Integral of invariant mass distribution in mass range


//A - rhoSib from mass band
//B - rhoRef D0 - not measured, but can be estimated
//C - rhoRef [kpi] - event mixing from mass band
//D - 1 + side band correlation



}