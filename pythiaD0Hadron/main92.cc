// This version of main92.cc is being used for D0-hadron correlations.
//
// Author: Alex Jentsch 
//
// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"

double calcInvMass(double, double, double, double, double, double, double, double);
int getPtBin(double);

using namespace Pythia8;

int main() {
  
    //Constants
    double piMass = .13957018;
    double kMass  = .493677;  
    double phiBinShift = (TMath::Pi()/12.0); 
    int NUM_MIXED_EVENTS = 20;
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 13;

    bool DEBUG = false;
    bool removeDStarContribution = false;

    // Create Pythia instance and set it up to generate hard QCD processes
    Pythia pythia;
    pythia.readFile("star_hf_tune.cmnd");
    pythia.init();

    // Set up the ROOT TFile and TTree.
    TFile *file = TFile::Open("pytree.root","recreate");
    Event *event = &pythia.event;
    //TTree *T = new TTree("T","ev1 Tree");
    //T->Branch("event",&event);

    //setup histograms for QA and testing
    TH1D * trackMult       = new TH1D("track_multiplicity", "track_multiplicity", 400, 1, 401);
    TH1D * hadronPtDist    = new TH1D("hadron_pt", "hadron_pt", 400, 0, 10);
    TH1D * hadronEtaDist   = new TH1D("hadron_eta", "hadron_eta", 400, -2, 2);
    TH1D * hadronPhiDist   = new TH1D("hadron_phi", "hadron_phi", 400, -TMath::TwoPi(), TMath::TwoPi());
    TH1D * D0PtDist        = new TH1D("D0Pt", "D0Pt", 100, 0, 10);
    TH1D * D0EtaDist       = new TH1D("D0Eta", "D0Eta", 100, -2, 2);
    TH1D * D0PhiDist       = new TH1D("D0Phi", "D0Phi", 100, -TMath::TwoPi(), TMath::TwoPi());

    
    TH1D * piPlusPtDist    = new TH1D("PiPlusPt","PiPlusPt", 100, 0, 10);
    TH1D * piMinusPtDist   = new TH1D("PiMiunsPt", "PiMinusPt", 100, 0, 10);
    TH1D * D0Mass          = new TH1D("D0Mass", "D0Mass", 100, 1.6, 2.2);
    TH1D * D0BarMass       = new TH1D("D0BarMass", "D0BarMass", 100, 1.6, 2.2);
    TH1D * DStarMass       = new TH1D("DStar_plus_mass", "DStar_plus_mass", 100, 1.6, 2.2);
    TH1D * DStarBarMass    = new TH1D("DStar_minus_mass", "DStar_minus_mass", 100, 1.6, 2.2);
    TH1D * kPiInvMass      = new TH1D("KPiInvMass", "KPiInvMass", 100, 1.7, 2.1);
    
    TH1D * totalD0Yield    = new TH1D("D0_Yield", "D0_Yield", 5, 1, 6);
    totalD0Yield->GetXaxis()->SetBinLabel(1, "0-1");
    totalD0Yield->GetXaxis()->SetBinLabel(2, "1-2");
    totalD0Yield->GetXaxis()->SetBinLabel(3, "2-3");
    totalD0Yield->GetXaxis()->SetBinLabel(4, "3-4");
    totalD0Yield->GetXaxis()->SetBinLabel(5, "4-10");
    
    
    TH1D * totalDStarToD0    = new TH1D("DStar_to_D0", "DStar_to_D0", 5, 1, 6);
    totalDStarToD0->GetXaxis()->SetBinLabel(1, "0-1");
    totalDStarToD0->GetXaxis()->SetBinLabel(2, "1-2");
    totalDStarToD0->GetXaxis()->SetBinLabel(3, "2-3");
    totalDStarToD0->GetXaxis()->SetBinLabel(4, "3-4");
    totalDStarToD0->GetXaxis()->SetBinLabel(5, "4-10");
    
    TH1D * kPiInvMassPtBin[5];

    //correlation histograms go here
    TH2D * sibDHadronCorr[5];

    sibDHadronCorr[0] = new TH2D("sibling_D0_Hadron_Corr_0_1GeV", "sibling_D0_Hadron_Corr_0_1_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorr[1] = new TH2D("sibling_D0_Hadron_Corr_1_2GeV", "sibling_D0_Hadron_Corr_1_2_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorr[2] = new TH2D("sibling_D0_Hadron_Corr_2_3GeV", "sibling_D0_Hadron_Corr_2_3_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorr[3] = new TH2D("sibling_D0_Hadron_Corr_3_4GeV", "sibling_D0_Hadron_Corr_3_4_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorr[4] = new TH2D("sibling_D0_Hadron_Corr_4_10GeV", "sibling_D0_Hadron_Corr_4_10_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH2D * sibDHadronCorrD0SoftPionOnly[5];

    sibDHadronCorrD0SoftPionOnly[0] = new TH2D("sibling_D0_Soft_Pion_Peak_0_1GeV", "sibling_D0_Soft_Pion_Peak_0_1_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrD0SoftPionOnly[1] = new TH2D("sibling_D0_Soft_Pion_Peak_1_2GeV", "sibling_D0_Soft_Pion_Peak_1_2_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrD0SoftPionOnly[2] = new TH2D("sibling_D0_Soft_Pion_Peak_2_3GeV", "sibling_D0_Soft_Pion_Peak_2_3_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrD0SoftPionOnly[3] = new TH2D("sibling_D0_Soft_Pion_Peak_3_4GeV", "sibling_D0_Soft_Pion_Peak_3_4_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrD0SoftPionOnly[4] = new TH2D("sibling_D0_Soft_Pion_Peak_4_10GeV", "sibling_D0_Soft_Pion_Peak_4_10_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    TH1D * sibDHadronCorrPhi[5];

    sibDHadronCorrPhi[0] = new TH1D("sibling_D0_Hadron_Corr_0_1GeVPhi", "sibling_D0_Hadron_Corr_0_1_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrPhi[1] = new TH1D("sibling_D0_Hadron_Corr_1_2GeVPhi", "sibling_D0_Hadron_Corr_1_2_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrPhi[2] = new TH1D("sibling_D0_Hadron_Corr_2_3GeVPhi", "sibling_D0_Hadron_Corr_2_3_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrPhi[3] = new TH1D("sibling_D0_Hadron_Corr_3_4GeVPhi", "sibling_D0_Hadron_Corr_3_4_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrPhi[4] = new TH1D("sibling_D0_Hadron_Corr_4_10GeVPhi", "sibling_D0_Hadron_Corr_4_10_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH2D * mixDHadronCorr[5];

    mixDHadronCorr[0] = new TH2D("mixed_D0_Hadron_Corr_0_1GeV", "mixed_D0_Hadron_Corr_0_1_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorr[1] = new TH2D("mixed_D0_Hadron_Corr_1_2GeV", "mixed_D0_Hadron_Corr_1_2_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorr[2] = new TH2D("mixed_D0_Hadron_Corr_2_3GeV", "mixed_D0_Hadron_Corr_2_3_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorr[3] = new TH2D("mixed_D0_Hadron_Corr_3_4GeV", "mixed_D0_Hadron_Corr_3_4_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorr[4] = new TH2D("mixed_D0_Hadron_Corr_4_10GeV", "mixed_D0_Hadron_Corr_4_10_GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    TH2D * ratioDHadronCorr[5];

    ratioDHadronCorr[0] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_0_1GeV", "RhoSib_Over_RhoMix_D0_Hadron_Corr_0_1GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorr[1] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_1_2GeV", "RhoSib_Over_RhoMix_D0_Hadron_Corr_1_2GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorr[2] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_2_3GeV", "RhoSib_Over_RhoMix_D0_Hadron_Corr_2_3GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorr[3] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_3_4GeV", "RhoSib_Over_RhoMix_D0_Hadron_Corr_3_4GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorr[4] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_4_10GeV", "RhoSib_Over_RhoMix_D0_Hadron_Corr_4_10GeV", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    //NO SOFT PION BELOW HERE
    
    TH2D * sibDHadronCorrNoSoftPion[5];

    sibDHadronCorrNoSoftPion[0] = new TH2D("sibling_D0_Hadron_Corr_0_1GeV_NoSoftPion", "sibling_D0_Hadron_Corr_0_1_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPion[1] = new TH2D("sibling_D0_Hadron_Corr_1_2GeV_NoSoftPion", "sibling_D0_Hadron_Corr_1_2_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPion[2] = new TH2D("sibling_D0_Hadron_Corr_2_3GeV_NoSoftPion", "sibling_D0_Hadron_Corr_2_3_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPion[3] = new TH2D("sibling_D0_Hadron_Corr_3_4GeV_NoSoftPion", "sibling_D0_Hadron_Corr_3_4_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPion[4] = new TH2D("sibling_D0_Hadron_Corr_4_10GeV_NoSoftPion", "sibling_D0_Hadron_Corr_4_10_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    TH1D * sibDHadronCorrNoSoftPionPhi[5];
    
    sibDHadronCorrNoSoftPionPhi[0] = new TH1D("sibling_D0_Hadron_Corr_NoSoftPion_0_1GeVPhi", "sibling_D0_Hadron_Corr_NoSoftPion_0_1_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPionPhi[1] = new TH1D("sibling_D0_Hadron_Corr_NoSoftPion_1_2GeVPhi", "sibling_D0_Hadron_Corr_NoSoftPion_1_2_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPionPhi[2] = new TH1D("sibling_D0_Hadron_Corr_NoSoftPion_2_3GeVPhi", "sibling_D0_Hadron_Corr_NoSoftPion_2_3_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPionPhi[3] = new TH1D("sibling_D0_Hadron_Corr_NoSoftPion_3_4GeVPhi", "sibling_D0_Hadron_Corr_NoSoftPion_3_4_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    sibDHadronCorrNoSoftPionPhi[4] = new TH1D("sibling_D0_Hadron_Corr_NoSoftPion_4_10GeVPhi", "sibling_D0_Hadron_Corr_NoSoftPion_4_10_GeVPhi", NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    TH2D * mixDHadronCorrNoSoftPion[5];

    mixDHadronCorrNoSoftPion[0] = new TH2D("mixed_D0_Hadron_Corr_0_1GeV_NoSoftPion", "mixed_D0_Hadron_Corr_0_1_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorrNoSoftPion[1] = new TH2D("mixed_D0_Hadron_Corr_1_2GeV_NoSoftPion", "mixed_D0_Hadron_Corr_1_2_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorrNoSoftPion[2] = new TH2D("mixed_D0_Hadron_Corr_2_3GeV_NoSoftPion", "mixed_D0_Hadron_Corr_2_3_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorrNoSoftPion[3] = new TH2D("mixed_D0_Hadron_Corr_3_4GeV_NoSoftPion", "mixed_D0_Hadron_Corr_3_4_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    mixDHadronCorrNoSoftPion[4] = new TH2D("mixed_D0_Hadron_Corr_4_10GeV_NoSoftPion", "mixed_D0_Hadron_Corr_4_10_GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    TH2D * ratioDHadronCorrNoSoftPion[5];

    ratioDHadronCorrNoSoftPion[0] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_0_1GeV_NoSoftPion", "RhoSib_Over_RhoMix_D0_Hadron_Corr_0_1GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorrNoSoftPion[1] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_1_2GeV_NoSoftPion", "RhoSib_Over_RhoMix_D0_Hadron_Corr_1_2GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorrNoSoftPion[2] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_2_3GeV_NoSoftPion", "RhoSib_Over_RhoMix_D0_Hadron_Corr_2_3GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorrNoSoftPion[3] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_3_4GeV_NoSoftPion", "RhoSib_Over_RhoMix_D0_Hadron_Corr_3_4GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    ratioDHadronCorrNoSoftPion[4] = new TH2D("RhoSib_Over_RhoMix_D0_Hadron_Corr_4_10GeV_NoSoftPion", "RhoSib_Over_RhoMix_D0_Hadron_Corr_4_10GeV_NoSoftPion", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);

    //SUM W2 HERE
    
    sibDHadronCorr[0]->Sumw2();
    sibDHadronCorr[1]->Sumw2();
    sibDHadronCorr[2]->Sumw2();
    sibDHadronCorr[3]->Sumw2();
    sibDHadronCorr[4]->Sumw2();
    sibDHadronCorrPhi[0]->Sumw2();
    sibDHadronCorrPhi[1]->Sumw2();
    sibDHadronCorrPhi[2]->Sumw2();
    sibDHadronCorrPhi[3]->Sumw2();
    sibDHadronCorrPhi[4]->Sumw2();
    mixDHadronCorr[0]->Sumw2();
    mixDHadronCorr[1]->Sumw2();
    mixDHadronCorr[2]->Sumw2();
    mixDHadronCorr[3]->Sumw2();
    mixDHadronCorr[4]->Sumw2();
    
    sibDHadronCorrD0SoftPionOnly[0]->Sumw2();
    sibDHadronCorrD0SoftPionOnly[1]->Sumw2();
    sibDHadronCorrD0SoftPionOnly[2]->Sumw2();
    sibDHadronCorrD0SoftPionOnly[3]->Sumw2();
    sibDHadronCorrD0SoftPionOnly[4]->Sumw2();
    sibDHadronCorrNoSoftPion[0]->Sumw2();
    sibDHadronCorrNoSoftPion[1]->Sumw2();
    sibDHadronCorrNoSoftPion[2]->Sumw2();
    sibDHadronCorrNoSoftPion[3]->Sumw2();
    sibDHadronCorrNoSoftPion[4]->Sumw2();
    sibDHadronCorrNoSoftPionPhi[0]->Sumw2();
    sibDHadronCorrNoSoftPionPhi[1]->Sumw2();
    sibDHadronCorrNoSoftPionPhi[2]->Sumw2();
    sibDHadronCorrNoSoftPionPhi[3]->Sumw2();
    sibDHadronCorrNoSoftPionPhi[4]->Sumw2();
    
    
    kPiInvMassPtBin[0] = new TH1D("0-1", "0-1", 100, 1.7, 2.1);
    kPiInvMassPtBin[1] = new TH1D("1-2", "1-2", 100, 1.7, 2.1);
    kPiInvMassPtBin[2] = new TH1D("2-3", "2-3", 100, 1.7, 2.1);
    kPiInvMassPtBin[3] = new TH1D("3-4", "3-4", 100, 1.7, 2.1);
    kPiInvMassPtBin[4] = new TH1D("4-10", "4-10", 100, 1.7, 2.1);

    std::vector<Vec4> hadrons;
    std::vector<Vec4> hadronsNoSoftPion;
    std::vector<Vec4> DMesons;
    std::vector<Vec4> DMesonsForMixing;
    std::vector<Vec4> DMesonsForMixingNoSoftPion;
    
    std::vector<Vec4> softPionOnly;

    Vec4 mixerHadrons[20][300];        // CHANGE ARRAY ARGUMENTS TO CHANGE # OF MIXED EVENTS
    int  numMixerHadrons[20];          //
    
    Vec4 mixerHadronsNoSoftPion[20][100];        // CHANGE ARRAY ARGUMENTS TO CHANGE # OF MIXED EVENTS
    int  numMixerHadronsNoSoftPion[20];          //
    
    std::vector<Vec4> piNeg;
    std::vector<Vec4> piPos;
    std::vector<Vec4> kNeg;
    std::vector<Vec4> kPos;
    std::vector<int> piNegMotherIdx;
    std::vector<int> piPosMotherIdx;
    std::vector<int> kNegMotherIdx;
    std::vector<int> kPosMotherIdx;

    std::vector<Vec4> charm;
    std::vector<Vec4> antiCharm;
    
    int nEvents = pythia.mode("Main:numberOfEvents");
   
    bool hasDtoKPi = false;
    bool dStarIsMother = false;
    bool isFinalState = false;

    bool mixerBufferReady = false;
    int mixerBufferIdx = 0;
    
    double p1x, p1y, p1z, p2x, p2y, p2z, pt, delEta, delPhi, delPhiCp;
    int ptBin;
    Vec4 fourMom;
    Vec4 temp;
    Vec4 fourMom1, fourMom2;
    temp.reset();
    int numD0ToKPi = 0;
    int D0Idx;
    int D0PtBin;
    int numSoftPion = 0;

    
    //D0 checking information
    int trkPID, motherPID, daughter1Idx, daughter2Idx, daughterPID1, daughterPID2, motherIdx, motherIdxD0;

    int numFinalStateTracks = 0;
    int numFinalStateTracksNoSoftPion = 0;

    //PID Codes
    //pi+ = 211
    //pi- = -211
    //K+  = 321
    //K-  = -321
    //p+  = 2212
    //P-  = -2212
    //D0  = 421
    //D*+ = 413

    double etaC, etaCBar, phiC, phiCBar;

    TH2D* CAndCBarCorr = new TH2D("c_cbar_corr", "c_cbar_corr", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    TH2D* DAndDBarCorr = new TH2D("D0_D0bar_corr", "D0_D0bar_corr", NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2()+phiBinShift, 3*TMath::PiOver2()+phiBinShift);
    
    
    //--------------------------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------ANALYSIS BEGINS HERE---------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  
    for(int iEvent = 0; iEvent < nEvents; ++iEvent) {//Event Loop BEGIN
        
        if (!pythia.next()) continue;

        //cout << "Processing event " << iEvent << "." << endl;
        numD0ToKPi = 0;
        D0Idx = 0;
        numSoftPion = 0;
        hasDtoKPi = false;
        dStarIsMother = false;
        
        for(int i = 0; i < pythia.event.size(); ++i){  //doing proper mother/daughter checks
        
            trkPID = pythia.event[i].id();
            ptBin = getPtBin(pythia.event[i].pT());
            if(ptBin < 0) { continue; }
           
            if(TMath::Abs(trkPID) == 421 && TMath::Abs(pythia.event[i].eta()) < 1.0){  //event has a D0
               
                daughter1Idx = pythia.event[i].daughter1();                              //since only one d0/event that decays this way, daughterIdx still useful
                daughter2Idx = pythia.event[i].daughter2();
                daughterPID1 = TMath::Abs(pythia.event[daughter1Idx].id());
                daughterPID2 = TMath::Abs(pythia.event[daughter2Idx].id());
               
                //if((TMath::Abs(pythia.event[daughter1Idx].eta()) > 1.0) || (TMath::Abs(pythia.event[daughter2Idx].eta()) > 1.0)) {continue;}  //makes sure daughters are in acceptance
                //if((TMath::Abs(pythia.event[daughter1Idx].pT()) > .15) || (TMath::Abs(pythia.event[daughter2Idx].pT()) > .15)) {continue;}    //makes sure daughters are in acceptance
               
                if((daughterPID1 == 321 && daughterPID2 == 211) || (daughterPID1 == 211 && daughterPID2 == 321)){ hasDtoKPi = true; D0Idx = i; D0PtBin = ptBin; numD0ToKPi++; }
                
            }
        }//end D0->kPi check
       
        if(hasDtoKPi && numD0ToKPi > 1) { continue; }
        //if(hasDtoKPi) { cout << "***************Event: " << iEvent << " ****************" << endl << "This event has a D0->KPi candidate." << endl; }
       
        numFinalStateTracks = 0;
        numFinalStateTracksNoSoftPion = 0;
        
        if(hasDtoKPi){//If there is a D0->kPi present, save this to buffer
            
            motherIdx = pythia.event[D0Idx].mother1();
            totalD0Yield->Fill(D0PtBin+1);
            
            if(TMath::Abs(pythia.event[motherIdx].id()) == 413) { dStarIsMother = true; }
            
            if(dStarIsMother) { totalDStarToD0->Fill(D0PtBin+1); }
            
            D0Mass->Fill(pythia.event[D0Idx].m());                  //this just checks for the existence of a D0
            D0PtDist->Fill(pythia.event[D0Idx].pT()); 
            D0EtaDist->Fill(pythia.event[D0Idx].eta()); 
            D0PhiDist->Fill(pythia.event[D0Idx].phi()); 
            
            DMesons.push_back(pythia.event[D0Idx].p());                    //store D0 4-mom into sibling container
            DMesonsForMixing.push_back(pythia.event[D0Idx].p());           //store D0 4-mom into mixed container
            
            if(dStarIsMother) { DMesonsForMixingNoSoftPion.push_back(pythia.event[D0Idx].p()); } //stores d*->D0+Pi in a separate buffer
        }
        
        for(int i = 0; i < pythia.event.size(); ++i){ //Track loop to generate hadron lists BEGIN -- stores "sibling" tracks regardless to simplify coding -- marginal inefficiency
            
                motherIdx = pythia.event[i].mother1();
                if(TMath::Abs(pythia.event[i].eta()) > 1.0 || pythia.event[i].pT() < .15) { continue; }  //acceptance filter
                if( !(TMath::Abs(pythia.event[i].id()) == 211 || TMath::Abs(pythia.event[i].id()) == 321 || TMath::Abs(pythia.event[i].id()) == 2212) ) {continue;} //checks if hadron
                if(hasDtoKPi && (i == daughter1Idx || i == daughter2Idx)) { continue; }                     //this ensures not-storing a d0 decay daughter
                if(pythia.event[i].daughter1() != 0) { continue; } //throws out non-final state particls
               
                if(hasDtoKPi){
                
                    hadrons.push_back(pythia.event[i].p());   //stores sibling event hadrons
                    hadronPtDist->Fill(pythia.event[i].pT());
                    hadronEtaDist->Fill(pythia.event[i].eta());
                    hadronPhiDist->Fill(pythia.event[i].phi());
                }
                if(mixerBufferIdx < NUM_MIXED_EVENTS && !hasDtoKPi) { mixerHadrons[mixerBufferIdx][numFinalStateTracks] = pythia.event[i].p(); }
                numFinalStateTracks++;
												
                if(dStarIsMother && TMath::Abs(pythia.event[motherIdx].id()) != 413) { hadronsNoSoftPion.push_back(pythia.event[i].p()); }  //stores all but soft pion from D*->D0pi
                if(dStarIsMother && TMath::Abs(pythia.event[motherIdx].id()) == 413) { softPionOnly.push_back(pythia.event[i].p()); }  //stores ONLY soft pion from D*->D0pi
                //if(dStarIsMother && TMath::Abs(pythia.event[motherIdx].id()) == 413) { numSoftPion++;}
        }//Track loop to generate hadron lists END
            
       
            
        if( mixerBufferIdx < NUM_MIXED_EVENTS && !hasDtoKPi) { 
        
            numMixerHadrons[mixerBufferIdx] = numFinalStateTracks; 
            mixerBufferIdx++; 
            if(DEBUG) { cout << "mixerBufferIdx: " << mixerBufferIdx << endl; }
    
        } //this tells us how many tracks per mixed event

		if(hasDtoKPi){//Sibling correlation loop begins here.
	        for(int dIndex = 0; dIndex < DMesons.size(); dIndex++){  //with soft pion
                for(int hIndex = 0; hIndex < hadrons.size(); hIndex++){

					delEta = DMesons[dIndex].eta() - hadrons[hIndex].eta();
                    delPhi = DMesons[dIndex].phi() - hadrons[hIndex].phi();
                    //delEta = TMath::Abs(delEta);
                
                    if(delPhi > TMath::TwoPi()) { delPhi = delPhi - TMath::TwoPi(); }
                    else if(delPhi < -TMath::TwoPi()) { delPhi = delPhi + TMath::TwoPi(); }

                    //delPhi = TMath::Abs(delPhi);
                    //delPhiCp = -delPhi;
              
                    if(delPhi > (3*TMath::PiOver2()+phiBinShift)) { delPhi = delPhi - TMath::TwoPi(); }
                    //if(delPhiCp < (-TMath::PiOver2()+phiBinShift)) { delPhiCp = delPhiCp + TMath::TwoPi(); }
                    if(delPhi < (-TMath::PiOver2()+phiBinShift)) { delPhi = delPhi + TMath::TwoPi(); }

                    sibDHadronCorr[D0PtBin]->Fill(delEta, delPhi);
                    //sibDHadronCorr[D0PtBin]->Fill(-delEta, delPhi);
                    //sibDHadronCorr[D0PtBin]->Fill(delEta, delPhiCp);
                    //sibDHadronCorr[D0PtBin]->Fill(-delEta, delPhiCp);
                    //if(!dStarIsMother) {

                      //  sibDHadronCorrNoSoftPion[D0PtBin]->Fill(delEta, delPhi);
                       // sibDHadronCorrNoSoftPion[D0PtBin]->Fill(-delEta, delPhi);
                        //sibDHadronCorrNoSoftPion[D0PtBin]->Fill(delEta, delPhiCp);
                        //sibDHadronCorrNoSoftPion[D0PtBin]->Fill(-delEta, delPhi);
                    //}

                }
            } //with soft pion
            if(dStarIsMother){
            
                for(int dIndex = 0; dIndex < DMesons.size(); dIndex++){ //SOFT PION ONLY
                    for(int hIndex = 0; hIndex < softPionOnly.size(); hIndex++){
                        
                        delEta = DMesons[dIndex].eta() - softPionOnly[hIndex].eta();
                        delPhi = DMesons[dIndex].phi() - softPionOnly[hIndex].phi();
                       
                        if(delPhi > TMath::TwoPi()) { delPhi = delPhi - TMath::TwoPi(); }
                        else if(delPhi < -TMath::TwoPi()) { delPhi = delPhi + TMath::TwoPi(); }

                        if(delPhi > (3*TMath::PiOver2()+phiBinShift)) { delPhi = delPhi - TMath::TwoPi(); }
                        if(delPhi < (-TMath::PiOver2()+phiBinShift)) { delPhi = delPhi + TMath::TwoPi(); }

                        sibDHadronCorrD0SoftPionOnly[D0PtBin]->Fill(delEta, delPhi);
			                        

                    }
                }
            
                for(int dIndex = 0; dIndex < DMesons.size(); dIndex++){ //no soft pion
                    for(int hIndex = 0; hIndex < hadronsNoSoftPion.size(); hIndex++){

                        delEta = DMesons[dIndex].eta() - hadronsNoSoftPion[hIndex].eta();
                        delPhi = DMesons[dIndex].phi() - hadronsNoSoftPion[hIndex].phi();
                        //delEta = TMath::Abs(delEta);
                    
                        if(delPhi > TMath::TwoPi()) { delPhi = delPhi - TMath::TwoPi(); }
                        else if(delPhi < -TMath::TwoPi()) { delPhi = delPhi + TMath::TwoPi(); }

                        //delPhi = TMath::Abs(delPhi);
                        //delPhiCp = -delPhi;
                  
                        if(delPhi > (3*TMath::PiOver2()+phiBinShift)) { delPhi = delPhi - TMath::TwoPi(); }
                        //if(delPhiCp < (-TMath::PiOver2()+phiBinShift)) { delPhiCp = delPhiCp + TMath::TwoPi(); }
                        if(delPhi < (-TMath::PiOver2()+phiBinShift)) { delPhi = delPhi + TMath::TwoPi(); }

                        sibDHadronCorrNoSoftPion[D0PtBin]->Fill(delEta, delPhi);
                        //sibDHadronCorrNoSoftPion[D0PtBin]->Fill(-delEta, delPhi);
                        //sibDHadronCorrNoSoftPion[D0PtBin]->Fill(delEta, delPhiCp);
                        //sibDHadronCorrNoSoftPion[D0PtBin]->Fill(-delEta, delPhiCp);
                    }
                }
            }
        }  //Sibling correlation loop ends here.        

		if(mixerBufferIdx == NUM_MIXED_EVENTS && DMesonsForMixing.size() > 0 ){
		
            //cout << "****************Event: " << iEvent << " *******************" << endl;
            //cout << "Event Mixing Begin....." << endl;
        
            for(int dIndex = 0; dIndex < DMesonsForMixing.size(); dIndex++){
                ptBin = getPtBin(DMesonsForMixing[dIndex].pT());
                if(ptBin < 0) {continue;}
        
                for(int bufferIdx = 0; bufferIdx < 5; bufferIdx++){
                    for(int hIndex = 0; hIndex < numMixerHadrons[bufferIdx]; hIndex++){
					
                        if(DEBUG) { cout << "hIndex: " << hIndex <<endl; }
                        
                        delEta = DMesonsForMixing[dIndex].eta() - mixerHadrons[bufferIdx][hIndex].eta();
                        delPhi = DMesonsForMixing[dIndex].phi() - mixerHadrons[bufferIdx][hIndex].phi();
                        //delEta = TMath::Abs(delEta);
                        
                        if(delPhi > TMath::TwoPi()) { delPhi = delPhi - TMath::TwoPi(); }
                        else if(delPhi < -TMath::TwoPi()) { delPhi = delPhi + TMath::TwoPi(); }

                        //delPhi = TMath::Abs(delPhi);
                        //delPhiCp = -delPhi;
                      
                        if(delPhi > (3*TMath::PiOver2()+phiBinShift)) { delPhi = delPhi - TMath::TwoPi(); }						
                        //if(delPhiCp < (-TMath::PiOver2()+phiBinShift)) { delPhiCp = delPhiCp + TMath::TwoPi(); }
                        if(delPhi < (-TMath::PiOver2()+phiBinShift)) { delPhi = delPhi + TMath::TwoPi(); }

                        mixDHadronCorr[ptBin]->Fill(delEta, delPhi);
                        //mixDHadronCorr[ptBin]->Fill(-delEta, delPhi);
                       // mixDHadronCorr[ptBin]->Fill(delEta, delPhiCp);
                        //mixDHadronCorr[ptBin]->Fill(-delEta, delPhiCp);
						
                    }
                }
			}
                
            mixerBufferIdx = 0;
            
            DMesonsForMixing.clear(); //clears d-meson mixing buffer 
            
		}

        
        
        piPos.clear();
        kNeg.clear();
        piPosMotherIdx.clear();
        kNegMotherIdx.clear();
        
        //clear same event containers
        hadrons.clear();
        hadronsNoSoftPion.clear();
        DMesons.clear();
        softPionOnly.clear();
  
    
        trackMult->Fill(numFinalStateTracks);
      // hasDMeson = false;
    } //Event Loop End

    double NSib = 1;
    double NMix = 1;
    double ratio = 1;
  
    for(int i = 0; i < 5; i++){  //calculate correlation ratios
  
        NSib = sibDHadronCorr[i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
        NMix = mixDHadronCorr[i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
        ratio = NSib/NMix;
        ratioDHadronCorr[i]->Divide(sibDHadronCorr[i], mixDHadronCorr[i], 1, ratio);
	  
	    NSib = sibDHadronCorrNoSoftPion[i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
	    NMix = mixDHadronCorr[i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
	    ratio = NSib/NMix;
	    ratioDHadronCorrNoSoftPion[i]->Divide(sibDHadronCorrNoSoftPion[i], mixDHadronCorr[i], 1, ratio);
    }
  
  
    // Statistics on event generation.
    pythia.stat();

    sibDHadronCorrPhi[0] = (TH1D*) sibDHadronCorr[0]->ProjectionY();
    sibDHadronCorrPhi[1] = (TH1D*) sibDHadronCorr[1]->ProjectionY();
    sibDHadronCorrPhi[2] = (TH1D*) sibDHadronCorr[2]->ProjectionY();
    sibDHadronCorrPhi[3] = (TH1D*) sibDHadronCorr[3]->ProjectionY();
    sibDHadronCorrPhi[4] = (TH1D*) sibDHadronCorr[4]->ProjectionY();
    sibDHadronCorrNoSoftPionPhi[0] = (TH1D*) sibDHadronCorrNoSoftPion[0]->ProjectionY();
    sibDHadronCorrNoSoftPionPhi[1] = (TH1D*) sibDHadronCorrNoSoftPion[1]->ProjectionY();
    sibDHadronCorrNoSoftPionPhi[2] = (TH1D*) sibDHadronCorrNoSoftPion[2]->ProjectionY();
    sibDHadronCorrNoSoftPionPhi[3] = (TH1D*) sibDHadronCorrNoSoftPion[3]->ProjectionY();
    sibDHadronCorrNoSoftPionPhi[4] = (TH1D*) sibDHadronCorrNoSoftPion[4]->ProjectionY();
    
    // Write tree.
    //T->Print();
    //T->Write();
    trackMult->Write();
    D0Mass->Write();
    D0PtDist->Write(); 
    D0EtaDist->Write();
    D0PhiDist->Write();
    totalD0Yield->Write();
    
    hadronPtDist->Write();
    hadronEtaDist->Write();
    hadronPhiDist->Write();
    sibDHadronCorr[0]->Write();
    sibDHadronCorr[1]->Write();
    sibDHadronCorr[2]->Write();
    sibDHadronCorr[3]->Write();
    sibDHadronCorr[4]->Write();
    sibDHadronCorrPhi[0]->Write();
    sibDHadronCorrPhi[1]->Write();
    sibDHadronCorrPhi[2]->Write();
    sibDHadronCorrPhi[3]->Write();
    sibDHadronCorrPhi[4]->Write();
    mixDHadronCorr[0]->Write();
    mixDHadronCorr[1]->Write();
    mixDHadronCorr[2]->Write();
    mixDHadronCorr[3]->Write();
    mixDHadronCorr[4]->Write();
    ratioDHadronCorr[0]->Write();
    ratioDHadronCorr[1]->Write();
    ratioDHadronCorr[2]->Write();
    ratioDHadronCorr[3]->Write();
    ratioDHadronCorr[4]->Write();

    
    sibDHadronCorrD0SoftPionOnly[0]->Write();
    sibDHadronCorrD0SoftPionOnly[1]->Write();
    sibDHadronCorrD0SoftPionOnly[2]->Write();
    sibDHadronCorrD0SoftPionOnly[3]->Write();
    sibDHadronCorrD0SoftPionOnly[4]->Write();
    sibDHadronCorrNoSoftPion[0]->Write();
    sibDHadronCorrNoSoftPion[1]->Write();
    sibDHadronCorrNoSoftPion[2]->Write();
    sibDHadronCorrNoSoftPion[3]->Write();
    sibDHadronCorrNoSoftPion[4]->Write();
    sibDHadronCorrNoSoftPionPhi[0]->Write();
    sibDHadronCorrNoSoftPionPhi[1]->Write();
    sibDHadronCorrNoSoftPionPhi[2]->Write();
    sibDHadronCorrNoSoftPionPhi[3]->Write();
    sibDHadronCorrNoSoftPionPhi[4]->Write();
    //mixDHadronCorrNoSoftPion[0]->Write();
    //mixDHadronCorrNoSoftPion[1]->Write();
    //mixDHadronCorrNoSoftPion[2]->Write();
    //mixDHadronCorrNoSoftPion[3]->Write();
    //mixDHadronCorrNoSoftPion[4]->Write();
    ratioDHadronCorrNoSoftPion[0]->Write();
    ratioDHadronCorrNoSoftPion[1]->Write();
    ratioDHadronCorrNoSoftPion[2]->Write();
    ratioDHadronCorrNoSoftPion[3]->Write();
    ratioDHadronCorrNoSoftPion[4]->Write();
    
    CAndCBarCorr->Write();
    DAndDBarCorr->Write();
    
    delete file;

    return 0;
}

double calcInvMass(double mass1, double mass2, double p1x, double p1y, double p1z, double p2x, double p2y, double p2z){

    return TMath::Sqrt((mass1*mass1)+(mass2*mass2)+ 2*TMath::Sqrt((mass1*mass1)+(p1x*p1x+p1y*p1y+p1z*p1z))*TMath::Sqrt((mass2*mass2)+(p2x*p2x+p2y*p2y+p2z*p2z))-2*(p1x*p2x+p1y*p2y+p1z*p2z));
    
}

int getPtBin(double pt){

    if(pt >= 0.0 && pt < 1.0)  { return 0; }
    if(pt >= 1.0 && pt < 2.0)  { return 1; }
    if(pt >= 2.0 && pt < 3.0)  { return 2; }
    if(pt >= 3.0 && pt < 4.0)  { return 3; }
    if(pt >= 4.0 && pt < 10.0) { return 4; }

    else return -1;
}
    
