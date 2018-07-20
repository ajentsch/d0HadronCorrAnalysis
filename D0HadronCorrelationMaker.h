//-------------------------------------------Functions------------------------------------------------------------------//    

void saveCutsToFile(TFile * file, int numFiles, TString path){

    double numEvents = 0;
        
    numEvents = (((TH1I*) file->Get("number_of_events_used"))->GetBinContent(2))/(1000000);
    
    ofstream cutFile;
    TString cutFileName = "Cuts_and_data_information.txt";
    TString cutOutputFile = path + cutFileName;
    cutFile.open (cutOutputFile);
    
    cutFile << "Important information about this data run" << endl << endl;
    cutFile << "Number of Events: " << numEvents << "M Events" << endl << endl;
    
    numFiles = ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(49);
    
    ((TH1D*) file->Get("HistOfCuts"))->Scale(1/numFiles);
    
    cutFile << "Trigger D0 Cuts" << endl << endl;
    cutFile << "D0 InvMass Signal Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(1) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(2) << endl;
    cutFile << "D0 SideBandLeft   Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(3) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(4) << endl;
    cutFile << "D0 SideBandRight  Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(5) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(6) << endl;
    cutFile << "D0 Pt Cuts             -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(7) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(8) << endl;
    cutFile << "D0 DecayLengthMax(no upper bound)      -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(11) << endl;
    
    cutFile << "D0 DaughterPtCut       -\t" << "Pi : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(9) << "\t" << "K  : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(10) << endl;
    
    cutFile << endl << "D0 pt 0-1 GeV/c cuts" << endl;
    
    cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(12) << endl;
    cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(13) << endl;
    cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(14) << endl;
    cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(15) << endl;
    cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(16) << endl;
    
    cutFile << endl << "D0 pt 1-2 GeV/c cuts" << endl;
    
    cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(17) << endl;
    cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(18) << endl;
    cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(19) << endl;
    cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(20) << endl;
    cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(21) << endl;
    
    cutFile << endl << "D0 pt 2-5 GeV/c cuts" << endl;
    
    cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(22) << endl;
    cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(23) << endl;
    cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(24) << endl;
    cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(25) << endl;
    cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(26) << endl;
    
    cutFile << endl << "D0 pt 5-10 GeV/c cuts" << endl;
    
    cutFile << "D0 DecayLengthMinimum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(27) << endl;
    cutFile << "D0 DaughterDCAMaximum  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(28) << endl;
    cutFile << "D0 DCA to PV Maximum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(29) << endl;
    cutFile << "Daughter Pion DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(30) << endl;
    cutFile << "Daughter Kaon DCA to PV minimum  -\t"  << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(31) << endl;
    
    cutFile << endl << "------------------------------------------" << endl << "Associated hadron Cuts" << endl << endl;
    cutFile << "associated hadron Pt   -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(37) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(38) << endl;
    cutFile << "Num of events used to mix  -\t" << ": " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(39) << endl;
    cutFile << "HFT Tracks only??  -\t" << "YES(1), NO(0): " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(42) << endl;
    cutFile << "Track Chi2 Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(43) << endl;
    cutFile << "Track DCA Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(44) << endl;
    
    cutFile.close();

    return;
}    

















void formatCorrHist(TH2D* hist) {

    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    //hist->GetXaxis()->SetFontSize();
    //hist->GetYaxis()->CenterTitle();
    
    hist->GetZaxis()->SetTitle("#frac{#Delta#rho}{#rho_{ref}}");
    
    hist->GetXaxis()->SetTitleSize(.065);
    hist->GetYaxis()->SetTitleSize(.065); 
    hist->GetZaxis()->SetTitleSize(.03);
    
    hist->GetXaxis()->SetTitleOffset(1.5);
    hist->GetYaxis()->SetTitleOffset(1.5); 
    hist->GetZaxis()->SetTitleOffset(1.62);
   
    return;
    
}

void formatCorrHist(TH1D* hist) {

    hist->GetXaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
   
   return;
    
}


void formatCorrHist(TH2D* hist, TString title) {

    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3); 
    hist->GetXaxis()->SetTitle("#Delta#eta");
    hist->GetYaxis()->SetTitle("#Delta#phi");
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3); 
    hist->SetNameTitle(title, title);
    
    return;
    
}


    
	
