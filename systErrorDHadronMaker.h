TH2D* errorsFromSecondaries(TH2D* hist){

    double errorScaleFactor = .015;
    
    hist->Scale(errorScaleFactor);

    return hist;

}


TH2D* errorsFromBFeedDown(TH2D* hist){

    double errorScaleFactor = .02;
    
    hist->Scale(errorScaleFactor);

    return hist;

}

TH2D* errorsFromPileup(TH2D* hist){

    return hist;
    
}


TH2D* errorsFromFBarFactor(TH2D* histNominal, TH2D* histOption1, TH2D* histOption2, int numEta, int numPhi){

    double nominalBin;
    double option1Bin;
    double option2Bin;
    double errorResult;
    
    TH2D *errorResultHist = new TH2D("","", numEta, -2, 2, numPhi, -TMath::PiOver2(), 3*TMath::PiOver2());

    for(int i = 1; i < numEta+1; i++){
        for(int j = 1; j < numPhi+1; j++){
            
            nominalBin = histNominal->GetBinContent(i, j);
            option1Bin = histOption1->GetBinContent(i, j);
            option2Bin = histNominal->GetBinContent(i, j);
            
            cout << "Bin : " << i << ", " << j << "  : " << (nominalBin-option1Bin) << endl;
            
            errorResult = .5*(TMath::Abs(nominalBin-option1Bin) + TMath::Abs(nominalBin-option2Bin));
            
            errorResultHist->SetBinContent(i, j, errorResult);
        }
    }
    
    return errorResultHist;
}

TH2D* errorsFromTopologyCut(TH2D* histNominal, TH2D* histOption1, TH2D* histOption2, int numEta, int numPhi){

    TH2D *errorResultHist = new TH2D("","", numEta, -2, 2, numPhi, -TMath::PiOver2(), 3*TMath::PiOver2());
   
    errorResultHist = (TH2D*) errorsFromFBarFactor(histNominal, histOption1, histOption2, numEta, numPhi);
    
    return errorResultHist;
}



TH2D* errorsFromBackgroundSubtraction(TH2D* histNominal, TH2D* histLS, int numEta, int numPhi){

    TH2D *errorResultHist = new TH2D("","", numEta, -2, 2, numPhi, -TMath::PiOver2(), 3*TMath::PiOver2());
    
    double nominalBin;
    double LSBin;
    double errorResult;
    
    for(int i = 1; i < numEta+1; i++){
        for(int j = 1; j < numPhi+1; j++){
            
            nominalBin = histNominal->GetBinContent(i, j);
            LSBin = histLS->GetBinContent(i, j);
            
            errorResult = .5*(TMath::Abs(nominalBin-LSBin));
            
            errorResultHist->SetBinContent(i, j, errorResult);
        }
    }

    return errorResultHist;
}



TH2D* calculateFinalSystematicErrors(TH2D* h1, TH2D* h2, TH2D* h3, TH2D* h4, TH2D* h5, int numEta, int numPhi){

    TH2D *result = new TH2D("", "", numEta, -2, 2, numPhi, -TMath::PiOver2(), 3*TMath::PiOver2());

    double e1;
    double e2;
    double e3;
    double e4;
    double e5;
    double final;
    
    for(int i = 1; i < numEta+1; i++){
        for(int j = 1; j < numPhi+1; j++){
        
            e1 = h1->GetBinContent(i, j);
            e2 = h2->GetBinContent(i, j);
            e3 = h3->GetBinContent(i, j);
            e4 = h4->GetBinContent(i, j);
            e5 = h5->GetBinContent(i, j);
            final = TMath::Sqrt((e1*e1)+(e2*e2)+(e3*e3)+(e4*e4)+(e5*e5));
            
            result->SetBinContent(i, j, final);
        }
    }



    return result;

}