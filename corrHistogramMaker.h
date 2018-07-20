//-------------------------------------------Functions------------------------------------------------------------------//    

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


    
	
