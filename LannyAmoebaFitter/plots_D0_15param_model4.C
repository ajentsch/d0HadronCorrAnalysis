void plots_D0_15param_model4(){

  TString str = "Amoeba_D0hadron_Cent0_pt2-3_model4.hist";
  ifstream in;
  in.open(str);
  cout << "Reading " << str << endl;

  Float_t chi2,nit,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
  // TNtuple(name,title, var list)
  TNtuple *nt  = new TNtuple("nt","NTuple","chi2:p1:p2:p3:p4:p5:p6:p7:p8:p9:p10:p11:p12:p13:p14");
  TNtuple *nt2 = new TNtuple("nt2","NTuple","chi2:p15");

  //cout.precision(4);   // print numbers to 4 sig figs
  cout.width(4);   // print numbers to 4 sig figs

  // Read text file
  while (1) {               // loop until done
  in >> chi2 >> nit >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> p9 >>p10 >>p11 >>p12 >>p13 >> p14 >> p15;
    if (!in.good()) break;  // stop if EOF or error
    nt->Fill(chi2,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14);
    nt2->Fill(chi2,p15);
  }
  in.close();

  gStyle->SetPalette(1); // use smooth color palette
  c2 = new TCanvas("c2","NTuple plots",1400,1500);
  c2->Divide(4,4);
  c2->cd(1);
  nt->Draw("chi2:p1","(chi2<5.0)");
  c2->cd(2);
  nt->Draw("chi2:p2","(chi2<5.0)");
  c2->cd(3);
  nt->Draw("chi2:p3","(chi2<5.0)");
  c2->cd(4);
  nt->Draw("chi2:p4","(chi2<5.0)");
  c2->cd(5);
  nt->Draw("chi2:p5","(chi2<5.0)");
  c2->cd(6);
  nt->Draw("chi2:p6","(chi2<5.0)");
  c2->cd(7);
  nt->Draw("chi2:p7","(chi2<5.0)");
  c2->cd(8);
  nt->Draw("chi2:p8","(chi2<5.0)");
  c2->cd(9);
  nt->Draw("chi2:p9","(chi2<5.0)");
  c2->cd(10);
  nt->Draw("chi2:p10","(chi2<5.0)");
  c2->cd(11);
  nt->Draw("chi2:p11","(chi2<5.0)");
  c2->cd(12);
  nt->Draw("chi2:p12","(chi2<5.0)");
  c2->cd(13);
  nt->Draw("chi2:p13","(chi2<5.0)");
  c2->cd(14);
  nt->Draw("chi2:p14","(chi2<5.0)");
  c2->cd(15);
  nt2->Draw("chi2:p15","(chi2<5.0)");


}
