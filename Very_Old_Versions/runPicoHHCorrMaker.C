/* **************************************************
 *  A macro to run StPicoHHCorrMaker
 *
 *  Authors:  Alex Jentsch
 *
 *  
 *
 * **************************************************
 */



void runPicoHHCorrMaker(TString d0list, TString outFileName, TString badRunListFileName = "/global/homes/a/ajentsch/myAnalysis/picoList_bad_MB.list")
//void runPicoD0AnaMaker(TString d0list = "test200k.list", TString outFileName = "test.root", TString badRunListFileName = "/global/homes/a/ajentsch/myAnalysis
//picoList_bad_MB.list")
{
   //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
   string SL_version = "SL16d";
   string env_SL = getenv("STAR");
   if (env_SL.find(SL_version) == string::npos)
   {
      cout << "Environment Star Library does not match the requested library in runPicoD0EventMaker.C. Exiting..." << endl;
      exit(1);
   }

   gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
   loadSharedLibraries();

   gSystem->Load("StBTofUtil");
   gSystem->Load("StPicoDstMaker");
   gSystem->Load("StPicoPrescales");
   gSystem->Load("StPicoCutsBase");
   gSystem->Load("StPicoCutsBase");
   gSystem->Load("StPicoD0EventMaker");
   gSystem->Load("StPicoHHCorrMaker");
   gSystem->Load("StPicoHFMaker");
   gSystem->Load("StRefMultCorr");
   gSystem->Load("StMixedEventBuffer");   

   chain = new StChain();

   // create list of picoDst files
   TString command = "sed 's/hft/picodsts/g' " + d0list + " >correspondingPico.list";
   gSystem->Exec(command.Data());
   command = "sed -i 's/picoD0/picoDst/g' correspondingPico.list";
   gSystem->Exec(command.Data());
   command = "sed -i 's/Pico16a/physics2/g' correspondingPico.list";
   gSystem->Exec(command.Data());
   command = "sed -i 's/D0//g' correspondingPico.list";
   gSystem->Exec(command.Data());
   
   
   StPicoDstMaker* picoDstMaker = new StPicoDstMaker(0, "correspondingPico.list", "picoDstMaker");
   
   //StPicoD0AnaMaker*  picoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", d0list, outFileName.Data(), picoDstMaker);
   StRefMultCorr* grefmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr();
   cout<<"here"<<endl;
   grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
   grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");
   for(Int_t i=0;i<6;i++){
      cout << i << " " << grefmultCorrUtil->get(i, 0) << endl;
   }
   
   //StPicoMixedEventMaker* mixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, grefmultCorrUtil, "Buffer Name"); 
   
   StPicoHHCorrMaker*  picoHHCorrMaker = new StPicoHHCorrMaker("picoHHCorrMaker", d0list, outFileName.Data(), picoDstMaker);
   
   //Need to check the input here -- may need to change the constructor to not have the other inputs   
   
   
                                                                
   StHFCuts* hadronCuts = new StHFCuts("hadronCuts");
   picoHHCorrMaker->setHFCuts(hadronCuts);

   // -------------- USER variables -------------------------

   // -- File name of bad run list
   hadronCuts->setBadRunListFileName(badRunListFileName);

   // add your cuts here.
   //
   hadronCuts->setCutRequireHFT(true); //BE SURE TO CHANGE IN ANA MAKER IF YOU CHANGE HERE FOR CUT FILE! 
   hadronCuts->setCutVzMax(6.0);
   hadronCuts->setCutVzVpdVzMax(3.0);
   hadronCuts->setCutPionPtRange(0.15, 20.0);
   hadronCuts->setCutKaonPtRange(0.15, 20.0);
   // tracking
   hadronCuts->setCutNHitsFitMin(20);
   hadronCuts->setCutNHitsFitnHitsMax(0.52);
   hadronCuts->setCutPrimaryDCAtoVtxMax(3.0);
   // pions
   hadronCuts->setCutTPCNSigmaPion(3.0);

   // kaons
   hadronCuts->setCutTPCNSigmaKaon(2.0);

   
   
   chain->Init();

   
   int per = 0;
   int nEntries = picoHHCorrMaker->getEntries();

   //int nEntries = 50000;
   for (int iEvent = 0; iEvent < nEntries; ++iEvent)
   {
      
    if(nEntries>=1000){
   
    	per = iEvent/(nEntries/100);
     
        if (iEvent%(nEntries/10)==0 && nEntries>99){cout<<"working on event: "<<iEvent<<" ("<<per<<"%)"<<endl;}   
    }

      chain->Clear();
      int iret = chain->Make();
      if (iret)
      {
         cout << "Bad return code!" << iret << endl;
         break;
      }
   }

   chain->Finish();
   delete chain;

   // delete list of picos
   command = "rm -f correspondingPico.list";
   gSystem->Exec(command.Data());

}
