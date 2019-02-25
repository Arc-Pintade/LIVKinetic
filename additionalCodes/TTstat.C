#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "include/const.hpp"
#include "include/NiqueTaMereRoot.h"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;


void createHistos(TString lum, TString cross, TString exp, TString XX){

//___________________________________________ SME preparation ___________________________________________//

    double cmunuL = 0.1;
    double cmunuR = 0;
    int time = 24;              //hours
    int bin = 24;               //bins

    TH1F* hBkgd = new TH1F("hBkgd", "hBkgd", 24, 0, 24);
    TH1F* hTTbarSM = new TH1F("hTTbarSM", "hTTbarSM", 24, 0, 24);
    TH1F* hAsimovNullHyp = new TH1F("hTTbarSM", "hTTbarSM", 24, 0, 24);
    TH1F* hSME = new TH1F("hSME", "hSME", 24, 0, 24);


   TFile* fOutput = new TFile("stats/statHistFactory.root","RECREATE");
   hBkgd->Write();
   hTTbarSM->Write();
   hAsimovNullHyp->Write();
   hSME->Write();
   fOutput->Write();
   fOutput->Close();

   return;
}


void TTstat(){

//_____________________ All needed parameters _____________________//


  createHistos();

   RooStats::HistFactory::Measurement meas("myMeasurement", "stat");

   meas.SetOutputFilePrefix("stats/myMeasurement");
   meas.SetExportOnly(false);

   //meas.SetPOI("SigXsecOverSM");
   meas.AddPOI("muSME");

   meas.SetLumi(1);
   meas.SetLumiRelErr(0.02);

   RooStats::HistFactory::Channel chan("channel");

   chan.SetData("hAsimovNullHyp", "stats/statHistFactory.root");

   RooStats::HistFactory::Sample ttbar("hTTbarSM", "hTTbarSM", "stats/statHistFactory.root");
   ttbar.AddNormFactor("SigXsecOverSM", 1, 0, 3);
   ttbar.AddOverallSys("syst1",  0.95, 1.05);
   chan.AddSample(ttbar); 

   RooStats::HistFactory::Sample bkgd("hBkgd", "hBkgd", "stats/statHistFactory.root");
   bkgd.AddOverallSys("syst2", 0.7, 1.3 );
   chan.AddSample(bkgd);        

   RooStats::HistFactory::Sample SME("hSME", "hSME", "stats/statHistFactory.root");
   SME.AddNormFactor("muSME", 1, -10, 10);
   chan.AddSample(SME); 
 

   meas.AddChannel(chan);
   meas.CollectHistograms();
   meas.PrintTree();
   MakeModelAndMeasurementFast(meas);

    cout<<

   return;

}
