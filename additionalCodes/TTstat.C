#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "../include/const.hpp"
#include "../include/kineticAnalyze.hpp"
#include "../src/kineticAnalyze.cpp"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;


void createHistos(int exp, double cmunuL, double cmunuR, double cmunu, double dmunu, double bkgd, double ttbar){

//___________________________________________ SME preparation ___________________________________________//

    int time = 24;              //hours
    int bin = 24;               //bins

    KineticAnalyze k;

    TH1F* hBkgd = k.statHistosConst("hBkgd", value);
    TH1F* hTTbarSM = statHistosConst("hTTbarSM", value);
    TH1F* hAsimovNullHyp = new statHistosConst("hAsimovNullHyp", value);
    TH1F* hSME = statHistosf("hSME", int exp, double cmunuL, double cmunuR, double cmunu, double dmunu, bool isXX);


   TFile* fOutput = new TFile("../stats/statgTT.root","RECREATE");
   hBkgd->Write();
   hTTbarSM->Write();
   hAsimovNullHyp->Write();
   hSME->Write();
   fOutput->Write();
   fOutput->Close();
}


void TTstat(){

//_____________________ All needed parameters _____________________//

/*
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

*/

   return;

}
