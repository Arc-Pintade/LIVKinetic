#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "../include/const.hpp"
#include "../include/kineticAnalyze.hpp"

using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;


void createHistos(int exp, double cmunuL, double cmunuR, double cmunu, double dmunu, double isXX, double bkgd, double ttbar){

//___________________________________________ SME preparation ___________________________________________//

    int time = 24;              //hours
    int bin = 24;               //bins

    KineticAnalyze k;

    TH1F* hBkgd = k.statHistosConst("hBkgd", bkgd);
    TH1F* hTTbarSM = k.statHistosConst("hTTbarSM", ttbar);
    TH1F* hAsimovNullHyp = k.statHistosConst("hAsimovNullHyp", bkgd+ttbar);
    TH1F* hSME = k.statHistosf("hSME", exp, cmunuL, cmunuR, cmunu, dmunu, isXX, bkgd, ttbar);


   TFile* fOutput = new TFile("../stats/statsgTT.root","RECREATE");
   hBkgd->Write();
   hTTbarSM->Write();
   hAsimovNullHyp->Write();
   hSME->Write();
   fOutput->Write();
   fOutput->Close();
}


void TTstat(){

//_____________________ All needed parameters _____________________//

double c[3]; c[0]=0.1; c[1]=0.01; c[2]=0.001;

    for(int i=0; i<1; i++){

       createHistos(3, c[i],0,0,0,true,2,2);
/*
       RooStats::HistFactory::Measurement meas("gTT", "stats");

       meas.SetOutputFilePrefix("../stats/statsgTT");
       meas.SetExportOnly(false);

       //meas.SetPOI("SigXsecOverSM");
       meas.AddPOI("muSME");

       meas.SetLumi(1);
       meas.SetLumiRelErr(0.02);

       RooStats::HistFactory::Channel chan("channel");

       chan.SetData("hAsimovNullHyp", "../stats/statsgTT.root");

       RooStats::HistFactory::Sample ttbar("hTTbarSM", "hTTbarSM", "../stats/statsgTT.root");
       ttbar.AddNormFactor("SigXsecOverSM", 1, 0, 3);
       ttbar.AddOverallSys("syst1",  0.95, 1.05);
       chan.AddSample(ttbar); 

       RooStats::HistFactory::Sample bkgd("hBkgd", "hBkgd", "../stats/statsgTT.root");
       bkgd.AddOverallSys("syst2", 0.7, 1.3 );
       chan.AddSample(bkgd);        

       RooStats::HistFactory::Sample SME("hSME", "hSME", "../stats/statsgTT.root");
       SME.AddNormFactor("muSME", 1, -10, 10);
       chan.AddSample(SME); 
     

       meas.AddChannel(chan);
       meas.CollectHistograms();
       meas.PrintTree();
       MakeModelAndMeasurementFast(meas);
*/
    }


   return;

}
