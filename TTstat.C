#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "include/kineticAnalyze.hpp"
#include "src/kineticAnalyze.cpp"


using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;


void createHistos(int exp, double cmunuGen, TString wilson, bool isZT, double bkgd, double ttbar){

//___________________________________________ SME preparation ___________________________________________//

    int time = 24;              //hours
    int bin = 24;               //bins

    KineticAnalyze k(cmunuGen, cmunuGen, cmunuGen, cmunuGen);

    TH1F* hBkgd = k.statHistosConst("hBkgd", bkgd);
    TH1F* hTTbarSM = k.statHistosConst("hTTbarSM", ttbar);
    TH1F* hAsimovNullHyp = k.statHistosConst("hAsimovNullHyp", bkgd+ttbar);
    TH1F* hSME = k.statHistosf("hSME", exp, wilson, isZT, bkgd, ttbar);
//    TH1F* hSME = k.statHistosConst("hSME", bkgd+ttbar);

   TFile* fOutput = new TFile("stats/"+TString::Format("%f",cmunuGen)+"statsgTT.root","RECREATE");
   hBkgd->Write();
   hTTbarSM->Write();
   hAsimovNullHyp->Write();
   hSME->Write();
   fOutput->Write();
   fOutput->Close();
}


void TTstat(){

//_____________________ All needed parameters _____________________//

    double ratioSB = 15.3576;
    double efficiency = 5.533;

    double bkgdVal[4], sigVal[4];
           bkgdVal[3] = 44045;      bkgdVal[2] = (970*efficiency*3000/ratioSB);     bkgdVal[1] = (3727*efficiency*15000/ratioSB);      bkgdVal[0] = (34810*efficiency*15000/ratioSB);
           sigVal[3] = 676431;      sigVal[2] = 970*efficiency*3000;     sigVal[1] = 3727*efficiency*15000;      sigVal[0] = 34810*efficiency*15000;

    double test[7];
            test[0] = 0.1;
            test[1] = 0.01;
            test[2] = 0.001;
            test[3] = 0.0001;
            test[4] = 0.00001;
            test[5] = 0.000001;
            test[6] = 0;

        for(int i=3; i<7; i++){

       createHistos(0, test[i], "cL",true,bkgdVal[0],sigVal[0]);
//       createHistos(2, 0,1,0,0,true,(970*efficiency*3000/ratioSB),(970*efficiency*3000));
//       createHistos(1, 0,1,0,0,true,(3727*efficiency*15000/ratioSB),(3727*efficiency*15000));
//       createHistos(0, 0,1,0,0,true,(34810*efficiency*15000/ratioSB),(34810*efficiency*15000));


            RooStats::HistFactory::Measurement meas("gTT", "stats");

            meas.SetOutputFilePrefix("stats/statsgTT");
            meas.SetExportOnly(false);

            //meas.SetPOI("SigXsecOverSM");
            meas.AddPOI("muSME");

            meas.SetLumi(1);
            meas.SetLumiRelErr(0.02);

            RooStats::HistFactory::Channel chan("channel");

            chan.SetData("hAsimovNullHyp", "stats/statsgTT.root");

            RooStats::HistFactory::Sample ttbar("hTTbarSM", "hTTbarSM", "stats/statsgTT.root");
            ttbar.AddNormFactor("SigXsecOverSM", 1, -1, 1);
            ttbar.AddOverallSys("syst1",  0.95, 1.05);
            chan.AddSample(ttbar); 

            RooStats::HistFactory::Sample bkgd("hBkgd", "hBkgd", "stats/statsgTT.root");
            bkgd.AddOverallSys("syst2", 0.7, 1.3 );
            chan.AddSample(bkgd);        

            RooStats::HistFactory::Sample SME("hSME", "hSME", "stats/statsgTT.root");
            SME.AddNormFactor("muSME", 1, -1, 1);
            chan.AddSample(SME);  


            meas.AddChannel(chan);
            meas.CollectHistograms();
            meas.PrintTree();
            MakeModelAndMeasurementFast(meas);
        }

   return;

}
