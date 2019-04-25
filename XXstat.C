#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "include/kineticAnalyze.hpp"
#include "src/kineticAnalyze.cpp"
#include <sstream>


using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;


void createHistos(int exp, double cmunuGen, TString wilson, bool isXX, double bkgd, double ttbar){

//___________________________________________ SME preparation ___________________________________________//

    int time = 24;              //hours
    int bin = 24;               //bins
    ostringstream sout;
    sout<<cmunuGen;
    TString s = sout.str();

    KineticAnalyze k(cmunuGen, cmunuGen, cmunuGen, cmunuGen);

    TH1F* hBkgd = k.statHistosConst("hBkgd", bkgd);
    TH1F* hTTbarSM = k.statHistosConst("hTTbarSM", ttbar);
    TH1F* hAsimovNullHyp = k.statHistosConst("hAsimovNullHyp", bkgd+ttbar);
    TH1F* hSME = k.statHistosf("hSME", exp, wilson, isXX, ttbar);
//    TH1F* hSME = k.statHistosConst("hSME", bkgd+ttbar);

   TFile* fOutput = new TFile("stats/cXX/"+s+"/statsfXX.root","RECREATE");
   hBkgd->Write();
   hTTbarSM->Write();
   hAsimovNullHyp->Write();
   hSME->Write();
   fOutput->Write();
   fOutput->Close();
}


void XXstat(){

//_____________________ All needed parameters _____________________//

    double ratioSB = 15.3576;
    double efficiency = 5.533;

    double bkgdVal[4], sigVal[4];
           bkgdVal[3] = 44045;      bkgdVal[2] = (970*efficiency*3000/ratioSB);     bkgdVal[1] = (3727*efficiency*15000/ratioSB);      bkgdVal[0] = (34810*efficiency*15000/ratioSB);
           sigVal[3] = 676431;      sigVal[2] = 970*efficiency*3000;     sigVal[1] = 3727*efficiency*15000;      sigVal[0] = 34810*efficiency*15000;

    double lumi[3];
            lumi[0]=0.02;
            lumi[1]=0.015;
            lumi[2]=0.01;

    double test[5];
            test[0] = 0.01;
            test[1] = 0.001;
            test[2] = 0.0001;
            test[3] = 0.00001;
            test[4] = 0.000001;

    for(int it=0; it<3; it++)
    for(int i=0; i<1; i++)
        for(int j=0; j<4; j++){

         createHistos(j, test[i], "D",false,bkgdVal[j],sigVal[j]);
//       createHistos(2, 0,1,0,0,true,(970*efficiency*3000/ratioSB),(970*efficiency*3000));
//       createHistos(1, 0,1,0,0,true,(3727*efficiency*15000/ratioSB),(3727*efficiency*15000));
//       createHistos(0, 0,1,0,0,true,(34810*efficiency*15000/ratioSB),(34810*efficiency*15000));

            RooStats::HistFactory::Measurement meas("fXX", "stats");

            ostringstream sout;
            sout<<"stats/cXX/"<<test[i]<<"/statsfXX";
            meas.SetOutputFilePrefix(sout.str());
            meas.SetExportOnly(false);

            //meas.SetPOI("SigXsecOverSM");
            meas.AddPOI("muSME");

            meas.SetLumi(1);
            meas.SetLumiRelErr(lumi[it]);

            RooStats::HistFactory::Channel chan("channel");

            chan.SetData("hAsimovNullHyp", sout.str()+".root");

            RooStats::HistFactory::Sample ttbar("hTTbarSM", "hTTbarSM", sout.str()+".root");
            ttbar.AddNormFactor("SigXsecOverSM", 1, -1, 1);
            ttbar.AddOverallSys("syst1",  0.95, 1.05);
            chan.AddSample(ttbar);

            RooStats::HistFactory::Sample bkgd("hBkgd", "hBkgd", sout.str()+".root");
            bkgd.AddOverallSys("syst2", 0.7, 1.3 );
            chan.AddSample(bkgd);

            RooStats::HistFactory::Sample SME("hSME", "hSME", sout.str()+".root");
            SME.AddNormFactor("muSME", 1, -1, 1);
            chan.AddSample(SME);


            meas.AddChannel(chan);
            meas.CollectHistograms();
            meas.PrintTree();
            MakeModelAndMeasurementFast(meas);
        }

   return;

}
