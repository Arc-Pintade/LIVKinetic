#include  "../include/kineticAnalyze.hpp"

#include <TCanvas.h>
#include <TF1.h>

#include <iostream>

using namespace std;

int main ()
{

    TCanvas *c = new TCanvas("", "", 450, 450);
    double wilson = 0.0001;
    string XX[4];
    XX[0] = "XX";    XX[1] = "XY";    XX[2] = "XZ";    XX[3] = "YZ";

    KineticAnalyze k(wilson,wilson,wilson,wilson);

    double ratioSB = 15.3576;
    double efficiency = 5.533;
    double range = 10;
    double step = 0.0001;

//doi:10.1140/epjc/s10052-017-4718-8.
    const double bkgdWitness = 44405;
    const double sigWitness  = 676431;
    //double bkgdWitness = 34810*efficiency*15000/ratioSB;
    //double sigWitness  = 34810*efficiency*15000;

    // 3 = LHC RunII
    TH1F* hSME = k.statHistosf("sigSME", 3, "L", true, sigWitness);

    TH1F* hXSquare = new TH1F("X^2", "", 2*range/step , -range, range);
    TF1* fitHist = new TF1("pol2","pol2", -range, range);

    double valueAzimov = bkgdWitness/24.;

    double systematicBkgd = 1;
    double systematicSignal = 1;
    double luminosity = 1.0;

    double systematicTotal = sqrt(luminosity*(valueAzimov*systematicBkgd));

    hXSquare->GetXaxis()->SetTitle("c_{#mu#nu}");
    hXSquare->GetYaxis()->SetTitle("#Delta#chi^{2}");
    hXSquare->GetXaxis()->CenterTitle(kTRUE);
    hXSquare->GetYaxis()->CenterTitle(kTRUE);
    hXSquare->SetStats(kFALSE);

    double divideSignal = 1.0/(systematicTotal*systematicTotal*24);

    int bin = 1;
    for (double i = -range ; i<range-0.0001; i+=0.0001){

      double dileptChannel = 0;
      for (int p=1; p<26; p++){
        dileptChannel += (valueAzimov
                        -luminosity*((i*hSME->GetBinContent(p)+1)*systematicSignal
                                        +systematicBkgd))
                         *(valueAzimov
                        -luminosity*((i*hSME->GetBinContent(p)+1)*systematicSignal
                                        +systematicBkgd));
      }
      hXSquare->SetBinContent(bin,dileptChannel*divideSignal+dileptChannel);
      bin++;
    }

    hXSquare->Fit(fitHist,"R");
    cout<<"Delta Chi^2 = "
        <<fitHist->GetX(hXSquare->GetBinContent(hXSquare->GetMinimumBin())+1)
        <<endl;
    hXSquare->Draw();
    c->Update();
    c->SaveAs("stats/Xsquare/X.png");

  return 0;
}
