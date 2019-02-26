#include "../include/kineticAnalyze.hpp"

#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TH1F.h>
#include <TSpline.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TString.h>

KineticAnalyze::KineticAnalyze(){

//_______________________________________________________//
//____________________initialisations____________________//
//_______________________________________________________//

    nExp = 6;
    time24 = 86400;
    expList= std::vector<double>(nExp);
    expList[0] = 100;    expList[1] = 27;     expList[2] = 14;
    expList[3] = 13;     expList[4] = 7;      expList[5] = 2;

    nF = std::vector<int>(nExp);
    nPqq = std::vector<int>(nExp);
    nPgg = std::vector<int>(nExp);

    AF = std::vector<TMatrixD>(nExp);
    APqq = std::vector<TMatrixD>(nExp);
    APgg = std::vector<TMatrixD>(nExp);
    Aani = std::vector<TMatrixD>(nExp);
    Afus = std::vector<TMatrixD>(nExp);
    Aprod = std::vector<TMatrixD>(nExp);
    Areal = std::vector<TMatrixD>(nExp);
    for(int i=0; i<nExp; i++){
        AF[i].ResizeTo(4,4,-1);
        APqq[i].ResizeTo(4,4,-1);
        APgg[i].ResizeTo(4,4,-1);
        Aani[i].ResizeTo(4,4,-1);
        Afus[i].ResizeTo(4,4,-1);
        Aprod[i].ResizeTo(4,4,-1);
        Areal[i].ResizeTo(4,4,-1);
    }

    aL = std::vector< std::vector<double> >(nExp);
    aR = std::vector< std::vector<double> >(nExp);
    bL = std::vector< std::vector<double> >(nExp);
    bR = std::vector< std::vector<double> >(nExp);
    for(int i=0; i<nExp; i++){
        aL[i] = std::vector<double>(6);
        aR[i] = std::vector<double>(6);
        bL[i] = std::vector<double>(4);
        bR[i] = std::vector<double>(4);
    }// first index "experiment", second index "a coefficient"


    fLList = std::vector< std::vector< std::vector<double> > >(nExp);
    fRList = std::vector< std::vector< std::vector<double> > >(nExp);
    fCList = std::vector< std::vector< std::vector<double> > >(nExp);
    fDList = std::vector< std::vector< std::vector<double> > >(nExp);
    gLList = std::vector< std::vector< std::vector<double> > >(nExp);
    gRList = std::vector< std::vector< std::vector<double> > >(nExp);
    gCList = std::vector< std::vector< std::vector<double> > >(nExp);
    gDList = std::vector< std::vector< std::vector<double> > >(nExp);
    amplitudeLListf = std::vector< std::vector<double> >(nExp);
    amplitudeRListf = std::vector< std::vector<double> >(nExp);
    amplitudeCListf = std::vector< std::vector<double> >(nExp);
    amplitudeDListf = std::vector< std::vector<double> >(nExp);
    amplitudeLListg = std::vector< std::vector<double> >(nExp);
    amplitudeRListg = std::vector< std::vector<double> >(nExp);
    amplitudeCListg = std::vector< std::vector<double> >(nExp);
    amplitudeDListg = std::vector< std::vector<double> >(nExp);
    for(int i=0; i<nExp; i++){
        amplitudeLListf[i] = std::vector<double>(4);
        amplitudeRListf[i] = std::vector<double>(4);
        amplitudeCListf[i] = std::vector<double>(4);
        amplitudeDListf[i] = std::vector<double>(4);
        amplitudeLListg[i] = std::vector<double>(4);
        amplitudeRListg[i] = std::vector<double>(4);
        amplitudeCListg[i] = std::vector<double>(4);
        amplitudeDListg[i] = std::vector<double>(4);
        fLList[i] = std::vector< std::vector<double> >(4);
        fRList[i] = std::vector< std::vector<double> >(4);
        fCList[i] = std::vector< std::vector<double> >(4);
        fDList[i] = std::vector< std::vector<double> >(4);
        gLList[i] = std::vector< std::vector<double> >(4);
        gRList[i] = std::vector< std::vector<double> >(4);
        gCList[i] = std::vector< std::vector<double> >(4);
        gDList[i] = std::vector< std::vector<double> >(4);
        for(int j=0; j<4; j++){
            fLList[i][j] = std::vector<double>(time24);
            fRList[i][j] = std::vector<double>(time24);
            fCList[i][j] = std::vector<double>(time24);
            fDList[i][j] = std::vector<double>(time24);
            gLList[i][j] = std::vector<double>(time24);
            gRList[i][j] = std::vector<double>(time24);
            gCList[i][j] = std::vector<double>(time24);
            gDList[i][j] = std::vector<double>(time24);
        }
    }// first index "experiment", second index "XX(0), XY(1), XZ(2), YZ(3)", third index "time stamp"

//_________________________________________________________________________________//
// read in txt files the values (for more informations go to line 70 of this file) //
//_________________________________________________________________________________//

// 100 TeV                                      //13 TeV
    nF[0] = readNumber("100TeVCMSF");           nF[3] = readNumber("13TeVCMSF");
    nPqq[0] = readNumber("100TeVCMSPqqbar");    nPqq[3] = readNumber("13TeVCMSPqqbar");
    nPgg[0] = readNumber("100TeVCMSP2g");       nPgg[3] = readNumber("13TeVCMSP2g");
    AF[0] = readMatrix("100TeVCMSF");           AF[3] = readMatrix("13TeVCMSF");
    APqq[0] = readMatrix("100TeVCMSPqqbar");    APqq[3] = readMatrix("13TeVCMSPqqbar");
    APgg[0] = readMatrix("100TeVCMSP2g");       APgg[3] = readMatrix("13TeVCMSP2g");

// 27 TeV                                       //7 TeV
    nF[1] = readNumber("27TeVCMSF");            nF[4] = readNumber("7TeVCMSF");
    nPqq[1] = readNumber("27TeVCMSPqqbar");     nPqq[4] = readNumber("7TeVCMSPqqbar");
    nPgg[1] = readNumber("27TeVCMSP2g");        nPgg[4] = readNumber("7TeVCMSP2g");
    AF[1] = readMatrix("27TeVCMSF");           AF[4] = readMatrix("7TeVCMSF");
    APqq[1] = readMatrix("27TeVCMSPqqbar");    APqq[4] = readMatrix("7TeVCMSPqqbar");
    APgg[1] = readMatrix("27TeVCMSP2g");       APgg[4] = readMatrix("7TeVCMSP2g");

// 14 TeV                                       //1.96TeV
    nF[2] = readNumber("14TeVCMSF");            nF[5] = readNumber("2TeVCMSF");
    nPqq[2] = readNumber("14TeVCMSPqqbar");     nPqq[5] = readNumber("2TeVCMSPqqbar");
    nPgg[2] = readNumber("14TeVCMSP2g");        nPgg[5] = readNumber("2TeVCMSP2g");
    AF[2] = readMatrix("14TeVCMSF");           AF[5] = readMatrix("2TeVCMSF");
    APqq[2] = readMatrix("14TeVCMSPqqbar");    APqq[5] = readMatrix("2TeVCMSPqqbar");
    APgg[2] = readMatrix("14TeVCMSP2g");       APgg[5] = readMatrix("2TeVCMSP2g");

    for(int i=0; i<nExp; i++){
        Aani[i] = calculateAverageMatrix(1, nPqq[i], nPgg[i], nF[i], APqq[i], APgg[i], AF[i]);
        Afus[i] = calculateAverageMatrix(2, nPqq[i], nPgg[i], nF[i], APqq[i], APgg[i], AF[i]);
        Aprod[i] = calculateAverageMatrix(3, nPqq[i], nPgg[i], nF[i], APqq[i], APgg[i], AF[i]);
        Areal[i] = calculateAverageMatrix(0, nPqq[i], nPgg[i], nF[i], APqq[i], APgg[i], AF[i]);
        for(int j=1; j<6; j++){ //carefull index start at 1
            aL[i][j] = calculateCoefficent_a(j, Areal[i], latitudeCMS, azimuthCMS, true);
            aR[i][j] = calculateCoefficent_a(j, Aprod[i], latitudeCMS, azimuthCMS, true);
        }
        for(int j=0; j<4; j++){
            bL[i][j] = calculateCoefficent_b(j+1, Areal[i], latitudeCMS, azimuthCMS);
            bR[i][j] = calculateCoefficent_b(j+1, Aprod[i], latitudeCMS, azimuthCMS);
        }
        aL[i][0] = (aL[i][1]-aL[i][2])/2.;
        aR[i][0] = (aR[i][1]-aR[i][2])/2.;
        for(int j=0; j<4; j++){
            fLList[i][j] = calculatef(time24, i, j, cmunuLConst, 0, 0, 0);
            fRList[i][j] = calculatef(time24, i, j, 0, cmunuRConst, 0, 0);
            fCList[i][j] = calculatef(time24, i, j, 0, 0, cmunuConst, 0);
            fDList[i][j] = calculatef(time24, i, j, 0, 0, 0, dmunuConst);
            gLList[i][j] = calculateg(time24, i, j, cmunuLConst, 0, 0, 0);
            gRList[i][j] = calculateg(time24, i, j, 0, cmunuRConst, 0, 0);
            gCList[i][j] = calculateg(time24, i, j, 0, 0, cmunuConst, 0);
            gDList[i][j] = calculateg(time24, i, j, 0, 0, 0, dmunuConst);
            amplitudeLListf[i][j] = calculateMax(fLList[i][j]);
            amplitudeRListf[i][j] = calculateMax(fRList[i][j]);
            amplitudeCListf[i][j] = calculateMax(fCList[i][j]);
            amplitudeDListf[i][j] = calculateMax(fDList[i][j]);
            amplitudeLListg[i][j] = calculateMax(gLList[i][j]);
            amplitudeRListg[i][j] = calculateMax(gRList[i][j]);
            amplitudeCListg[i][j] = calculateMax(gCList[i][j]);
            amplitudeDListg[i][j] = calculateMax(gDList[i][j]);
        }
    }
/*
for(int i=0; i<nExp; i++)
    for(int j=0; j<time24; j++)
        if(j%500==0)
            std::cout<<fLList[i][1][j]<<std::endl;
*/
}



//_________________________________________________//
//__________________ Read Values __________________//
//_________________________________________________//

// Text input type : 
/*
________________________________
| Title                        |
|                              |
| Events n= 'number of events' |
|                              |
| m00 m01 m02 m03              |
| m10 m11 m12 m13              |
| m20 m21 m22 m23              |
| m30 m31 m32 m33              |
|______________________________|

*/

TMatrixD KineticAnalyze::readMatrix(TString s){
    TMatrixD foo(4,4);
    std::string tmp;
    int tmp2;
    std::fstream f("data/matrix/"+s+".txt", std::ios::in);
    f>>tmp>>tmp>>tmp;
    f>>tmp2;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            f>>foo(i,j);
    return foo;
}

int KineticAnalyze::readNumber(TString s){
    int foo;
    std::string tmp;
    std::fstream f("data/matrix/"+s+".txt", std::ios::in);
    f>>tmp>>tmp>>tmp;
    f>>foo;
    return foo;
}

double* KineticAnalyze::extractArray(TString s, int n){
    double* foo = new double[n]{0};
    std::fstream f("data/pT/"+s+".txt", std::ios::in);
    for(int i=0; i<n; i++)
        f>>foo[i];
    return foo;
}

//_________________________________________________//
//______________ Calculated Values ________________//
//_________________________________________________//

// Matrix A
TMatrixD KineticAnalyze::calculateAverageMatrix(int option, int nP1_user, int nP2_user, int nF_user, TMatrixD AP1_user, TMatrixD AP2_user, TMatrixD AF_user){
    // option : (1 = 1+decay), (2 = 2+decay), (0 = total)
    TMatrixD foo(4,4);
    if(option == 1)
        foo = 0.5*AP1_user + AF_user;
    else if(option == 2)
        foo = 0.5*AP2_user + AF_user;
    else if(option == 3)
        foo = ((double)nP1_user / (double)nF_user) * 0.5*AP1_user + ((double)nP2_user / (double)nF_user) * 0.5*AP2_user;
    else if(option == 0)
        foo = ((double)nP1_user / (double)nF_user) * 0.5*AP1_user + ((double)nP2_user / (double)nF_user) * 0.5*AP2_user + AF_user;
    else
        std::cout<<"Error with matrix option"<<std::endl;
    return foo;
}

double KineticAnalyze::calculateCoefficent_a(int i, TMatrixD m_user, double latitude_user, double azimuth_user, bool isCMS){
    double foo = 0;
    if(isCMS){
        if(i==1)
            foo = (sin(latitude_user)*sin(latitude_user)*sin(azimuth_user)*sin(azimuth_user) + cos(latitude_user)*cos(latitude_user)) * m_user(0,0) + (sin(latitude_user)*sin(latitude_user)*cos(azimuth_user)*cos(azimuth_user)) * m_user(2,2);
        else if(i==2)
            foo = (cos(azimuth_user)*cos(azimuth_user)) * m_user(0,0) + (sin(azimuth_user)*sin(azimuth_user)) * m_user(2,2);
        else if(i==3)
            foo = (sin(latitude_user)*cos(azimuth_user)*sin(azimuth_user)) * (m_user(2,2) - m_user(0,0));
        else if(i==4)
            foo = (cos(latitude_user)*sin(latitude_user)*cos(azimuth_user)*cos(azimuth_user)) * (m_user(2,2) - m_user(0,0));
        else if(i==5)
            foo = (cos(azimuth_user)*cos(latitude_user)*sin(azimuth_user)) * (m_user(2,2) - m_user(0,0));
        else{
            std::cout<<"error with coefficient"<<std::endl;
        }
    }
    else{
        if(i==1)
            foo = (c1D0*c1D0*c2D0*c2D0 + s1D0*s1D0) * m_user(0,0) + (c1D0*c1D0*s2D0*s2D0) * m_user(2,2);
        else if(i==2)
            foo = (s2D0*s2D0) * m_user(0,0) + (c2D0*c2D0) * m_user(2,2);
        else if(i==3)
            foo = (c1D0*c2D0*s2D0) * (m_user(2,2) - m_user(0,0));
        else if(i==4)
            foo = (c1D0*s1D0*s2D0*s2D0) * (m_user(0,0) + m_user(2,2));
        else if(i==5)
            foo = (c2D0*s1D0*s2D0) * (m_user(0,0) - m_user(2,2));
        else{
            std::cout<<"error with coefficient D0"<<std::endl;
        }
    }
    return foo;
}

double KineticAnalyze::calculateCoefficent_b(int i, TMatrixD m_user, double latitude_user, double azimuth_user){
    double foo = 0;
        if(i==1)
            foo = -sin(latitudeCMS)*sin(azimuthCMS)*m_user(3,0) + cos(latitudeCMS)*m_user(3,1) + sin(latitudeCMS)*cos(azimuthCMS)*m_user(3,2);
        else if(i==2)
            foo = cos(azimuthCMS)*m_user(3,0) + sin(azimuthCMS)*m_user(3,2);
        else if(i==3)
            foo = -cos(latitudeCMS)*sin(azimuthCMS)*m_user(3,0) - sin(latitudeCMS)*m_user(3,1) + cos(latitudeCMS)*cos(azimuthCMS)*m_user(3,2);
        else if(i==4)
            foo = (cos(latitudeCMS)*cos(latitudeCMS)*sin(azimuthCMS)*sin(azimuthCMS) + sin(latitudeCMS)*sin(latitudeCMS))*m_user(0,0) + cos(latitudeCMS)*cos(latitudeCMS)*cos(azimuthCMS)*cos(azimuthCMS)*m_user(2,2) + m_user(3,3);
        else
            std::cout<<"error with coefficient"<<std::endl;
    return foo;
}

std::vector<double> KineticAnalyze::calculatef(int n, int exp,int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu){
    std::vector<double> foo(n);
    for(int i=0; i<n; i++){
        if(cmunu==0 and dmunu==0){
            if(wilson==0)
                foo[i] = cmunuL *  2 * (aL[exp][0] * cos(2*omega*i) + aL[exp][3] * sin(2*omega*i)) + cmunuR * 2 * (aR[exp][0] * cos(2*omega*i) + aR[exp][3] * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunuL *  2 * (aL[exp][0] * sin(2*omega*i) - aL[exp][3] * cos(2*omega*i)) + cmunuR * 2 * (aR[exp][0] * sin(2*omega*i) - aR[exp][3] * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunuL *  2 * (aL[exp][4] * cos(omega*i) + aL[exp][5] * sin(omega*i)) + cmunuR * 2 * (aR[exp][4] * cos(omega*i) + aR[exp][5] * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunuL *  2 * (aL[exp][4] * sin(omega*i) - aL[exp][5] * cos(omega*i)) + cmunuR * 2 * (aR[exp][4] * sin(omega*i) - aR[exp][5] * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else if(cmunuL==0 and cmunuR==0){
            if(wilson==0)
                foo[i] = cmunu * 2 * ((aL[exp][0] + aR[exp][0]) * cos(2*omega*i) + (aL[exp][3] + aR[exp][3]) * sin(2*omega*i))
                + dmunu * 2 * ((aL[exp][0] - aR[exp][0]) * cos(2*omega*i) + (aL[exp][3] - aR[exp][3]) * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunu * 2 * ((aL[exp][0] + aR[exp][0]) * sin(2*omega*i) - (aL[exp][3] + aR[exp][3]) * cos(2*omega*i))
                    + dmunu * 2 * ((aL[exp][0] - aR[exp][0]) * sin(2*omega*i) - (aL[exp][3] - aR[exp][3]) * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunu * 2 * ((aL[exp][4] + aR[exp][4]) * cos(omega*i) + (aL[exp][5] + aR[exp][5]) * sin(omega*i))
                    + dmunu * 2 * ((aL[exp][4] - aR[exp][4]) * cos(omega*i) + (aL[exp][5] - aR[exp][5]) * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunu * 2 * ((aL[exp][4] + aR[exp][4]) * sin(omega*i) - (aL[exp][5] + aR[exp][5]) * cos(omega*i))
                    + dmunu * 2 * ((aL[exp][4] - aR[exp][4]) * sin(omega*i) - (aL[exp][5] - aR[exp][5]) * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else
            std::cout<<"error with cmunu coeff with f function"<<std::endl;
    }
    return foo;
}

std::vector<double> KineticAnalyze::calculateg(int n, int exp,int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu){
    std::vector<double> foo(n);
    for(int i=0; i<n; i++){
        if(cmunu==0 and dmunu==0){
            if(wilson==0)
                foo[i] = cmunuL * 2 * (bL[exp][0] * cos(omega*i) + bL[exp][1] * sin(omega*i)) + cmunuR * 2 * (bR[exp][0] * cos(omega*i) + bR[exp][1] * sin(omega*i));
            else if(wilson==1)
                foo[i] = cmunuL * 2 * (bL[exp][0] * sin(omega*i) - bL[exp][1] * cos(omega*i)) + cmunuR * 2 * (bR[exp][0] * sin(omega*i) - bR[exp][1] * cos(omega*i));
            else if(wilson==2)
                foo[i] = cmunuL * 2 * bL[exp][2] + cmunuR * 2 * bR[exp][2];
            else if(wilson==3)
                foo[i] = cmunuL * 2 * bL[exp][3] + cmunuR * 2 * bR[exp][3];
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else if(cmunuL==0 and cmunuR==0){
            if(wilson==0)
                foo[i] = 2 * cmunu * ((bL[exp][0] * cos(omega*i) + bL[exp][1] * sin(omega*i)) + (bR[exp][0] * cos(omega*i) + bR[exp][1] * sin(omega*i)))
                + 2 * dmunu * ((bL[exp][0] * cos(omega*i) + bL[exp][1] * sin(omega*i)) - (bR[exp][0] * cos(omega*i) + bR[exp][1] * sin(omega*i)));
            else if(wilson==1)
                foo[i] = 2 * cmunu * ((bL[exp][0] * sin(omega*i) - bL[exp][1] * cos(omega*i)) + (bR[exp][0] * sin(omega*i) - bR[exp][1] * cos(omega*i)))
                + 2 * dmunu * ((bL[exp][0] * sin(omega*i) - bL[exp][1] * cos(omega*i)) - (bR[exp][0] * sin(omega*i) - bR[exp][1] * cos(omega*i)));
            else if(wilson==2)
                foo[i] = 2 * cmunu * (bL[exp][2] + bR[exp][2]) + 2 * dmunu * (bL[exp][2] - bR[exp][2]);
            else if(wilson==3)
                foo[i] = 2 * cmunu * (bL[exp][3] + bR[exp][3]) + 2 * dmunu * (bL[exp][3] - bR[exp][3]);
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else
            std::cout<<"error with cmunu coeff with f function"<<std::endl;
    }
    return foo;
}

double KineticAnalyze::calculateMax(std::vector<double> vector_user){
    double foo = 0;
    for(int i=0; i<time24; i++){
        if(vector_user[i] >= foo)
            foo = vector_user[i];
    }
    return foo;
}

//_________________________________________________//
//________________ Histos and graphs ______________//
//_________________________________________________//


void KineticAnalyze::fTime(int munu, int exp){
    // histograms of fSME in time for each coefficients cmunu, Wilson is the choice of cmunu left/right or cmunu/dmunu
    TCanvas* w = new TCanvas("","",200,10,800,600);
    TH1F* h1 = new TH1F("","", time24, 0, 24);
    TH1F* h2 = new TH1F("","", time24, 0, 24);
    TH1F* h3 = new TH1F("","", time24, 0, 24);
    TH1F* h4 = new TH1F("","", time24, 0, 24);
    for (int k=0; k<86400; k++){
        if(munu==0){
            h1->SetBinContent(k+1, fLList[exp][0][k]);
            h2->SetBinContent(k+1, fLList[exp][1][k]);
            h3->SetBinContent(k+1, fLList[exp][2][k]);
            h4->SetBinContent(k+1, fLList[exp][3][k]);
        }
        else if(munu==1){
            h1->SetBinContent(k+1, fRList[exp][0][k]);
            h2->SetBinContent(k+1, fRList[exp][1][k]);
            h3->SetBinContent(k+1, fRList[exp][2][k]);
            h4->SetBinContent(k+1, fRList[exp][3][k]);
        }
        else if(munu==2){
            h1->SetBinContent(k+1, fCList[exp][0][k]);
            h2->SetBinContent(k+1, fCList[exp][1][k]);
            h3->SetBinContent(k+1, fCList[exp][2][k]);
            h4->SetBinContent(k+1, fCList[exp][3][k]);
        }
        else if(munu==3){
            h1->SetBinContent(k+1, fDList[exp][0][k]);
            h2->SetBinContent(k+1, fDList[exp][1][k]);
            h3->SetBinContent(k+1, fDList[exp][2][k]);
            h4->SetBinContent(k+1, fDList[exp][3][k]);
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
    }

    h1->Draw();          h1->Write();   h1->SetLineWidth(2);   h1->SetLineColor(kRed);
    h2->Draw("SAME");    h2->Write();   h2->SetLineWidth(2);   h2->SetLineColor(kMagenta);
    h3->Draw("SAME");    h3->Write();   h3->SetLineWidth(2);   h3->SetLineColor(kBlue);
    h4->Draw("SAME");    h4->Write();   h4->SetLineWidth(2);   h4->SetLineColor(kGreen);

    h1->GetYaxis()->SetTitle("f_{SME}(t)");
    h1->GetXaxis()->SetTitle("sideral time #hat{t} (in h)");
    h1->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) functions histogram ","C"); // option "C" allows to center the header
    legend->AddEntry(h1,"c_{XX} = - c_{YY}","l"); // option "l" is for line (form of legend)
    legend->AddEntry(h2,"c_{XY} = c_{YX}","l");
    legend->AddEntry(h3,"c_{XZ} = c_{ZX}","l");
    legend->AddEntry(h4,"c_{YZ} = c_{ZY}","l");
    legend->Draw();

    w->Update();
    if(munu==0){
        h1->SetTitle("c_{L#mu#nu} = 0.1");
        w->SaveAs("results/fComparaisonL.png");
    }
    else if(munu==1){
        h1->SetTitle("c_{R#mu#nu} = 0.1");
        w->SaveAs("results/fComparaisonR.png");
    }
    else if(munu==2){
        h1->SetTitle("c_{#mu#nu} = 0.1");
        w->SaveAs("results/fComparaisonC.png");
    }
    else if(munu==3){
        h1->SetTitle("d_{#mu#nu} = 0.1");
        w->SaveAs("results/fComparaisonD.png");
    }
    else
        std::cout<<"error with cmunu"<<std::endl;
}

void KineticAnalyze::gTime(int munu, int exp){
    // histograms of fSME in time for each coefficients cmunu, Wilson is the choice of cmunu left/right or cmunu/dmunu
    TCanvas* w = new TCanvas("","",200,10,800,600);
    TH1F* h1 = new TH1F("","", time24, 0, 24);
    TH1F* h2 = new TH1F("","", time24, 0, 24);
    TH1F* h3 = new TH1F("","", time24, 0, 24);
    for (int k=0; k<86400; k++){
        if(munu==0){
            h1->SetBinContent(k+1, gLList[exp][0][k]);
            h2->SetBinContent(k+1, gLList[exp][1][k]);
            h3->SetBinContent(k+1, gLList[exp][2][k]);
        }
        else if(munu==1){
            h1->SetBinContent(k+1, gRList[exp][0][k]);
            h2->SetBinContent(k+1, gRList[exp][1][k]);
            h3->SetBinContent(k+1, gRList[exp][2][k]);
        }
        else if(munu==2){
            h1->SetBinContent(k+1, gCList[exp][0][k]);
            h2->SetBinContent(k+1, gCList[exp][1][k]);
            h3->SetBinContent(k+1, gCList[exp][2][k]);
        }
        else if(munu==3){
            h1->SetBinContent(k+1, gDList[exp][0][k]);
            h2->SetBinContent(k+1, gDList[exp][1][k]);
            h3->SetBinContent(k+1, gDList[exp][2][k]);
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
    }

    h1->Draw();          h1->Write();   h1->SetLineWidth(2);   h1->SetLineColor(kRed);
    h2->Draw("SAME");    h2->Write();   h2->SetLineWidth(2);   h2->SetLineColor(kMagenta);
    h3->Draw("SAME");    h3->Write();   h3->SetLineWidth(2);   h3->SetLineColor(kBlue);

    h1->GetYaxis()->SetTitle("f_{SME}(t)");
    h1->GetXaxis()->SetTitle("sideral time #hat{t} (in h)");
    h1->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) functions histogram ","C"); // option "C" allows to center the header
    legend->AddEntry(h1,"c_{TX} = c_{XT}","l"); // option "l" is for line (form of legend)
    legend->AddEntry(h2,"c_{TY} = c_{YT}","l");
    legend->AddEntry(h3,"c_{TZ} = c_{ZT}","l");
    legend->Draw();

    w->Update();
    if(munu==0){
        h1->SetTitle("c_{L#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonL.png");
    }
    else if(munu==1){
        h1->SetTitle("c_{R#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonR.png");
    }
    else if(munu==2){
        h1->SetTitle("c_{#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonC.png");
    }
    else if(munu==3){
        h1->SetTitle("d_{#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonD.png");
    }
    else
        std::cout<<"error with cmunu"<<std::endl;
}

void KineticAnalyze::gTimeTT(int munu, int exp){
    // histograms of fSME in time for each coefficients cmunu, Wilson is the choice of cmunu left/right or cmunu/dmunu
    TCanvas* w = new TCanvas("","",200,10,800,600);
    TH1F* h1 = new TH1F("","", time24, 0, 24);
    for (int k=0; k<86400; k++){
        if(munu==0)
            h1->SetBinContent(k+1, gLList[exp][3][k]);
        else if(munu==1)
            h1->SetBinContent(k+1, gRList[exp][3][k]);
        else if(munu==2)
            h1->SetBinContent(k+1, gCList[exp][3][k]);
        else if(munu==3)
            h1->SetBinContent(k+1, gDList[exp][3][k]);
        else
            std::cout<<"error with cmunu"<<std::endl;
    }

    h1->Draw();          h1->Write();   h1->SetLineWidth(2);   h1->SetLineColor(kRed);

    h1->GetYaxis()->SetTitle("f_{SME}(t)");
    h1->GetXaxis()->SetTitle("sideral time #hat{t} (in h)");
    h1->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("f_{SME}(t) functions histogram ","C"); // option "C" allows to center the header
    legend->AddEntry(h1,"c_{TT} = c_{ZZ}","l"); // option "l" is for line (form of legend)
    legend->Draw();

    w->Update();
    if(munu==0){
        h1->SetTitle("c_{L#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonLTT.png");
    }
    else if(munu==1){
        h1->SetTitle("c_{R#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonRTT.png");
    }
    else if(munu==2){
        h1->SetTitle("c_{#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonCTT.png");
    }
    else if(munu==3){
        h1->SetTitle("d_{#mu#nu} = 0.1");
        w->SaveAs("results/gComparaisonDTT.png");
    }
    else
        std::cout<<"error with cmunu"<<std::endl;
}

void KineticAnalyze::amplEnergy(){
    TString name[4];    TString name2[4];
    name[0]="cL";    name[1]="cR";    name[2]="c";    name[3]="d";
    name2[0]="_{XX}";    name2[1]="_{XY}";    name2[2]="_{XZ}";    name2[3]="_{YZ}";

    double* x = new double[nExp];
    for(int i=0; i<nExp; i++)
        x[i] = expList[i];

    std::vector< std::vector<double*> > y(4);
    for(int i=0; i<4; i++){
        y[i] = std::vector<double*>(4);
        for(int j=0; j<4; j++){
            y[i][j] = new double[nExp];
            for(int k=0; k<nExp; k++)
                y[i][j][k] = amplitudeLListf[k][j];
        }
    }

    TCanvas*** w = new TCanvas**[4];
    TGraph*** g = new TGraph**[4];
    for(int i=0; i<4; i++){
        w[i] = new TCanvas*[4];
        g[i] = new TGraph*[4];
        for(int j=0; j<4; j++){
            w[i][j] = new TCanvas(name[i]+name2[j],"",200,10,800,600);
            g[i][j] = new TGraph(nExp, x, y[i][j]);
            g[i][j]->SetTitle(name[i]+name2[j]);
            g[i][j]->Draw("A*");
            g[i][j]->GetYaxis()->SetTitle("f_{SME}(t)");
            g[i][j]->GetXaxis()->SetTitle("sideral time #hat{t} (in h)");
            w[i][j]->SetLogx();
            w[i][j]->SaveAs("results/amplEnergy"+name[i]+name2[j]+".png");
        }
    }
}


//_________________________________________________//
//________________ statistical stuff ______________//
//_________________________________________________//

TH1F* KineticAnalyze::statHistosConst(TString name, double value){
    int bin = 24;
    TH1F* foo = new TH1F(name, name, bin, 0, bin);
    for(int i=0; i<24; i++)
        foo->SetBinContent(i+1, value/((double)bin));
    return foo;
}

TH1F* KineticAnalyze::statHistosf(TString name, int exp, double cmunuL, double cmunuR, double cmunu, double dmunu, bool isXX){
    int bin=24;
    TH1F* foo = new TH1F(name, name, bin, 0, bin);
    if(cmunuL!=0 and cmunuR==0 and cmunu==0 and dmunu==0)//L
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, fLList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, fLList[exp][2][i*3600]);
        }
    else if(cmunuL==0 and cmunuR!=0 and cmunu==0 and dmunu==0)//R
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, fRList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, fRList[exp][2][i*3600]);
        }
    else if(cmunuL==0 and cmunuR==0 and cmunu!=0 and dmunu==0)//C
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, fCList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, fCList[exp][2][i*3600]);
        }
    else if(cmunuL==0 and cmunuR==0 and cmunu==0 and dmunu!=0)//D
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, fDList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, fDList[exp][2][i*3600]);
        }
    else
        std::cout<<"error with stats"<<std::endl;

    return foo;
}


TH1F* KineticAnalyze::statHistosg(TString name, int exp, double cmunuL, double cmunuR, double cmunu, double dmunu, bool isTZ){}
TH1F* KineticAnalyze::statHistosgTT(TString name, int exp, double cmunuL, double cmunuR, double cmunu, double dmunu){}

