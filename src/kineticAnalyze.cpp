#include "../include/kineticAnalyze.hpp"

#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TSpline.h>
#include <TMultiGraph.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TH2F.h>

KineticAnalyze::KineticAnalyze(double cmunuL_user, double cmunuR_user, double cmunu_user, double dmunu_user){

//_______________________________________________________//
//____________________initialisations____________________//
//_______________________________________________________//

    cmunuL_class = cmunuL_user;
    cmunuR_class = cmunuR_user;
    cmunu_class = cmunu_user;
    dmunu_class = dmunu_user;

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
//-------------------------------------------- Article ---------------------------------------//
    aLTEV = std::vector< std::vector<double> >(nExp);
    aRTEV = std::vector< std::vector<double> >(nExp);
    aLgg = std::vector< std::vector<double> >(nExp);
    aRgg = std::vector< std::vector<double> >(nExp);
    bLgg = std::vector< std::vector<double> >(nExp);
    bRgg = std::vector< std::vector<double> >(nExp);
    aLqq = std::vector< std::vector<double> >(nExp);
    aRqq = std::vector< std::vector<double> >(nExp);
    bLqq = std::vector< std::vector<double> >(nExp);
    bRqq = std::vector< std::vector<double> >(nExp);
//-------------------------------------------------------------------------------------------//
    for(int i=0; i<nExp; i++){
        aL[i] = std::vector<double>(6);
        aR[i] = std::vector<double>(6);
        bL[i] = std::vector<double>(4);
        bR[i] = std::vector<double>(4);
//-------------------------------------------- Article ---------------------------------------//
        aLTEV[i] = std::vector<double>(6);
        aRTEV[i] = std::vector<double>(6);
        aLgg[i] = std::vector<double>(6);
        aRgg[i] = std::vector<double>(6);
        bLgg[i] = std::vector<double>(4);
        bRgg[i] = std::vector<double>(4);
        aLqq[i] = std::vector<double>(6);
        aRqq[i] = std::vector<double>(6);
        bLqq[i] = std::vector<double>(4);
        bRqq[i] = std::vector<double>(4);
//-------------------------------------------------------------------------------------------//
    }// first index "experiment", second index "a coefficient"

    fLList = std::vector< std::vector< std::vector<double> > >(nExp);
    fRList = std::vector< std::vector< std::vector<double> > >(nExp);
    fCList = std::vector< std::vector< std::vector<double> > >(nExp);
    fDList = std::vector< std::vector< std::vector<double> > >(nExp);
    gLList = std::vector< std::vector< std::vector<double> > >(nExp);
    gRList = std::vector< std::vector< std::vector<double> > >(nExp);
    gCList = std::vector< std::vector< std::vector<double> > >(nExp);
    gDList = std::vector< std::vector< std::vector<double> > >(nExp);
//-------------------------------------------- Article ---------------------------------------//
    fLListTEV = std::vector< std::vector< std::vector<double> > >(nExp);
    fRListTEV = std::vector< std::vector< std::vector<double> > >(nExp);
    fCListTEV = std::vector< std::vector< std::vector<double> > >(nExp);
    fDListTEV = std::vector< std::vector< std::vector<double> > >(nExp);
    fLListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    fRListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    fCListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    fDListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    gLListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    gRListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    gCListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    gDListgg = std::vector< std::vector< std::vector<double> > >(nExp);
    fLListqq = std::vector< std::vector< std::vector<double> > >(nExp);
    fRListqq = std::vector< std::vector< std::vector<double> > >(nExp);
    fCListqq = std::vector< std::vector< std::vector<double> > >(nExp);
    fDListqq = std::vector< std::vector< std::vector<double> > >(nExp);
    gLListqq = std::vector< std::vector< std::vector<double> > >(nExp);
    gRListqq = std::vector< std::vector< std::vector<double> > >(nExp);
    gCListqq = std::vector< std::vector< std::vector<double> > >(nExp);
    gDListqq = std::vector< std::vector< std::vector<double> > >(nExp);
//-------------------------------------------------------------------------------------------//
    amplitudeLListf = std::vector< std::vector<double> >(nExp);
    amplitudeRListf = std::vector< std::vector<double> >(nExp);
    amplitudeCListf = std::vector< std::vector<double> >(nExp);
    amplitudeDListf = std::vector< std::vector<double> >(nExp);
    amplitudeLListg = std::vector< std::vector<double> >(nExp);
    amplitudeRListg = std::vector< std::vector<double> >(nExp);
    amplitudeCListg = std::vector< std::vector<double> >(nExp);
    amplitudeDListg = std::vector< std::vector<double> >(nExp);
//-------------------------------------------- Article ---------------------------------------//
    amplitudeLListfTEV = std::vector< std::vector<double> >(nExp);
    amplitudeRListfTEV = std::vector< std::vector<double> >(nExp);
    amplitudeCListfTEV = std::vector< std::vector<double> >(nExp);
    amplitudeDListfTEV = std::vector< std::vector<double> >(nExp);
//-------------------------------------------------------------------------------------------//
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
//-------------------------------------------- Article ---------------------------------------//
        fLListTEV[i] = std::vector< std::vector<double> >(4);
        fRListTEV[i] = std::vector< std::vector<double> >(4);
        fCListTEV[i] = std::vector< std::vector<double> >(4);
        fDListTEV[i] = std::vector< std::vector<double> >(4);
        fLListgg[i] = std::vector< std::vector<double> >(4);
        fRListgg[i] = std::vector< std::vector<double> >(4);
        fCListgg[i] = std::vector< std::vector<double> >(4);
        fDListgg[i] = std::vector< std::vector<double> >(4);
        gLListgg[i] = std::vector< std::vector<double> >(4);
        gRListgg[i] = std::vector< std::vector<double> >(4);
        gCListgg[i] = std::vector< std::vector<double> >(4);
        gDListgg[i] = std::vector< std::vector<double> >(4);
        fLListqq[i] = std::vector< std::vector<double> >(4);
        fRListqq[i] = std::vector< std::vector<double> >(4);
        fCListqq[i] = std::vector< std::vector<double> >(4);
        fDListqq[i] = std::vector< std::vector<double> >(4);
        gLListqq[i] = std::vector< std::vector<double> >(4);
        gRListqq[i] = std::vector< std::vector<double> >(4);
        gCListqq[i] = std::vector< std::vector<double> >(4);
        gDListqq[i] = std::vector< std::vector<double> >(4);
        amplitudeLListfTEV[i] = std::vector<double>(4);
        amplitudeRListfTEV[i] = std::vector<double>(4);
        amplitudeCListfTEV[i] = std::vector<double>(4);
        amplitudeDListfTEV[i] = std::vector<double>(4);
//-------------------------------------------------------------------------------------------//
        for(int j=0; j<4; j++){
            fLList[i][j] = std::vector<double>(time24);
            fRList[i][j] = std::vector<double>(time24);
            fCList[i][j] = std::vector<double>(time24);
            fDList[i][j] = std::vector<double>(time24);
            gLList[i][j] = std::vector<double>(time24);
            gRList[i][j] = std::vector<double>(time24);
            gCList[i][j] = std::vector<double>(time24);
            gDList[i][j] = std::vector<double>(time24);
//-------------------------------------------- Article ---------------------------------------//
            fLListTEV[i][j] = std::vector<double>(time24);
            fRListTEV[i][j] = std::vector<double>(time24);
            fCListTEV[i][j] = std::vector<double>(time24);
            fDListTEV[i][j] = std::vector<double>(time24);
            fLListgg[i][j] = std::vector<double>(time24);
            fRListgg[i][j] = std::vector<double>(time24);
            fCListgg[i][j] = std::vector<double>(time24);
            fDListgg[i][j] = std::vector<double>(time24);
            gLListgg[i][j] = std::vector<double>(time24);
            gRListgg[i][j] = std::vector<double>(time24);
            gCListgg[i][j] = std::vector<double>(time24);
            gDListgg[i][j] = std::vector<double>(time24);
            fLListqq[i][j] = std::vector<double>(time24);
            fRListqq[i][j] = std::vector<double>(time24);
            fCListqq[i][j] = std::vector<double>(time24);
            fDListqq[i][j] = std::vector<double>(time24);
            gLListqq[i][j] = std::vector<double>(time24);
            gRListqq[i][j] = std::vector<double>(time24);
            gCListqq[i][j] = std::vector<double>(time24);
            gDListqq[i][j] = std::vector<double>(time24);
//-------------------------------------------------------------------------------------------//
        }
    }// first index "experiment", second index "XX(0), XY(1), XZ(2), YZ(3)", third index "time stamp"

//_________________________________________________________________________________//
// read in txt files the values (for more informations go to line 179 of this file) //
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
            aL[i][j] = calculateCoefficent_a(j, Areal[i], latitudeCMS, azimuthCMS, true); //bool isCMS or D0
            aR[i][j] = calculateCoefficent_a(j, Aprod[i], latitudeCMS, azimuthCMS, true);
//-------------------------------------------- Article ---------------------------------------//
            aLTEV[i][j] = calculateCoefficent_a(j, Areal[i], latitudeCMS, azimuthCMS, false);
            aRTEV[i][j] = calculateCoefficent_a(j, Aprod[i], latitudeCMS, azimuthCMS, false);
            aLgg[i][j] = calculateCoefficent_a(j, Afus[i], latitudeCMS, azimuthCMS, true);
            aRgg[i][j] = calculateCoefficent_a(j, APgg[i], latitudeCMS, azimuthCMS, true);
            aLqq[i][j] = calculateCoefficent_a(j, Aani[i], latitudeCMS, azimuthCMS, true);
            aRqq[i][j] = calculateCoefficent_a(j, APqq[i], latitudeCMS, azimuthCMS, true);
//--------------------------------------------------------------------------------------------//
        }
        for(int j=0; j<4; j++){
            bL[i][j] = calculateCoefficent_b(j+1, Areal[i], latitudeCMS, azimuthCMS);
            bR[i][j] = calculateCoefficent_b(j+1, Aprod[i], latitudeCMS, azimuthCMS);
//-------------------------------------------- Article ---------------------------------------//
            bLgg[i][j] = calculateCoefficent_b(j+1, Afus[i], latitudeCMS, azimuthCMS);
            bRgg[i][j] = calculateCoefficent_b(j+1, APgg[i], latitudeCMS, azimuthCMS);
            bLqq[i][j] = calculateCoefficent_b(j+1, Aani[i], latitudeCMS, azimuthCMS);
            bRqq[i][j] = calculateCoefficent_b(j+1, APqq[i], latitudeCMS, azimuthCMS);
//--------------------------------------------------------------------------------------------//
        }
        aL[i][0] = (aL[i][1]-aL[i][2])/2.;
        aR[i][0] = (aR[i][1]-aR[i][2])/2.;
//-------------------------------------------- Article ---------------------------------------//
        aLTEV[i][0] = (aLTEV[i][1]-aLTEV[i][2])/2.;
        aRTEV[i][0] = (aRTEV[i][1]-aRTEV[i][2])/2.;
        aLgg[i][0] = (aLgg[i][1]-aLgg[i][2])/2.;
        aRgg[i][0] = (aRgg[i][1]-aRgg[i][2])/2.;
        aLqq[i][0] = (aLqq[i][1]-aLqq[i][2])/2.;
        aRqq[i][0] = (aRqq[i][1]-aRqq[i][2])/2.;
//--------------------------------------------------------------------------------------------//
        for(int j=0; j<4; j++){
            fLList[i][j] = calculatef(time24, i, j, cmunuL_class, 0, 0, 0);
            fRList[i][j] = calculatef(time24, i, j, 0, cmunuR_class, 0, 0);
            fCList[i][j] = calculatef(time24, i, j, 0, 0, cmunu_class, 0);
            fDList[i][j] = calculatef(time24, i, j, 0, 0, 0, dmunu_class);
            gLList[i][j] = calculateg(time24, i, j, cmunuL_class, 0, 0, 0);
            gRList[i][j] = calculateg(time24, i, j, 0, cmunuR_class, 0, 0);
            gCList[i][j] = calculateg(time24, i, j, 0, 0, cmunu_class, 0);
            gDList[i][j] = calculateg(time24, i, j, 0, 0, 0, dmunu_class);
//-------------------------------------------- Article ---------------------------------------//
            fLListTEV[i][j] = calculatefTEV(time24, i, j, cmunuL_class, 0, 0, 0);
            fRListTEV[i][j] = calculatefTEV(time24, i, j, 0, cmunuR_class, 0, 0);
            fCListTEV[i][j] = calculatefTEV(time24, i, j, 0, 0, cmunu_class, 0);
            fDListTEV[i][j] = calculatefTEV(time24, i, j, 0, 0, 0, dmunu_class);
            fLListgg[i][j] = calculatefgg(time24, i, j, cmunuL_class, 0, 0, 0);
            fRListgg[i][j] = calculatefgg(time24, i, j, 0, cmunuR_class, 0, 0);
            fCListgg[i][j] = calculatefgg(time24, i, j, 0, 0, cmunu_class, 0);
            fDListgg[i][j] = calculatefgg(time24, i, j, 0, 0, 0, dmunu_class);
            gLListgg[i][j] = calculateggg(time24, i, j, cmunuL_class, 0, 0, 0);
            gRListgg[i][j] = calculateggg(time24, i, j, 0, cmunuR_class, 0, 0);
            gCListgg[i][j] = calculateggg(time24, i, j, 0, 0, cmunu_class, 0);
            gDListgg[i][j] = calculateggg(time24, i, j, 0, 0, 0, dmunu_class);
            fLListqq[i][j] = calculatefqq(time24, i, j, cmunuL_class, 0, 0, 0);
            fRListqq[i][j] = calculatefqq(time24, i, j, 0, cmunuR_class, 0, 0);
            fCListqq[i][j] = calculatefqq(time24, i, j, 0, 0, cmunu_class, 0);
            fDListqq[i][j] = calculatefqq(time24, i, j, 0, 0, 0, dmunu_class);
            gLListqq[i][j] = calculategqq(time24, i, j, cmunuL_class, 0, 0, 0);
            gRListqq[i][j] = calculategqq(time24, i, j, 0, cmunuR_class, 0, 0);
            gCListqq[i][j] = calculategqq(time24, i, j, 0, 0, cmunu_class, 0);
            gDListqq[i][j] = calculategqq(time24, i, j, 0, 0, 0, dmunu_class);
//--------------------------------------------------------------------------------------------//
            amplitudeLListf[i][j] = calculateMax(fLList[i][j]);
            amplitudeRListf[i][j] = calculateMax(fRList[i][j]);
            amplitudeCListf[i][j] = calculateMax(fCList[i][j]);
            amplitudeDListf[i][j] = calculateMax(fDList[i][j]);
            amplitudeLListg[i][j] = calculateMax(gLList[i][j]);
            amplitudeRListg[i][j] = calculateMax(gRList[i][j]);
            amplitudeCListg[i][j] = calculateMax(gCList[i][j]);
            amplitudeDListg[i][j] = calculateMax(gDList[i][j]);
//-------------------------------------------- Article ---------------------------------------//
            amplitudeLListfTEV[i][j] = calculateMax(fLListTEV[i][j]);
            amplitudeRListfTEV[i][j] = calculateMax(fRListTEV[i][j]);
            amplitudeCListfTEV[i][j] = calculateMax(fCListTEV[i][j]);
            amplitudeDListfTEV[i][j] = calculateMax(fDListTEV[i][j]);
//--------------------------------------------------------------------------------------------//
        }
    }
/*
for(int i=0; i<nExp; i++)
    for(int j=0; j<time24; j++)
        if(j%500==0)
calculatef(time24, i, j, cmunuLConst, 0, 0, 0);
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
//--------------------------------------------------------------------- Article -----------------------------------------------------------------------------//
std::vector<double> KineticAnalyze::calculatefTEV(int n, int exp,int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu){
    std::vector<double> foo(n);
    for(int i=0; i<n; i++){
        if(cmunu==0 and dmunu==0){
            if(wilson==0)
                foo[i] = cmunuL *  2 * (aLTEV[exp][0] * cos(2*omega*i) + aLTEV[exp][3] * sin(2*omega*i)) + cmunuR * 2 * (aRTEV[exp][0] * cos(2*omega*i) + aRTEV[exp][3] * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunuL *  2 * (aLTEV[exp][0] * sin(2*omega*i) - aLTEV[exp][3] * cos(2*omega*i)) + cmunuR * 2 * (aRTEV[exp][0] * sin(2*omega*i) - aRTEV[exp][3] * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunuL *  2 * (aLTEV[exp][4] * cos(omega*i) + aLTEV[exp][5] * sin(omega*i)) + cmunuR * 2 * (aRTEV[exp][4] * cos(omega*i) + aRTEV[exp][5] * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunuL *  2 * (aLTEV[exp][4] * sin(omega*i) - aLTEV[exp][5] * cos(omega*i)) + cmunuR * 2 * (aRTEV[exp][4] * sin(omega*i) - aRTEV[exp][5] * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else if(cmunuL==0 and cmunuR==0){
            if(wilson==0)
                foo[i] = cmunu * 2 * ((aLTEV[exp][0] + aRTEV[exp][0]) * cos(2*omega*i) + (aLTEV[exp][3] + aRTEV[exp][3]) * sin(2*omega*i))
                + dmunu * 2 * ((aLTEV[exp][0] - aRTEV[exp][0]) * cos(2*omega*i) + (aLTEV[exp][3] - aRTEV[exp][3]) * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunu * 2 * ((aLTEV[exp][0] + aRTEV[exp][0]) * sin(2*omega*i) - (aLTEV[exp][3] + aRTEV[exp][3]) * cos(2*omega*i))
                    + dmunu * 2 * ((aLTEV[exp][0] - aRTEV[exp][0]) * sin(2*omega*i) - (aLTEV[exp][3] - aRTEV[exp][3]) * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunu * 2 * ((aLTEV[exp][4] + aRTEV[exp][4]) * cos(omega*i) + (aLTEV[exp][5] + aRTEV[exp][5]) * sin(omega*i))
                    + dmunu * 2 * ((aLTEV[exp][4] - aRTEV[exp][4]) * cos(omega*i) + (aLTEV[exp][5] - aRTEV[exp][5]) * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunu * 2 * ((aLTEV[exp][4] + aRTEV[exp][4]) * sin(omega*i) - (aLTEV[exp][5] + aRTEV[exp][5]) * cos(omega*i))
                    + dmunu * 2 * ((aLTEV[exp][4] - aRTEV[exp][4]) * sin(omega*i) - (aLTEV[exp][5] - aRTEV[exp][5]) * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else
            std::cout<<"error with cmunu coeff with f function"<<std::endl;
    }
    return foo;
}
std::vector<double> KineticAnalyze::calculatefgg(int n, int exp,int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu){
    std::vector<double> foo(n);
    for(int i=0; i<n; i++){
        if(cmunu==0 and dmunu==0){
            if(wilson==0)
                foo[i] = cmunuL *  2 * (aLgg[exp][0] * cos(2*omega*i) + aLgg[exp][3] * sin(2*omega*i)) + cmunuR * 2 * (aRgg[exp][0] * cos(2*omega*i) + aRgg[exp][3] * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunuL *  2 * (aLgg[exp][0] * sin(2*omega*i) - aLgg[exp][3] * cos(2*omega*i)) + cmunuR * 2 * (aRgg[exp][0] * sin(2*omega*i) - aRgg[exp][3] * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunuL *  2 * (aLgg[exp][4] * cos(omega*i) + aLgg[exp][5] * sin(omega*i)) + cmunuR * 2 * (aRgg[exp][4] * cos(omega*i) + aRgg[exp][5] * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunuL *  2 * (aLgg[exp][4] * sin(omega*i) - aLgg[exp][5] * cos(omega*i)) + cmunuR * 2 * (aRgg[exp][4] * sin(omega*i) - aRgg[exp][5] * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else if(cmunuL==0 and cmunuR==0){
            if(wilson==0)
                foo[i] = cmunu * 2 * ((aLgg[exp][0] + aRgg[exp][0]) * cos(2*omega*i) + (aLgg[exp][3] + aRgg[exp][3]) * sin(2*omega*i))
                + dmunu * 2 * ((aLgg[exp][0] - aRgg[exp][0]) * cos(2*omega*i) + (aLgg[exp][3] - aRgg[exp][3]) * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunu * 2 * ((aLgg[exp][0] + aRgg[exp][0]) * sin(2*omega*i) - (aLgg[exp][3] + aRgg[exp][3]) * cos(2*omega*i))
                    + dmunu * 2 * ((aLgg[exp][0] - aRgg[exp][0]) * sin(2*omega*i) - (aLgg[exp][3] - aRgg[exp][3]) * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunu * 2 * ((aLgg[exp][4] + aRgg[exp][4]) * cos(omega*i) + (aLgg[exp][5] + aRgg[exp][5]) * sin(omega*i))
                    + dmunu * 2 * ((aLgg[exp][4] - aRgg[exp][4]) * cos(omega*i) + (aLgg[exp][5] - aRgg[exp][5]) * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunu * 2 * ((aLgg[exp][4] + aRgg[exp][4]) * sin(omega*i) - (aLgg[exp][5] + aRgg[exp][5]) * cos(omega*i))
                    + dmunu * 2 * ((aLgg[exp][4] - aRgg[exp][4]) * sin(omega*i) - (aLgg[exp][5] - aRgg[exp][5]) * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else
            std::cout<<"error with cmunu coeff with f function"<<std::endl;
    }
    return foo;
}

std::vector<double> KineticAnalyze::calculatefqq(int n, int exp,int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu){
    std::vector<double> foo(n);
    for(int i=0; i<n; i++){
        if(cmunu==0 and dmunu==0){
            if(wilson==0)
                foo[i] = cmunuL *  2 * (aLqq[exp][0] * cos(2*omega*i) + aLqq[exp][3] * sin(2*omega*i)) + cmunuR * 2 * (aRqq[exp][0] * cos(2*omega*i) + aRqq[exp][3] * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunuL *  2 * (aLqq[exp][0] * sin(2*omega*i) - aLqq[exp][3] * cos(2*omega*i)) + cmunuR * 2 * (aRqq[exp][0] * sin(2*omega*i) - aRqq[exp][3] * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunuL *  2 * (aLqq[exp][4] * cos(omega*i) + aLqq[exp][5] * sin(omega*i)) + cmunuR * 2 * (aRqq[exp][4] * cos(omega*i) + aRqq[exp][5] * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunuL *  2 * (aLqq[exp][4] * sin(omega*i) - aLqq[exp][5] * cos(omega*i)) + cmunuR * 2 * (aRqq[exp][4] * sin(omega*i) - aRqq[exp][5] * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else if(cmunuL==0 and cmunuR==0){
            if(wilson==0)
                foo[i] = cmunu * 2 * ((aLqq[exp][0] + aRqq[exp][0]) * cos(2*omega*i) + (aLqq[exp][3] + aRqq[exp][3]) * sin(2*omega*i))
                + dmunu * 2 * ((aLqq[exp][0] - aRqq[exp][0]) * cos(2*omega*i) + (aLqq[exp][3] - aRqq[exp][3]) * sin(2*omega*i));
            else if(wilson==1)
                foo[i] = cmunu * 2 * ((aLqq[exp][0] + aRqq[exp][0]) * sin(2*omega*i) - (aLqq[exp][3] + aRqq[exp][3]) * cos(2*omega*i))
                    + dmunu * 2 * ((aLqq[exp][0] - aRqq[exp][0]) * sin(2*omega*i) - (aLqq[exp][3] - aRqq[exp][3]) * cos(2*omega*i));
            else if(wilson==2)
                foo[i] = cmunu * 2 * ((aLqq[exp][4] + aRqq[exp][4]) * cos(omega*i) + (aLqq[exp][5] + aRqq[exp][5]) * sin(omega*i))
                    + dmunu * 2 * ((aLqq[exp][4] - aRqq[exp][4]) * cos(omega*i) + (aLqq[exp][5] - aRqq[exp][5]) * sin(omega*i));
            else if(wilson==3)
                foo[i] = cmunu * 2 * ((aLqq[exp][4] + aRqq[exp][4]) * sin(omega*i) - (aLqq[exp][5] + aRqq[exp][5]) * cos(omega*i))
                    + dmunu * 2 * ((aLqq[exp][4] - aRqq[exp][4]) * sin(omega*i) - (aLqq[exp][5] - aRqq[exp][5]) * cos(omega*i));
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else
            std::cout<<"error with cmunu coeff with f function"<<std::endl;
    }
    return foo;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
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
//--------------------------------------------------------------------- Article -----------------------------------------------------------------------------//
std::vector<double> KineticAnalyze::calculateggg(int n, int exp,int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu){
    std::vector<double> foo(n);
    for(int i=0; i<n; i++){
        if(cmunu==0 and dmunu==0){
            if(wilson==0)
                foo[i] = cmunuL * 2 * (bLgg[exp][0] * cos(omega*i) + bLgg[exp][1] * sin(omega*i)) + cmunuR * 2 * (bRgg[exp][0] * cos(omega*i) + bRgg[exp][1] * sin(omega*i));
            else if(wilson==1)
                foo[i] = cmunuL * 2 * (bLgg[exp][0] * sin(omega*i) - bLgg[exp][1] * cos(omega*i)) + cmunuR * 2 * (bRgg[exp][0] * sin(omega*i) - bRgg[exp][1] * cos(omega*i));
            else if(wilson==2)
                foo[i] = cmunuL * 2 * bLgg[exp][2] + cmunuR * 2 * bRgg[exp][2];
            else if(wilson==3)
                foo[i] = cmunuL * 2 * bLgg[exp][3] + cmunuR * 2 * bRgg[exp][3];
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else if(cmunuL==0 and cmunuR==0){
            if(wilson==0)
                foo[i] = 2 * cmunu * ((bLgg[exp][0] * cos(omega*i) + bLgg[exp][1] * sin(omega*i)) + (bRgg[exp][0] * cos(omega*i) + bRgg[exp][1] * sin(omega*i)))
                + 2 * dmunu * ((bLgg[exp][0] * cos(omega*i) + bLgg[exp][1] * sin(omega*i)) - (bRgg[exp][0] * cos(omega*i) + bRgg[exp][1] * sin(omega*i)));
            else if(wilson==1)
                foo[i] = 2 * cmunu * ((bLgg[exp][0] * sin(omega*i) - bLgg[exp][1] * cos(omega*i)) + (bRgg[exp][0] * sin(omega*i) - bRgg[exp][1] * cos(omega*i)))
                + 2 * dmunu * ((bLgg[exp][0] * sin(omega*i) - bLgg[exp][1] * cos(omega*i)) - (bRgg[exp][0] * sin(omega*i) - bRgg[exp][1] * cos(omega*i)));
            else if(wilson==2)
                foo[i] = 2 * cmunu * (bLgg[exp][2] + bRgg[exp][2]) + 2 * dmunu * (bLgg[exp][2] - bRgg[exp][2]);
            else if(wilson==3)
                foo[i] = 2 * cmunu * (bLgg[exp][3] + bRgg[exp][3]) + 2 * dmunu * (bLgg[exp][3] - bRgg[exp][3]);
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else
            std::cout<<"error with cmunu coeff with f function"<<std::endl;
    }
    return foo;
}
std::vector<double> KineticAnalyze::calculategqq(int n, int exp,int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu){
    std::vector<double> foo(n);
    for(int i=0; i<n; i++){
        if(cmunu==0 and dmunu==0){
            if(wilson==0)
                foo[i] = cmunuL * 2 * (bLqq[exp][0] * cos(omega*i) + bLqq[exp][1] * sin(omega*i)) + cmunuR * 2 * (bRqq[exp][0] * cos(omega*i) + bRqq[exp][1] * sin(omega*i));
            else if(wilson==1)
                foo[i] = cmunuL * 2 * (bLqq[exp][0] * sin(omega*i) - bLqq[exp][1] * cos(omega*i)) + cmunuR * 2 * (bRqq[exp][0] * sin(omega*i) - bRqq[exp][1] * cos(omega*i));
            else if(wilson==2)
                foo[i] = cmunuL * 2 * bLqq[exp][2] + cmunuR * 2 * bRqq[exp][2];
            else if(wilson==3)
                foo[i] = cmunuL * 2 * bLqq[exp][3] + cmunuR * 2 * bRqq[exp][3];
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else if(cmunuL==0 and cmunuR==0){
            if(wilson==0)
                foo[i] = 2 * cmunu * ((bLqq[exp][0] * cos(omega*i) + bLqq[exp][1] * sin(omega*i)) + (bRqq[exp][0] * cos(omega*i) + bRqq[exp][1] * sin(omega*i)))
                + 2 * dmunu * ((bLqq[exp][0] * cos(omega*i) + bLqq[exp][1] * sin(omega*i)) - (bRqq[exp][0] * cos(omega*i) + bRqq[exp][1] * sin(omega*i)));
            else if(wilson==1)
                foo[i] = 2 * cmunu * ((bLqq[exp][0] * sin(omega*i) - bLqq[exp][1] * cos(omega*i)) + (bRqq[exp][0] * sin(omega*i) - bRqq[exp][1] * cos(omega*i)))
                + 2 * dmunu * ((bLqq[exp][0] * sin(omega*i) - bLqq[exp][1] * cos(omega*i)) - (bRqq[exp][0] * sin(omega*i) - bRqq[exp][1] * cos(omega*i)));
            else if(wilson==2)
                foo[i] = 2 * cmunu * (bLqq[exp][2] + bRqq[exp][2]) + 2 * dmunu * (bLqq[exp][2] - bRqq[exp][2]);
            else if(wilson==3)
                foo[i] = 2 * cmunu * (bLqq[exp][3] + bRqq[exp][3]) + 2 * dmunu * (bLqq[exp][3] - bRqq[exp][3]);
            else
                std::cout<<"error f fonction"<<std::endl;
        }
        else
            std::cout<<"error with cmunu coeff with f function"<<std::endl;
    }
    return foo;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------//
double KineticAnalyze::calculateMax(std::vector<double> vector_user){
    double foo = 0;
    for(int i=0; i<time24; i++){
        if(vector_user[i] >= foo)
            foo = vector_user[i];
    }
    if(foo==0)
        for(int i=0; i<time24; i++){
            if(vector_user[i] <= foo)
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
    w->SetFillColor(kWhite);
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
    TString format[2];    format[0]=".png";    format[1]=".eps";
    for(int i=0; i<2; i++){
        if(munu==0){
            h1->SetTitle("c_{L#mu#nu} = "+TString::Format("%f",cmunuL_class));
            w->SaveAs("results/fComparaisonL"+format[i]);
        }
        else if(munu==1){
            h1->SetTitle("c_{R#mu#nu} = "+TString::Format("%f",cmunuR_class));
            w->SaveAs("results/fComparaisonR"+format[i]);
        }
        else if(munu==2){
            h1->SetTitle("c_{#mu#nu} = "+TString::Format("%f",cmunu_class));
            w->SaveAs("results/fComparaisonC"+format[i]);
        }
        else if(munu==3){
            h1->SetTitle("d_{#mu#nu} = "+TString::Format("%f",dmunu_class));
            w->SaveAs("results/fComparaisonD"+format[i]);
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
    }
}

void KineticAnalyze::gTime(int munu, int exp){
    // histograms of fSME in time for each coefficients cmunu, Wilson is the choice of cmunu left/right or cmunu/dmunu
    TCanvas* w = new TCanvas("","",200,10,800,600);
    w->SetFillColor(kWhite);
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
    TString format[2];    format[0]=".png";    format[1]=".eps";
    for(int i=0; i<2; i++){
        if(munu==0){
            h1->SetTitle("c_{L#mu#nu} = "+TString::Format("%f",cmunuL_class));
            w->SaveAs("results/gComparaisonL"+format[i]);
        }
        else if(munu==1){
            h1->SetTitle("c_{R#mu#nu} = "+TString::Format("%f",cmunuR_class));
            w->SaveAs("results/gComparaisonR"+format[i]);
        }
        else if(munu==2){
            h1->SetTitle("c_{#mu#nu} = "+TString::Format("%f",cmunu_class));
            w->SaveAs("results/gComparaisonC"+format[i]);
        }
        else if(munu==3){
            h1->SetTitle("d_{#mu#nu} = "+TString::Format("%f",dmunu_class));
            w->SaveAs("results/gComparaisonD"+format[i]);
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
    }
}

void KineticAnalyze::gTimeTT(int munu, int exp){
    // histograms of fSME in time for each coefficients cmunu, Wilson is the choice of cmunu left/right or cmunu/dmunu
    TCanvas* w = new TCanvas("","",200,10,800,600);
    w->SetFillColor(kWhite);
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
    TString format[2];    format[0]=".png";    format[1]=".eps";
    for(int i=0; i<2; i++){
        if(munu==0){
            h1->SetTitle("c_{L#mu#nu} = "+TString::Format("%f",cmunuL_class));
            w->SaveAs("results/gComparaisonLTT"+format[i]);
        }
        else if(munu==1){
            h1->SetTitle("c_{R#mu#nu} = "+TString::Format("%f",cmunuR_class));
            w->SaveAs("results/gComparaisonRTT"+format[i]);
        }
        else if(munu==2){
            h1->SetTitle("c_{#mu#nu} = "+TString::Format("%f",cmunu_class));
            w->SaveAs("results/gComparaisonCTT"+format[i]);
        }
        else if(munu==3){
            h1->SetTitle("d_{#mu#nu} = "+TString::Format("%f",dmunu_class));
            w->SaveAs("results/gComparaisonDTT"+format[i]);
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
    }
}

void KineticAnalyze::amplEnergyComparaison(bool isBenchmark){
    TString name[4];    TString name2[4];
    if(isBenchmark){
        name[0]="cL";    name[1]="cR";    name[2]="c";    name[3]="d";
        name2[0]="_{XX}";    name2[1]="_{XY}";    name2[2]="_{XZ}";    name2[3]="_{YZ}";
    }    
    else{
        name[0]="cL";    name[1]="cR";    name[2]="c";    name[3]="d";
        name2[0]="_{TX}";    name2[1]="_{TY}";    name2[2]="_{TZ}";    name2[3]="_{TT}";
    }

    std::vector<TCanvas*> c(4);
    std::vector<TMultiGraph*> m(4);
    std::vector<TLegend*> legend(4);
    std::vector< std::vector<TGraph*> > g(4);
    std::vector< std::vector<double*> > y(4);
    for(int i=0; i<4; i++){
        m[i] = new TMultiGraph();
        y[i] = std::vector<double*>(4);
        g[i] = std::vector<TGraph*>(4);
        for(int j=0; j<4; j++){
            y[i][j] = new double[nExp];
        }
    }

    double* x = new double[nExp];
    for(int i=0; i<nExp; i++)
        x[i] = expList[i];

    for(int j=0; j<4; j++){
        for(int k=0; k<nExp; k++){
            if(isBenchmark){
                y[0][j][k] = amplitudeLListf[k][j];
                y[1][j][k] = amplitudeRListf[k][j];
                y[2][j][k] = amplitudeCListf[k][j];
                y[3][j][k] = amplitudeDListf[k][j];
            }
            else{
                y[0][j][k] = amplitudeLListg[k][j];
                y[1][j][k] = amplitudeRListg[k][j];
                y[2][j][k] = amplitudeCListg[k][j];
                y[3][j][k] = amplitudeDListg[k][j];                
            }
        }       
    }

    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            g[i][j] = new TGraph(nExp, x, y[i][j]);
            g[i][j]->SetMarkerStyle(20);
            g[i][j]->SetMarkerSize(1);
        }
        g[i][0]->SetMarkerColor(i+1);
        g[i][1]->SetMarkerColor(i+1);
        g[i][2]->SetMarkerColor(i+1);
        g[i][3]->SetMarkerColor(i+1);
        for(int j=0; j<4; j++){
            m[j]->Add(g[i][j]);
            legend[j] = new TLegend(0.1,0.7,0.48,0.9);
            legend[j]->SetHeader("Amplitude","C"); // option "C" allows to center the header
            legend[j]->AddEntry(g[0][j],"c_{L}"); // option "l" is for line (form of legend)
            legend[j]->AddEntry(g[1][j],"c_{R}");
            legend[j]->AddEntry(g[2][j],"c");
            legend[j]->AddEntry(g[3][j],"d");
        }
    }
    for(int i=0; i<4; i++){
        c[i] = new TCanvas(name[i],"",200,10,800,600);
        c[i]->SetFillColor(kWhite);
        m[i]->Draw("AP");
        m[i]->SetTitle("");
        legend[i]->Draw();
        c[i]->SetLogx();
        c[i]->SaveAs("results/amplitude/amplEnergy"+name2[i]+".png");
        c[i]->SaveAs("results/amplitude/amplEnergy"+name2[i]+".eps");
    }
}

void KineticAnalyze::compareFusAni(int munu, int exp){
    // histograms of fSME in time for each coefficients cmunu, Wilson is the choice of cmunu left/right or cmunu/dmunu
    TCanvas* w = new TCanvas("","",200,10,800,600);
    w->SetFillColor(kWhite);
    TH1F* h1 = new TH1F("","", time24, 0, 24);
    TH1F* h2 = new TH1F("","", time24, 0, 24);
    for (int k=0; k<86400; k++){
        if(munu==0){
            h1->SetBinContent(k+1, fLListqq[exp][0][k]);
            h2->SetBinContent(k+1, fLListgg[exp][0][k]);
        }
        else if(munu==1){
            h2->SetBinContent(k+1, fRListqq[exp][0][k]);
            h1->SetBinContent(k+1, fRListgg[exp][0][k]);
        }
        else if(munu==2){
            h1->SetBinContent(k+1, fCListqq[exp][0][k]);
            h2->SetBinContent(k+1, fCListgg[exp][0][k]);
        }
        else if(munu==3){
            h2->SetBinContent(k+1, fDListqq[exp][0][k]);
            h1->SetBinContent(k+1, fDListgg[exp][0][k]);
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
    }

    h1->Draw();          h1->Write();   h1->SetLineWidth(2);   h1->SetLineColor(kRed);
    h2->Draw("SAME");    h2->Write();   h2->SetLineWidth(2);   h2->SetLineColor(kBlue);

    if(munu==1 or munu==3){
        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);
    }
    else{
        h1->SetLineColor(kRed);
        h2->SetLineColor(kBlue);
    }

    h1->GetYaxis()->SetTitle("f_{SME}(t)");
    h1->GetXaxis()->SetTitle("sideral time #hat{t} (in h)");
    h1->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    if(munu==1 or munu==3){
        legend->AddEntry(h2,"quark/antiquark anihilation only","l"); // option "l" is for line (form of legend)
        legend->AddEntry(h1,"gluon fusion only","l");
    }
    else{
        legend->AddEntry(h1,"quark/antiquark anihilation only","l"); // option "l" is for line (form of legend)
        legend->AddEntry(h2,"gluon fusion only","l");
    }

    w->Update();
    TString format[2];    format[0]=".png";    format[1]=".eps";
    for(int i=0; i<2; i++){
        if(munu==0){            
            legend->SetHeader("f_{SME}(t) functions histogram for c_{LXX} = - c_{LYY} = 0.01", "C"); // option "C" allows to center the header
            legend->Draw();
            w->SaveAs("results/compaFusAni/fComparaisonL"+format[i]);
        }
        else if(munu==1){
            legend->SetHeader("f_{SME}(t) functions histogram for c_{RXX} = - c_{RYY} = 0.01","C"); // option "C" allows to center the header
            legend->Draw();
            w->SaveAs("results/compaFusAni/fComparaisonR"+format[i]);
        }
        else if(munu==2){
            legend->SetHeader("f_{SME}(t) functions histogram for c_{XX} = - c_{YY} = 0.01","C"); // option "C" allows to center the header
            legend->Draw();
            w->SaveAs("results/compaFusAni/fComparaisonC"+format[i]);
        }
        else if(munu==3){
            legend->SetHeader("f_{SME}(t) functions histogram for d_{XX} = - d_{YY} = 0.01","C"); // option "C" allows to center the header
            legend->Draw();
            w->SaveAs("results/compaFusAni/fComparaisonD"+format[i]);
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
    }

}

void KineticAnalyze::compareCMSD0(int munu){
    TCanvas* w = new TCanvas("","",200,10,800,600);
    w->SetFillColor(kWhite);
    TH1F* h1 = new TH1F("","", 4, 0, 4);
    TH1F* h2 = new TH1F("","", 4, 0, 4);
    std::vector<std::vector<TMarker*> > m(2);
    for(int i=0; i<2; i++)
        m[i] = std::vector<TMarker*>(4);


    h1->GetXaxis()->SetBinLabel(1,"c_{XX} = -c_{YY} #neq 0"); 
    h1->GetXaxis()->SetBinLabel(2,"c_{XY} = c_{YX} #neq 0"); 
    h1->GetXaxis()->SetBinLabel(3,"c_{XZ} = -c_{ZX} #neq 0"); 
    h1->GetXaxis()->SetBinLabel(4,"c_{YZ} = -c_{ZY} #neq 0"); 

    for (int k=0; k<4; k++){
        if(munu==0){
            h1->SetBinContent(k+1, amplitudeLListf[5][k]);
            h1->Draw("P SAME");
            m[0][k] = new TMarker(k+0.5, amplitudeLListf[5][k], 20);
            m[0][k]->SetMarkerColor(kBlue);
            m[0][k]->SetMarkerSize(1);
            m[0][k]->Draw("SAME");
            h2->SetBinContent(k+1, amplitudeLListfTEV[5][k]);
            h2->Draw("P SAME");
            m[1][k] = new TMarker(k+0.5, amplitudeLListfTEV[5][k], 20);
            m[1][k]->SetMarkerColor(kRed);
            m[1][k]->SetMarkerSize(1);
            m[1][k]->Draw("SAME");
        }
        else if(munu==1){
            h1->SetBinContent(k+1, amplitudeRListf[5][k]);
            h1->Draw("P SAME");
            m[0][k] = new TMarker(k+0.5, amplitudeRListf[5][k], 20);
            m[0][k]->SetMarkerColor(kBlue);
            m[0][k]->SetMarkerSize(1);
            m[0][k]->Draw("SAME");
            h2->SetBinContent(k+1, amplitudeRListfTEV[5][k]);
            h2->Draw("P SAME");
            m[1][k] = new TMarker(k+0.5, amplitudeRListfTEV[5][k], 20);
            m[1][k]->SetMarkerColor(kRed);
            m[1][k]->SetMarkerSize(1);
            m[1][k]->Draw("SAME");
        }
        else if(munu==2){
            h1->SetBinContent(k+1, amplitudeCListf[5][k]);
            h1->Draw("P SAME");
            m[0][k] = new TMarker(k+0.5, amplitudeCListf[5][k], 20);
            m[0][k]->SetMarkerColor(kBlue);
            m[0][k]->SetMarkerSize(1);
            m[0][k]->Draw("SAME");
            h2->SetBinContent(k+1, amplitudeCListfTEV[5][k]);
            h2->Draw("P SAME");
            m[1][k] = new TMarker(k+0.5, amplitudeCListfTEV[5][k], 20);
            m[1][k]->SetMarkerColor(kRed);
            m[1][k]->SetMarkerSize(1);
            m[1][k]->Draw("SAME");
        }
        else if(munu==3){
            h1->SetBinContent(k+1, amplitudeDListf[5][k]);
            h1->Draw("P SAME");
            m[0][k] = new TMarker(k+0.5, amplitudeDListf[5][k], 20);
            m[0][k]->SetMarkerColor(kBlue);
            m[0][k]->SetMarkerSize(1);
            m[0][k]->Draw("SAME");
            h2->SetBinContent(k+1, amplitudeDListfTEV[5][k]);
            h2->Draw("P SAME");
            m[1][k] = new TMarker(k+0.5, amplitudeDListfTEV[5][k], 20);
            m[1][k]->SetMarkerColor(kRed);
            m[1][k]->SetMarkerSize(1);
            m[1][k]->Draw("SAME");
        }
        else
            std::cout<<"error with cmunu"<<std::endl;
   }
         h1->Write();   h1->SetLineWidth(2);    h1->SetFillColor(kBlue);
         h2->Write();   h2->SetLineWidth(2);    h2->SetFillColor(kRed);

    h1->GetYaxis()->SetTitle("Amplitude f(t)");
    h1->GetXaxis()->SetTitle("");
    h1->SetMaximum(0.02);
    if(munu==2)
        h1->SetMaximum(0.006);
    if(munu==3)
        h1->SetMaximum(0.035);
    h1->SetMinimum(0);
    h1->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Experiment","C"); // option "C" allows to center the header
    legend->AddEntry(m[0][0],"CMS 1.96 TEV","p"); // option "l" is for line (form of legend)
    legend->AddEntry(m[1][0],"D0 1.96 TEV","p");
    legend->Draw();

    w->Update();
    TString format[2];    format[0]=".png";    format[1]=".eps";
    for(int i=0; i<2; i++){
        if(munu==0)
            w->SaveAs("results/compaFusAni/compaCMSTEVL"+format[i]);
        else if(munu==1)
            w->SaveAs("results/compaFusAni/compaCMSTEVR"+format[i]);
        else if(munu==2)
            w->SaveAs("results/compaFusAni/compaCMSTEVC"+format[i]);
        else if(munu==3)
            w->SaveAs("results/compaFusAni/compaCMSTEVD"+format[i]);
        else
            std::cout<<"error with cmunu"<<std::endl;
    }

}

void KineticAnalyze::earthSignal(bool isXX){

    double pas = 1000;
    double cmunuEarthSignal = 0.1;
    double tmp = M_PI/pas;                  // latitude
    double tmp2 = 2*(M_PI/pas);             // azimuth
    double a0, a1, a2, a3, a4, a5, arg1, arg2;
    TCanvas* c = new TCanvas("max f_{SME}(#lambda, #theta) f","", 10,10,800,600);
    TH2F* h = new TH2F("", "", pas, 0, M_PI, pas, 0, 2*M_PI);

    for(int i = 0; i<pas; i++)
        for(int j=0; j<pas; j++){
            if(isXX){
                a1 = calculateCoefficent_a(1, Areal[3], i*tmp, j*tmp2, true);
                a2 = calculateCoefficent_a(2, Areal[3], i*tmp, j*tmp2, true);
                a0 = (a1-a2)/2.;
                a3 = calculateCoefficent_a(3, Areal[3], i*tmp, j*tmp2, true);
                if(a0*cos(2*omega*arg1+0.001) + a3*sin(2*omega*arg1+0.001)<0){
                    if(a0>0)
                        arg1 = atan(a3/a0);
                    else if(a0<0)
                        arg1 = -atan(a3/abs(a0)) + M_PI;
                    else if(a0==0)
                        arg1 = 0;
                    h->SetBinContent(i+1, j+1, 2*cmunuEarthSignal*(a0*cos(2*omega*arg1) + a3*sin(2*omega*arg1)) );
                }
            }
            else if(!isXX){
                a4 = calculateCoefficent_a(4, Areal[3], i*tmp, j*tmp2, true);
                a5 = calculateCoefficent_a(5, Areal[3], i*tmp, j*tmp2, true);
                if(a4*cos(omega*arg1+0.001) + a5*sin(omega*arg1+0.001)<0){
                    if(a4>0)
                        arg1 = atan(a5/a4);
                    else if(a4<0)
                        arg1 = -atan(a5/abs(a4)) + M_PI;
                    else if(a4==0)
                        arg1 = 0;
                    h->SetBinContent(i+1, j+1, 2*cmunuEarthSignal*(a4*cos(omega*arg1) + a5*sin(omega*arg1)) );
                }
            }
        }

    if(isXX)
        h->SetTitle("Amplitude f_{SME}(#lambda, #theta)) fXX");
    else
        h->SetTitle("Amplitude f_{SME}(#lambda, #theta)) fXZ");
    h->GetYaxis()->SetTitle("Azimuth #theta (in rad)");
    h->GetXaxis()->SetTitle("Latitude #lambda (in rad)");

    TMarker* pointCMS = new TMarker(latitudeCMS, azimuthCMS, 41);
    TMarker* pointD0 = new TMarker(latitD0, azimD0Redifine, 41);
    TMarker* pointATLAS = new TMarker(latitudeATLAS, azimuthATLAS, 41);
    pointCMS->SetMarkerSize(3);
    pointCMS->SetMarkerColor(kRed);
    pointD0->SetMarkerSize(3);
    pointD0->SetMarkerColor(kGreen);
    pointATLAS->SetMarkerSize(3);
    pointATLAS->SetMarkerColor(kMagenta);

    h->Write();
    h->Draw("colz");
    pointCMS->Draw("same");
    pointD0->Draw("same");
    pointATLAS->Draw("same");
    h->SetStats(0);

    TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);
    legend->AddEntry(pointCMS,"CMS","p");
    legend->AddEntry(pointATLAS,"ATLAS","p");
    legend->AddEntry(pointD0,"D0","p");
    legend->Draw();

    if(isXX){
        c->SaveAs("results/earthSignal/XX.png");
        c->SaveAs("results/earthSignal/XX.eps");
    }
    else{
        c->SaveAs("results/earthSignal/XZ.png");
        c->SaveAs("results/earthSignal/XZ.eps");
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

TH1F* KineticAnalyze::statHistosf(TString name, int exp, TString wilson, bool isXX, double ttbar){
    int bin=24;
    double s12 = ttbar/((double)bin);
    TH1F* foo = new TH1F(name, name, bin, 0, bin);
    if(wilson=="L")//L
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, s12*fLList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, s12*fLList[exp][2][i*3600]);
        }
    else if(wilson=="R")//R
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, s12*fRList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, s12*fRList[exp][2][i*3600]);
        }
    else if(wilson=="C")//C
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, s12*fCList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, s12*fCList[exp][2][i*3600]);
        }
    else if(wilson=="D")//D
        for(int i=0; i<bin; i++){
            if(isXX)
                foo->SetBinContent(i+1, s12*fDList[exp][0][i*3600]);
            else
                foo->SetBinContent(i+1, s12*fDList[exp][2][i*3600]);
        }
    else
        std::cout<<"error with stats"<<std::endl;

    return foo;
}


TH1F* KineticAnalyze::statHistosg(TString name, int exp, TString wilson, bool isTZ, double ttbar){
    int bin=24;
    double s12 = ttbar/((double)bin);
    TH1F* foo = new TH1F(name, name, bin, 0, bin);
    if(wilson=="L")//L
        for(int i=0; i<bin; i++){
            if(!isTZ)
                foo->SetBinContent(i+1, (s12*gLList[exp][0][i*3600]));
            else
                foo->SetBinContent(i+1, (s12*gLList[exp][2][i*3600]));
        }
    else if(wilson=="R")//R
        for(int i=0; i<bin; i++){
            if(!isTZ)
                foo->SetBinContent(i+1, (s12*gRList[exp][0][i*3600]));
            else
                foo->SetBinContent(i+1, (s12*gRList[exp][2][i*3600]));
        }
    else if(wilson=="C")//C
        for(int i=0; i<bin; i++){
            if(!isTZ)
                foo->SetBinContent(i+1, (s12*gCList[exp][0][i*3600]));
            else
                foo->SetBinContent(i+1, (s12*gCList[exp][2][i*3600]));
        }
    else if(wilson=="D")//D
        for(int i=0; i<bin; i++){
            if(!isTZ)
                foo->SetBinContent(i+1, (s12*gDList[exp][0][i*3600]));
            else
                foo->SetBinContent(i+1, (s12*gDList[exp][2][i*3600]));
        }
    else
        std::cout<<"error with stats"<<std::endl;

    return foo;
}



TH1F* KineticAnalyze::statHistosgTT(TString name, int exp, TString wilson, double ttbar){
    int bin=24;
    double s12 = ttbar/((double)bin);
    TH1F* foo = new TH1F(name, name, bin, 0, bin);

    if(wilson=="L")//L
        for(int i=0; i<bin; i++)
            foo->SetBinContent(i+1, (s12*gLList[exp][3][i*3600]));
    else if(wilson=="R")//R
        for(int i=0; i<bin; i++)
            foo->SetBinContent(i+1, (s12*gRList[exp][3][i*3600]));
    else if(wilson=="C")//C
        for(int i=0; i<bin; i++)
            foo->SetBinContent(i+1, (s12*gCList[exp][3][i*3600]));
    else if(wilson=="D")//D
        for(int i=0; i<bin; i++)
            foo->SetBinContent(i+1, (s12*gDList[exp][3][i*3600]));
    else
        std::cout<<"error with stats"<<std::endl;
    return foo;
}

