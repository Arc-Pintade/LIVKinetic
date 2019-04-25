#ifndef KineticAnalyze_h
#define KineticAnalyze_h

#include "../include/const.hpp"
#include <TMatrixD.h>
#include <TString.h>
#include <TH1F.h>

class KineticAnalyze{

    private:


//______________ cmunu _____________________//

        double cmunuL_class;
        double cmunuR_class;
        double cmunu_class;
        double dmunu_class;

        int nExp;
        int time24;
        std::vector<double> expList;

        std::vector<int> nF;
        std::vector<int> nPqq;
        std::vector<int> nPgg;
        std::vector<TMatrixD> AF;
        std::vector<TMatrixD> APqq;
        std::vector<TMatrixD> APgg;

        std::vector<TMatrixD> Aani;
        std::vector<TMatrixD> Afus;
        std::vector<TMatrixD> Aprod;
        std::vector<TMatrixD> Areal;

        std::vector< std::vector<double> > aL;
        std::vector< std::vector<double> > aR;
        std::vector< std::vector<double> > bL;
        std::vector< std::vector<double> > bR;

        std::vector< std::vector< std::vector<double> > > fLList;
        std::vector< std::vector< std::vector<double> > > fRList;
        std::vector< std::vector< std::vector<double> > > fCList;
        std::vector< std::vector< std::vector<double> > > fDList;
        std::vector< std::vector< std::vector<double> > > gLList;
        std::vector< std::vector< std::vector<double> > > gRList;
        std::vector< std::vector< std::vector<double> > > gCList;
        std::vector< std::vector< std::vector<double> > > gDList;

        std::vector< std::vector<double> > amplitudeLListf;
        std::vector< std::vector<double> > amplitudeRListf;
        std::vector< std::vector<double> > amplitudeCListf;
        std::vector< std::vector<double> > amplitudeDListf;
        std::vector< std::vector<double> > amplitudeLListg;
        std::vector< std::vector<double> > amplitudeRListg;
        std::vector< std::vector<double> > amplitudeCListg;
        std::vector< std::vector<double> > amplitudeDListg;

        TMatrixD calculateAverageMatrix(int option, int nP1_user, int nP2_user, int nF_user, TMatrixD AP1_user, TMatrixD AP2_user, TMatrixD AF_user);
        double calculateCoefficent_a(int i, TMatrixD m_user, double latitude_user, double azimuth_user, bool isD0);
        double calculateCoefficent_b(int i, TMatrixD m_user, double latitude_user, double azimuth_user);
        std::vector<double> calculatef(int n, int exp, int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu);
        std::vector<double> calculateg(int n, int exp, int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu);

        double calculateMax(std::vector<double> vector_user);

//--------------------------------------------------------------- Article ---------------------------------------------------------------//
        std::vector< std::vector<double> > aLTEV;
        std::vector< std::vector<double> > aRTEV;
        std::vector< std::vector<double> > aLgg;
        std::vector< std::vector<double> > aRgg;
        std::vector< std::vector<double> > bLgg;
        std::vector< std::vector<double> > bRgg;
        std::vector< std::vector<double> > aLqq;
        std::vector< std::vector<double> > aRqq;
        std::vector< std::vector<double> > bLqq;
        std::vector< std::vector<double> > bRqq;

        std::vector< std::vector< std::vector<double> > > fLListTEV;
        std::vector< std::vector< std::vector<double> > > fRListTEV;
        std::vector< std::vector< std::vector<double> > > fCListTEV;
        std::vector< std::vector< std::vector<double> > > fDListTEV;
        std::vector< std::vector< std::vector<double> > > fLListgg;
        std::vector< std::vector< std::vector<double> > > fRListgg;
        std::vector< std::vector< std::vector<double> > > fCListgg;
        std::vector< std::vector< std::vector<double> > > fDListgg;
        std::vector< std::vector< std::vector<double> > > gLListgg;
        std::vector< std::vector< std::vector<double> > > gRListgg;
        std::vector< std::vector< std::vector<double> > > gCListgg;
        std::vector< std::vector< std::vector<double> > > gDListgg;
        std::vector< std::vector< std::vector<double> > > fLListqq;
        std::vector< std::vector< std::vector<double> > > fRListqq;
        std::vector< std::vector< std::vector<double> > > fCListqq;
        std::vector< std::vector< std::vector<double> > > fDListqq;
        std::vector< std::vector< std::vector<double> > > gLListqq;
        std::vector< std::vector< std::vector<double> > > gRListqq;
        std::vector< std::vector< std::vector<double> > > gCListqq;
        std::vector< std::vector< std::vector<double> > > gDListqq;

        std::vector< std::vector<double> > amplitudeLListfTEV;
        std::vector< std::vector<double> > amplitudeRListfTEV;
        std::vector< std::vector<double> > amplitudeCListfTEV;
        std::vector< std::vector<double> > amplitudeDListfTEV;

        std::vector<double> calculatefTEV(int n, int exp, int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu);
        std::vector<double> calculatefgg(int n, int exp, int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu);
        std::vector<double> calculateggg(int n, int exp, int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu);
        std::vector<double> calculatefqq(int n, int exp, int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu);
        std::vector<double> calculategqq(int n, int exp, int wilson, double cmunuL, double cmunuR, double cmunu, double dmunu);
//---------------------------------------------------------------------------------------------------------------------------------------//

    public:

//_________________________________________________//
//__________________ Constructor __________________//
//_________________________________________________//

        KineticAnalyze();
        KineticAnalyze(double cmunuL_user, double cmunuR_user, double cmunu_user, double dmunu_user);
        ~KineticAnalyze(){};
//_________________________________________________//
//__________________ Read Values __________________//
//_________________________________________________//

        // Text file input, text file with a number of events and the value of the average matrix assossiated
        static void showMatrix(TMatrixD m);
        static int readNumber(TString s);
        static TMatrixD readMatrix(TString s);
        static std::vector<TMatrixD> readVecMatrix(TString s, int n);
        static double* extractArray(TString s, int n);

//_________________________________________________//
//________________ Histos and graphs ______________//
//_________________________________________________//

        void fTime(int munu, int exp);
        void gTime(int munu, int exp);
        void gTimeTT(int munu, int exp);
        void amplEnergy(bool isBenchmark);
        void amplEnergyComparaison(bool isBenchmark);
        void compareFusAni(int munu, int exp);
        void compareCMSD0(int munu);
        void earthSignal(TString XX);
        void amunuHist();
        void amunuHistSolo(int row, int column, bool isComparaison);

//_________________________________________________//
//________________ statistical stuff ______________//
//_________________________________________________//

        TH1F* statHistosConst(TString name, double value);
        TH1F* statHistosf(TString name, int exp, TString wilson, bool isXX, double ttbar);
        TH1F* statHistosg(TString name, int exp, TString wilson, bool isTZ, double ttbar);
        TH1F* statHistosgTT(TString name, int exp, TString wilson, double ttbar);

};

#endif
