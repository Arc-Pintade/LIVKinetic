#ifndef KineticAnalyze_h
#define KineticAnalyze_h

#include "const.hpp"
#include <TMatrixD.h>
#include <TString.h>

class KineticAnalyze{

    private:

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

    public:

//_________________________________________________//
//__________________ Constructor __________________//
//_________________________________________________//

        KineticAnalyze();
        ~KineticAnalyze(){};
//_________________________________________________//
//__________________ Read Values __________________//
//_________________________________________________//

        // Text file input, text file with a number of events and the value of the average matrix assossiated
        static int readNumber(TString s);
        static TMatrixD readMatrix(TString s);
        static double* extractArray(TString s, int n);

//_________________________________________________//
//________________ Histos and graphs ______________//
//_________________________________________________//

        void fTime(int munu, int exp);
        void gTime(int munu, int exp);
        void gTimeTT(int munu, int exp);
        void amplEnergy();

};

#endif

