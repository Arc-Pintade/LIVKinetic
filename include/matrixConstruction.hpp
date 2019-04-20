////////////////////////////////////////////////////////
//                                                    //
// Program Created by Aurélien CARLE (CMS Experiment) //
//    Thanks Xavier VALCARCE and Gregoire UHLRICH     //
//                                                    //
////////////////////////////////////////////////////////

#ifndef matrixElement_h
#define matrixElement_h
#define matrixElement_cxx

#include <TLorentzVector.h>
#include <TString.h>
#include <TMatrixD.h>

#ifndef __CINT__


//___________ Matrix Element Class SM and SME ____________//

class MatrixElement {

    private :

//********* Initiales Values **********//

        // P
        TLorentzVector pBeam1;
        TLorentzVector pBeam2;
        TLorentzVector pt;
        TLorentzVector ptbar;
        // F
        TLorentzVector pb;
        TLorentzVector pDecay1;
        TLorentzVector pAntiDecay1;
        TLorentzVector pbbar;
        TLorentzVector pDecay2;
        TLorentzVector pAntiDecay2;


//******** Calculated Values **********//

        bool isFusion;
        double time;

        // For SM
        double Pqqbar;       double calculatePqqbar();
        double P2g;          double calculateP2g();
        double F;            double calculateF();
        double Fbar;         double calculateFbar();
        // For SME
        TMatrixD Pqqbar_munu;
        TMatrixD P2g_munu;
        TMatrixD F_munu;
        TMatrixD Fbar_munu;
        TMatrixD calculatePqqbar_munu();
        TMatrixD calculateP2g_munu();
        TMatrixD calculateF_munu();
        TMatrixD calculateFbar_munu();

    public :

//************ Constructor *************//
    // Takes 4-momentum and creat MatrixElement, isFusion is bool which mark digluon fusion or quark/anti-quark anihilation

        MatrixElement();
        ~MatrixElement(){};

        MatrixElement(TLorentzVector pBeam1_user, TLorentzVector pBeam2_user, TLorentzVector pt_user, TLorentzVector ptbar_user, TLorentzVector pb_user, TLorentzVector pDecay1_user, TLorentzVector pAntiDecay1_user, TLorentzVector pbbar_user, TLorentzVector pDecay2_user, TLorentzVector pAntiDecay2_user, bool isFusion_user);

        MatrixElement(TLorentzVector pBeam1_user, TLorentzVector pBeam2_user, TLorentzVector pt_user, TLorentzVector ptbar_user, TLorentzVector pb_user, TLorentzVector pDecay1_user, TLorentzVector pAntiDecay1_user, TLorentzVector pbbar_user, TLorentzVector pDecay2_user, TLorentzVector pAntiDecay2_user, bool isFusion_user, double time_user);


//****** Recover values previously obtained ******//

        double getPqqbar();
        double getP2g();
        double getF();
        double getFbar();
        TMatrixD getPqqbar_munu();
        TMatrixD getP2g_munu();
        TMatrixD getF_munu();
        TMatrixD getFbar_munu();

//___________________ Some Tools ______________________//

        static double normMinkowski2(TLorentzVector p);
        static double scalarProduct(TLorentzVector p1, TLorentzVector p2);
        static void writeDownVecMatrix(double*** Vecmatrix, TString s, double nEvents);
        static void writeNumMatrix(TString s, int n, TMatrixD m);
        static void writeArray(TString s, int n, double* array);

    //******************* Matrix Tools ********************//

        static TMatrixD matrixA(TMatrixD M, double m);
        static TMatrixD matrixAF(TMatrixD AF, TMatrixD AFbar, double F, double Fbar);
        static TMatrixD matrixGen(TLorentzVector p1, TLorentzVector p2);
        static TMatrixD matrixGenvMss(TLorentzVector p1, TLorentzVector p2);
        static TMatrixD calculateMatrix(int ntot, int n, double*** vecMat);
        static void showMatrix(TMatrixD m);

};


#endif
#endif
