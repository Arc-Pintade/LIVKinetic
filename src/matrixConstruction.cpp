////////////////////////////////////////////////////////
//                                                    //
// Program Created by Aurélien CARLE (CMS Experiment) //
//    Thanks Xavier VALCARCE and Gregoire UHLRICH     //
//                                                    //
////////////////////////////////////////////////////////

#include "../include/matrixConstruction.hpp"

//___________________ MATRIX ELEMENT ___ SM ______________________//

//________________________ Constructor  _________________________//


MatrixElement::MatrixElement(TLorentzVector pBeam1_user, TLorentzVector pBeam2_user, TLorentzVector pt_user, TLorentzVector ptbar_user, TLorentzVector pb_user, TLorentzVector pDecay1_user, TLorentzVector pAntiDecay1_user, TLorentzVector pbbar_user, TLorentzVector pDecay2_user, TLorentzVector pAntiDecay2_user, bool isFusion_user){

    isFusion = isFusion_user;
    // Creation
    pBeam1 = pBeam1_user;
    pBeam2 = pBeam2_user;
    pt = pt_user;
    ptbar = ptbar_user;
    // Decay W⁺
    pb = pb_user;
    pDecay1 = pDecay1_user;
    pAntiDecay1 = pAntiDecay1_user;
    // Decay W⁻
    pbbar = pbbar_user;
    pDecay2 = pDecay2_user;
    pAntiDecay2 = pAntiDecay2_user;
    if(isFusion)
        P2g = calculateP2g();
    else
    Pqqbar = calculatePqqbar();
    F = calculateF();
    Fbar = calculateFbar();
}

//_________________ Functions Pqqbar P2g F Fbar __________________//

// Calcul Pqqbar
double MatrixElement::calculatePqqbar(){
    double s = normMinkowski2(pt + ptbar);

    double foo = ((16 * gS4) / (18 * s * s)) * (scalarProduct(pBeam1, pt) * scalarProduct(pBeam2, ptbar));
    foo += ((16 * gS4) / (18 * s * s)) * (scalarProduct(pBeam1, ptbar) * scalarProduct(pBeam2, pt));
    foo += ((16 * gS4) / (18 * s * s)) * (scalarProduct(pBeam1, pBeam2) * mt2);
    return foo;
}

// Calcul P2g
double MatrixElement::calculateP2g(){
    // Intermediate Element
    double s = normMinkowski2(pt + ptbar);
    double t = normMinkowski2(pBeam2 - pt);
    double u = normMinkowski2(pBeam1 - pt);
    double tmm = t - mt2;
    double tpm = t + mt2;
    double umm = u - mt2;
    double upm = u + mt2;
    double tmu = t - u;
    double sm4m = s - 4*mt2;

    double foo = ((3 * gS4)/4) * ((tmm * umm) / (s * s));
    foo += (gS4/6) * (((tmm * umm) - (2*mt2 * tpm)) / (tmm * tmm));
    foo += (gS4/6) * (((umm * tmm) - (2*mt2 * upm)) / (umm * umm));
    foo += ((3*gS4)/8) * (((tmm * umm) - (mt2 * tmu)) / (s * tmm));
    foo += ((3*gS4)/8) * (((umm * tmm) + (mt2 * tmu)) / (s * umm));
    foo += (-gS4/24) * ((mt2 * sm4m) / (tmm * umm));

    return foo;
}

// Calcul F
double MatrixElement::calculateF(){
    double foo = (-4*gW4) / (mt2 * gammat2 * ((pow((normMinkowski2(pDecay1 + pAntiDecay1) - mW2),2)) + (mW2 * gammaW2)));
    foo *= scalarProduct(pDecay1, pb) * scalarProduct(pAntiDecay1, pt);
    return foo;
}

// Calcul Fbar
double MatrixElement::calculateFbar(){
    double foo = (-4*gW4) / (mt2 * gammat2 * ((pow((normMinkowski2(pDecay2 + pAntiDecay2) - mW2),2)) + (mW2 * gammaW2)));
    foo *= scalarProduct(pDecay2, pbbar) * scalarProduct(pAntiDecay2, ptbar);
    return foo;
}

/***************************************************************/

//___________________ MATRIX ELEMENT ___ SME ______________________//

//__________________________ Constructor  ________________________//


MatrixElement::MatrixElement(TLorentzVector pBeam1_user, TLorentzVector pBeam2_user, TLorentzVector pt_user, TLorentzVector ptbar_user, TLorentzVector pb_user, TLorentzVector pDecay1_user, TLorentzVector pAntiDecay1_user, TLorentzVector pbbar_user, TLorentzVector pDecay2_user, TLorentzVector pAntiDecay2_user, bool isFusion_user, Double_t time_user) {

    P2g_munu.ResizeTo(4,4,-1);
    Pqqbar_munu.ResizeTo(4,4,-1);
    F_munu.ResizeTo(4,4,-1);
    Fbar_munu.ResizeTo(4,4,-1);

    isFusion = isFusion_user;
    pBeam1 = pBeam1_user;
    pBeam2 = pBeam2_user;
    // creation
    pt = pt_user;
    ptbar = ptbar_user;
    // decay W⁺
    pb = pb_user;
    pDecay1 = pDecay1_user;
    pAntiDecay1 = pAntiDecay1_user;
    // decay W⁻
    pbbar = pbbar_user;
    pDecay2 = pDecay2_user;
    pAntiDecay2 = pAntiDecay2_user;
    // time
    time = time_user;
    if(isFusion)
        P2g_munu = calculateP2g_munu();
    else
    Pqqbar_munu = calculatePqqbar_munu();
    F_munu = calculateF_munu();
    Fbar_munu = calculateFbar_munu();
}

//________________________ Pqqbar SME ____________________________//

TMatrixD MatrixElement::calculatePqqbar_munu(){
    // Intermediate Element
    double s = normMinkowski2(pt + ptbar);

    // Ppqqbar
    TMatrixD foo(4,4);
    foo += ((16 * gS4) / (18 * s * s)) * (scalarProduct(pBeam1, ptbar) * (matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T()));
    foo += ((16 * gS4) / (18 * s * s)) * (scalarProduct(pBeam2, pt) * (matrixGen(ptbar, pBeam1) + matrixGen(ptbar, pBeam1).T()));
    foo += ((16 * gS4) / (18 * s * s)) * (scalarProduct(pBeam2, ptbar) * (matrixGen(pt, pBeam1) + matrixGen(pt, pBeam1).T()));
    foo += -((16 * gS4) / (18 * s * s)) * (scalarProduct(pBeam1, pBeam2) * (matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T()));
    foo += -((16 * gS4) / (18 * s * s)) * ((scalarProduct(pt, ptbar) + mt2) * (matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T()));

    return foo;
}

//___________________________ P2g SME ___________________________//

TMatrixD MatrixElement::calculateP2g_munu(){
// Intermediate Element
    double s = normMinkowski2(pt + ptbar);
    double t = normMinkowski2(pBeam2 - pt);
    double u = normMinkowski2(pBeam1 - pt);

//*************************** deltap ****************************//

// deltapMss
    TMatrixD foo1(4,4);
    foo1  = ((3*gS4) / (4*s*s)) * (s * (matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T())
    + (t - mt2) * (matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2).T())
    + (u - mt2) * (matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2).T()));
// deltapMtt
    foo1 += (gS4 / (6*pow((t-mt2),3))) * ((-t * t + t * u - mt2 * u + 3 * mt2 * t - 10*mt2 * mt2) * (matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T())
    + (t * t + t * u - mt2 * u - 9*mt2 * t) * (matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T())
    + (-t * t - t * u + mt2 * u + 5*mt2 * t + 4*mt2 * mt2) * (matrixGen(pt, pBeam1) + matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2) + matrixGen(ptbar, pBeam2).T()));
// deltapMuu
    foo1 += (gS4 / (6*pow((u-mt2),3))) * ((-u * u + t * u - mt2 * t + 3*mt2 * u - 10*mt2 * mt2) * (matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T()) + (u * u + t * u - mt2 * t - 9*mt2 * u) * (matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T()) + (-u * u - t * u + mt2 * t + 5*mt2 * u + 4*mt2 * mt2) * (matrixGen(ptbar, pBeam1) + matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T()));
// deltapMst
    foo1 += ((3*gS4) / (32*s * pow((t-mt2),2))) * ((2*t*s - (t + mt2)*(u - mt2)) * (matrixGen(pt, pBeam1) + matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2) + matrixGen(ptbar, pBeam2).T())
    + (t-mt2) * ((3*t - 5*mt2) * (matrixGen(ptbar, pBeam1) + matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T()) + (t + 3*u - 8*mt2) * (matrixGen(pBeam1, pBeam2) - matrixGen(pBeam1, pBeam2).T() + matrixGen(ptbar, pBeam1) - matrixGen(ptbar, pBeam1).T() - matrixGen(ptbar, pBeam2) + matrixGen(ptbar, pBeam2).T()))
    - 2*(8*mt2*mt2 + (t-3*mt2)*(3*t+u)) * (matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T())
    + 4*(2*t*u - 3*mt2*t - mt2*u + 2*mt2*mt2) * (matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T()));
// deltapMsu
    foo1 += ((3*gS4) / (32*s * pow((u-mt2),2))) * ((2*u*s - (u + mt2)*(t - mt2)) * (matrixGen(ptbar, pBeam1) + matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T())
    + (u-mt2) * ((3*u - 5*mt2) * (matrixGen(pt, pBeam1) + matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2) + matrixGen(ptbar, pBeam2).T()) + (u + 3*t - 8*mt2) * (matrixGen(pBeam1, pBeam2) - matrixGen(pBeam1, pBeam2).T() + matrixGen(pt, pBeam1) - matrixGen(pt, pBeam1).T() - matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T()))
    - 2*(8*mt2*mt2 + (u-3*mt2)*(3*u+t)) * (matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T())
    + 4*(2*t*u - 3*mt2*u - mt2*t + 2*mt2*mt2) * (matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T()));
// deltapMtu
    foo1 += (gS4 / (24*pow((u-mt2),2)*pow((t-mt2),2))) * ((2*s+mt2)*(t-mt2)*(u-mt2)*(matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T() - matrixGen(pt, ptbar) - matrixGen(pt, ptbar).T())
    + mt2*((s*s - 7*mt2*s - 3*t*u + 3*mt2*mt2)*(matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T() + matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T())
    - (t-mt2)*(t-u+4*mt2)*(matrixGen(ptbar, pBeam1) + matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T())
    + (u-mt2)*(t-u-4*mt2)*(matrixGen(pt, pBeam1) + matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2) + matrixGen(ptbar, pBeam2).T())));
//**************************** deltav ***************************//

// deltavMss
    TMatrixD foo2(4,4);
    foo2  = ((3*gS4) / (4*s*s)) * (t*(matrixGen(ptbar, pBeam1) + matrixGen(pt, pBeam2) - matrixGen(pBeam1, pBeam2) - matrixGen(pBeam1, pBeam2).T())
    + u*(matrixGen(ptbar, pBeam2) + matrixGen(pt, pBeam1) - matrixGen(pBeam1, pBeam2) - matrixGen(pBeam1, pBeam2).T())
    - mt2*(matrixGenvMss(pBeam1, pBeam2)));
// deltavMtt
    foo2 += (gS4 / (3*pow(3*(t-mt2),2))) * ((t-3*mt2)*(matrixGen(pt, pBeam1) + matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2) + matrixGen(ptbar, pBeam2).T())
    + 4*mt2*(matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T()));
// deltavMuu
    foo2 += (gS4 / (3*pow(3*(u-mt2),2))) * ((u-3*mt2)*(matrixGen(ptbar, pBeam1) + matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T())
    + 4*mt2*(matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T()));
// deltavMst
    foo2 += ((3*gS4) / (32*s*(t-mt2))) * (2*(s+4*mt2)*(matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T())
    + (4*t+3*u-13*mt2)*(matrixGen(pt, pBeam1) + matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2) + matrixGen(ptbar, pBeam2).T())
    + 4*(t-u)*(matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T())
    + (-2*t-3*u+7*mt2)*(matrixGen(ptbar, pBeam1) + matrixGen(pt, pBeam2))
    + (3*u-9*mt2)*(matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2).T()));
// deltavMsu
    foo2 += ((3*gS4) / (32*s*(u-mt2))) * (2*(s+4*mt2)*(matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T())
    + (4*u+3*t-13*mt2)*(matrixGen(ptbar, pBeam1) + matrixGen(ptbar, pBeam1).T() + matrixGen(pt, pBeam2) + matrixGen(pt, pBeam2).T())
    + 4*(u-t)*(matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T())
    + (-2*u-3*t+7*mt2)*(matrixGen(pt, pBeam1) + matrixGen(ptbar, pBeam2))
    + (3*t-9*mt2)*(matrixGen(pt, pBeam1).T() + matrixGen(ptbar, pBeam2).T()));
// deltavMtu
    foo2 += (gS4 / (6*(t-mt2)*(u-mt2)))*(mt2*(matrixGen(pBeam1, pBeam2) + matrixGen(pBeam1, pBeam2).T()) - 2.*(matrixGen(pt, ptbar) + matrixGen(pt, ptbar).T()));

    return foo1 + foo2;
}

//_______________________________ F ______________________________//

TMatrixD MatrixElement::calculateF_munu(){
    TMatrixD foo(4,4);
    foo = ((2*gW4) / (mt2 * gammat2 * ((pow((normMinkowski2(pDecay1 + pAntiDecay1) - mW2),2)) + (mW2 * gammaW2)))) *
    (scalarProduct(pb, pt) * (matrixGen(pDecay1, pAntiDecay1) + matrixGen(pDecay1, pAntiDecay1).T())
    + scalarProduct(pb, pDecay1) * (matrixGen(pt, pAntiDecay1) + matrixGen(pt, pAntiDecay1).T())
    - scalarProduct(pb, pAntiDecay1) * (matrixGen(pt, pDecay1) + matrixGen(pt, pDecay1).T())
    - scalarProduct(pt, pDecay1) * (matrixGen(pb, pAntiDecay1) + matrixGen(pb, pAntiDecay1).T())
    + scalarProduct(pt, pAntiDecay1) * (matrixGen(pb, pDecay1) + matrixGen(pb, pDecay1).T())
    + scalarProduct(pDecay1, pAntiDecay1) * (matrixGen(pb, pt) + matrixGen(pb, pt).T()));
    return foo;
}

//_____________________________ Fbar ______________________________//

TMatrixD MatrixElement::calculateFbar_munu(){
    TMatrixD foo(4,4);
    foo = ((2*gW4) / (mt2 * gammat2 * ((pow((normMinkowski2(pDecay2 + pAntiDecay2) - mW2),2)) + (mW2 * gammaW2)))) *
    (scalarProduct(pbbar, ptbar) * (matrixGen(pDecay2, pAntiDecay2) + matrixGen(pDecay2, pAntiDecay2).T())
    + scalarProduct(ptbar, pDecay2) * (matrixGen(pbbar, pAntiDecay2) + matrixGen(pbbar, pAntiDecay2).T())
    - scalarProduct(ptbar, pAntiDecay2) * (matrixGen(pbbar, pDecay2) + matrixGen(pbbar, pDecay2).T())
    - scalarProduct(pbbar, pDecay2) * (matrixGen(ptbar, pAntiDecay2) + matrixGen(ptbar, pAntiDecay2).T())
    + scalarProduct(pbbar, pAntiDecay2) * (matrixGen(ptbar, pDecay2) + matrixGen(ptbar, pDecay2).T())
    + scalarProduct(pDecay2, pAntiDecay2) * (matrixGen(pbbar, ptbar) + matrixGen(pbbar, ptbar).T()));
    return foo;
}

/***************************************************************/

//____________________ Show outstream ___________________//

// Show SM
double MatrixElement::getPqqbar(){return Pqqbar;}
double MatrixElement::getP2g(){return P2g;}
double MatrixElement::getF(){return F;}
double MatrixElement::getFbar(){return Fbar;}

// Show SME
TMatrixD MatrixElement::getPqqbar_munu(){return Pqqbar_munu;}
TMatrixD MatrixElement::getP2g_munu(){return P2g_munu;}
TMatrixD MatrixElement::getF_munu(){return F_munu;}
TMatrixD MatrixElement::getFbar_munu(){return Fbar_munu;}

//______________________ Some tools _____________________//

double MatrixElement::normMinkowski2(TLorentzVector p){
    return p.Mag2();
}
double MatrixElement::scalarProduct(TLorentzVector p1, TLorentzVector p2){
    return p1*p2;
}

void MatrixElement::writeDownVecMatrix(double*** Vecmatrix, TString s, double nEvents){
    std::ofstream f(""+s, std::ios::out);
    for(int i=0; i<nEvents; i++){
        for(int j =0; j<4; j++){
            f<<Vecmatrix[i][j][0]<<" "<<Vecmatrix[i][j][1]<<" "<<Vecmatrix[i][j][2]<<" "<<Vecmatrix[i][j][3]<<std::endl;
        }
    }
    f.close();
}

void MatrixElement::writeNumMatrix(TString s, int n, TMatrixD m){
    std::ofstream f("data/matrix/"+s+".txt", std::ios::out);
    f<<s<<" "<<std::endl<<std::endl;
    f<<"Events n= "<<n<<std::endl<<std::endl;
    for(int j =0; j<4; j++)
        f<<m(j,0)<<" "<<m(j,1)<<" "<<m(j,2)<<" "<<m(j,3)<<std::endl;
    f.close();
}

void MatrixElement::writeArray(TString s, int n, double* array){
    std::ofstream f("data/pT/"+s+".txt", std::ios::out);
    for(int j =0; j<n; j++)
        f<<array[j]<<std::endl;
    f.close();
}

void MatrixElement::showMatrix(TMatrixD m){
   for(int i=0; i<4; i++)
      std::cout<<m(i,0)<<" "<<m(i,1)<<" "<<m(i,2)<<" "<<m(i,3)<<std::endl;
   std::cout<<std::endl;
}

//*****************  Matrix tools ******************//

TMatrixD MatrixElement::matrixAF(TMatrixD AF, TMatrixD AFbar, double F, double Fbar){
    TMatrixD foo(4,4);
    foo = (1/F)*AF;
    foo += (1/Fbar)*AFbar;
    return foo;
}

TMatrixD MatrixElement::matrixA(TMatrixD M, double m){
    TMatrixD foo(4,4);
    foo = (1/m)*M;
    return foo;
}

TMatrixD MatrixElement::matrixGen(TLorentzVector p1, TLorentzVector p2){
    TMatrixD foo(4,4);
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++)
            foo(i,j) = p1[i] * p2[j];
    }
    return foo;
}

TMatrixD MatrixElement::matrixGenvMss(TLorentzVector p1, TLorentzVector p2){
    TMatrixD foo(4,4);
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++)
            foo(i,j) = (p1-p2)[i] * (p1-p2)[j];
    }
    return foo;
}

TMatrixD MatrixElement::calculateMatrix(int ntot, int n, double*** vecMat){
    TMatrixD foo(4,4);
    for(int i=0; i<ntot; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                foo(j,k) += vecMat[i][j][k];
    foo = (1/((double)n)) * foo;
    return foo;
}
