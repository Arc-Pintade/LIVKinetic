#include "../include/rootAnalyze.hpp"
#include "../include/matrixConstruction.hpp"

#include <iostream>

using namespace std;

int main(){

    string cuts, simu;
    do{
        cout<<"Enter a simulation ('CMS13TeV', 'CMS7TeV', 'CMS2TeV' or 'D02TeV') : ";
        cin>>simu;
    }
    while(simu != "CMS13TeV" && simu != "CMS7TeV" && simu != "CMS2TeV" && simu != "D02TeV");
    cout<<endl;

    Analyze* an;

    if(simu == "CMS13TeV")
        an = new Analyze("EVENTS_LHE_ttbar_dilep_5M.root");
    else if(simu == "CMS7TeV")
        an = new Analyze("unweighted_events_ttbardilepLHC7TeV.root");
    else if(simu == "CMS2TeV")
        an = new Analyze("unweighted_events_ttbardilepLHC2TeV.root");
    else if(simu == "D02TeV")
        an = new Analyze("unweighted_events_ttbardilepTevatron.root");
    else{
        cout<<"error with simulation"<<endl;
        return 0;
    }
    an->Loop();
    int Events = an->fChain->GetEntriesFast();

    do{
        cout<<"Enter a cut parameters ('none', 'CMS', 'CMStop' or 'D0') : ";
        cin>>cuts;
    }
    while(cuts != "CMS" && cuts != "D0" && cuts != "none" && cuts != "CMStop");
    cout<<endl;

//__________________________ A matrices declarations _________________________//

    double*** APvecmatrixqqbar = new double** [n13TeV]{0};
    double*** APvecmatrix2g = new double** [n13TeV]{0};
    double*** AFvecmatrix = new double** [n13TeV]{0};
    for(int i=0; i<n13TeV; i++){
        APvecmatrixqqbar[i] = new double*[4]{0};
        APvecmatrix2g[i] = new double*[4]{0};
        AFvecmatrix[i] = new double*[4]{0};
            for (int j=0; j<4; j++){
                APvecmatrixqqbar[i][j] = new double[4]{0};
                APvecmatrix2g[i][j] = new double[4]{0};
                AFvecmatrix[i][j] = new double[4]{0};
            }
    }

//_________________________ Loop to creat .txt files _________________________//
//___________________________ with different cut _____________________________//

    int n2g=0;            int nqqbar=0;           int nF=0;
    int n2gtot=0;         int nqqbartot=0;        int nFtot=0;

    for (int i=0; i<Events; i++){
        if(i%10000==0)
            cout<<i<<endl;
        if(an->isFusion[i]){
            if(cuts == "none"){
                MatrixElement *M = new MatrixElement(an->g1[i], an->g2[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], true);
                MatrixElement *ME = new MatrixElement(an->g1[i], an->g2[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], true, i);
                for(int j=0; j<4; j++){
                    for(int k=0; k<4; k++){
                        APvecmatrix2g[i][j][k] = MatrixElement::matrixA(ME->getP2g_munu(), M->getP2g())[j][k];
                        AFvecmatrix[i][j][k] = MatrixElement::matrixAF(ME->getF_munu(), ME->getFbar_munu(), M->getF(),M->getFbar())[j][k];
                    }
                }
                n2g++;
                n2gtot++;
                nF++;
                nFtot++;
                delete M;
                delete ME;
            }
            else if (cuts == "CMS"){
                if(an->etal[i]<2.4 && an->etal[i]>-2.4 && an->etalbar[i]<2.4 && an->etalbar[i]>-2.4 && an->etab[i]<2.4 && an->etab[i]>-2.4 && an->etabbar[i]<2.4 && an->etabbar[i]>-2.4 && an->pTl[i]>20 && an->pTlbar[i]>20 && an->pTb[i]>30 && an->pTbbar[i]>30){
                    MatrixElement *M = new MatrixElement(an->g1[i], an->g2[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], true);
                    MatrixElement *ME = new MatrixElement(an->g1[i], an->g2[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], true, i);
                    for(int j=0; j<4; j++){
                        for(int k=0; k<4; k++){
                            APvecmatrix2g[i][j][k] = MatrixElement::matrixA(ME->getP2g_munu(), M->getP2g())[j][k];
                            AFvecmatrix[i][j][k] = MatrixElement::matrixAF(ME->getF_munu(), ME->getFbar_munu(), M->getF(),M->getFbar())[j][k];
                        }
                    }
                    n2g++;
                    n2gtot++;
                    nF++;
                    nFtot++;
                    delete M;
                    delete ME;
                }
                else{
                    n2gtot++;
                    nFtot++;
                }
            }
            else if ( cuts == "CMStop"){
                if(an->etal[i]<2.4 && an->etal[i]>-2.4 && an->etalbar[i]<2.4 && an->etalbar[i]>-2.4 && an->etab[i]<2.4 && an->etab[i]>-2.4 && an->etabbar[i]<2.4 && an->etabbar[i]>-2.4 && an->pTl[i]>20 && an->pTlbar[i]>20 && an->pTb[i]>30 && an->pTbbar[i]>30 && an->pTt[i]>50 && an->pTtbar[i]>50){
                    MatrixElement *M = new MatrixElement(an->g1[i], an->g2[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], true);
                    MatrixElement *ME = new MatrixElement(an->g1[i], an->g2[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], true, i);
                    for(int j=0; j<4; j++){
                        for(int k=0; k<4; k++){
                            APvecmatrix2g[i][j][k] = MatrixElement::matrixA(ME->getP2g_munu(), M->getP2g())[j][k];
                            AFvecmatrix[i][j][k] = MatrixElement::matrixAF(ME->getF_munu(), ME->getFbar_munu(), M->getF(),M->getFbar())[j][k];
                        }
                    }
                    n2g++;
                    n2gtot++;
                    nF++;
                    nFtot++;
                    delete M;
                    delete ME;
                }
                else{
                    n2gtot++;
                    nFtot++;
                }
            }
            else if (cuts == "D0"){

            }
            else
                cout<<"error with cuts"<<endl;
        }
        else{
           if(cuts == "none"){
                MatrixElement *M = new MatrixElement(an->q[i], an->qbar[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], false);
                MatrixElement *ME = new MatrixElement(an->q[i], an->qbar[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], false, i);
                for(int j=0; j<4; j++){
                    for(int k=0; k<4; k++){
                        APvecmatrixqqbar[i][j][k] = MatrixElement::matrixA(ME->getPqqbar_munu(), M->getPqqbar())[j][k];
                        AFvecmatrix[i][j][k] = MatrixElement::matrixAF(ME->getF_munu(), ME->getFbar_munu(), M->getF(),M->getFbar())[j][k];
                    }
                }
                nqqbar++;
                nqqbartot++;
                nF++;
                nFtot++;
                delete M;
                delete ME;
            }
            else if (cuts == "CMS"){
                if(an->etal[i]<2.4 && an->etal[i]>-2.4 && an->etalbar[i]<2.4 && an->etalbar[i]>-2.4 && an->etab[i]<2.4 && an->etab[i]>-2.4 && an->etabbar[i]<2.4 && an->etabbar[i]>-2.4 && an->pTl[i]>20 && an->pTlbar[i]>20 && an->pTb[i]>30 && an->pTbbar[i]>30){
                    MatrixElement *M = new MatrixElement(an->q[i], an->qbar[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], false);
                    MatrixElement *ME = new MatrixElement(an->q[i], an->qbar[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], false, i);
                    for(int j=0; j<4; j++){
                        for(int k=0; k<4; k++){
                            APvecmatrixqqbar[i][j][k] = MatrixElement::matrixA(ME->getPqqbar_munu(), M->getPqqbar())[j][k];
                            AFvecmatrix[i][j][k] = MatrixElement::matrixAF(ME->getF_munu(), ME->getFbar_munu(), M->getF(),M->getFbar())[j][k];
                        }
                    }
                    nqqbar++;
                    nqqbartot++;
                    nF++;
                    nFtot++;
                    delete M;
                    delete ME;
                }
                else{
                    nqqbartot++;
                    nFtot++;
                }
            }
            else if ( cuts == "CMStop"){
                if(an->etal[i]<2.4 && an->etal[i]>-2.4 && an->etalbar[i]<2.4 && an->etalbar[i]>-2.4 && an->etab[i]<2.4 && an->etab[i]>-2.4 && an->etabbar[i]<2.4 && an->etabbar[i]>-2.4 && an->pTl[i]>20 && an->pTlbar[i]>20 && an->pTb[i]>30 && an->pTbbar[i]>30 && an->pTt[i]>50 && an->pTtbar[i]>50){
                    MatrixElement *M = new MatrixElement(an->q[i], an->qbar[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], false);
                    MatrixElement *ME = new MatrixElement(an->q[i], an->qbar[i], an->t[i], an->tbar[i], an->b[i], an->nu[i], an->lbar[i], an->bbar[i], an->l[i], an->nubar[i], false, i);
                    for(int j=0; j<4; j++){
                        for(int k=0; k<4; k++){
                            APvecmatrixqqbar[i][j][k] = MatrixElement::matrixA(ME->getPqqbar_munu(), M->getPqqbar())[j][k];
                            AFvecmatrix[i][j][k] = MatrixElement::matrixAF(ME->getF_munu(), ME->getFbar_munu(), M->getF(),M->getFbar())[j][k];
                        }
                    }
                    nqqbar++;
                    nqqbartot++;
                    nF++;
                    nFtot++;
                    delete M;
                    delete ME;
                }
                else{
                    nqqbartot++;
                    nFtot++;
                }
            }
            else if (cuts == "D0"){

            }
            else
                cout<<"error with cuts"<<endl;
        }
    }

//________________________ Creation txt files _____________________//

    MatrixElement::writeDownVecMatrix(APvecmatrixqqbar, "VMAPqqbar"+cuts+".txt", n13TeV);
    cout<<"VMAPqqbar is wrote"<<endl;
    MatrixElement::writeDownVecMatrix(APvecmatrix2g, "VMAP2g"+cuts+".txt", n13TeV);
    cout<<"VMAP2g is wrote"<<endl;
    MatrixElement::writeDownVecMatrix(AFvecmatrix, "VMAF"+cuts+".txt",  n13TeV);
    cout<<"VMAF is wrote"<<endl;

//____________________________ Resume ____________________________//

    cout<<"Totaux :"<<endl;
    cout<<"cut number fusion : "<<n2g<<endl;
    cout<<"total number fusion : "<<n2gtot<<endl;
    cout<<"cut number anihilation : "<<nqqbar<<endl;
    cout<<"total number anihilation : "<<nqqbartot<<endl;
    cout<<"cut number decay : "<<nF<<endl;
    cout<<"total number decay : "<<nFtot<<endl;

    cout<<"End Bye Bro !"<<endl;
    return 0;
   //Bye.
}
