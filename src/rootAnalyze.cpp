#define rootAnalyze_cxx
#include "../include/rootAnalyze.hpp"
#include <iostream>


//_______________________ 4-momentum generating ____________________//

void Analyze::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Analyze.C
//      root> Analyze t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
//   Long64_t nentries = 100000;                /* just for test */

   q        = new TLorentzVector[nentries];
   qbar     = new TLorentzVector[nentries];
   g1       = new TLorentzVector[nentries];
   g2       = new TLorentzVector[nentries];
   t        = new TLorentzVector[nentries];
   tbar     = new TLorentzVector[nentries];
   b        = new TLorentzVector[nentries];
   lbar     = new TLorentzVector[nentries];
   nu       = new TLorentzVector[nentries];
   bbar     = new TLorentzVector[nentries];
   l        = new TLorentzVector[nentries];
   nubar    = new TLorentzVector[nentries];
   isFusion = new bool[nentries];

   etab     = new double[nentries];
   pTb      = new double[nentries];
   etalbar  = new double[nentries];
   pTlbar   = new double[nentries];
   etabbar  = new double[nentries];
   pTbbar   = new double[nentries];
   etal     = new double[nentries];
   pTl      = new double[nentries];

   pTt      = new double[nentries];
   pTtbar   = new double[nentries];

   for(int i=0; i<nentries; i++)
      isFusion[i] = false;
   bool is_g1;

//_______________________ 4-momentum filling ______________________//

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      is_g1 = false;
      for(int part = 0; part<kMaxParticle; part++){

//_______________________ Beam digluon fusion _____________________//

         if(Particle_PID[part]==21) {
            if(!is_g1){
               g1[jentry].SetPx(Particle_Px[part]);
               g1[jentry].SetPy(Particle_Py[part]);
               g1[jentry].SetPz(Particle_Pz[part]);
               g1[jentry].SetE(Particle_E[part]);
               isFusion[jentry] = true;
               is_g1 = true;
            }
            else{
               g2[jentry].SetPx(Particle_Px[part]);
               g2[jentry].SetPy(Particle_Py[part]);
               g2[jentry].SetPz(Particle_Pz[part]);
               g2[jentry].SetE(Particle_E[part]);
               isFusion[jentry] = true;
            }
         }

//________________________ Beam anihilation ______________________//

         if(Particle_PID[part]==1 || Particle_PID[part]==2 || Particle_PID[part]==3 || Particle_PID[part]==4) {
            q[jentry].SetPx(Particle_Px[part]);
            q[jentry].SetPy(Particle_Py[part]);
            q[jentry].SetPz(Particle_Pz[part]);
            q[jentry].SetE(Particle_E[part]);
         }
         if(Particle_PID[part]==-1 || Particle_PID[part]==-2 || Particle_PID[part]==-3 || Particle_PID[part]==-4) {
            qbar[jentry].SetPx(Particle_Px[part]);
            qbar[jentry].SetPy(Particle_Py[part]);
            qbar[jentry].SetPz(Particle_Pz[part]);
            qbar[jentry].SetE(Particle_E[part]);
         }

//________________________ top and anti-top ______________________//

         if(Particle_PID[part]==6) {
            t[jentry].SetPx(Particle_Px[part]);
            t[jentry].SetPy(Particle_Py[part]);
            t[jentry].SetPz(Particle_Pz[part]);
            t[jentry].SetE(Particle_E[part]);
            pTt[jentry]=Particle_PT[part];
         }
         if(Particle_PID[part]==-6) {
            tbar[jentry].SetPx(Particle_Px[part]);
            tbar[jentry].SetPy(Particle_Py[part]);
            tbar[jentry].SetPz(Particle_Pz[part]);
            tbar[jentry].SetE(Particle_E[part]);
            pTtbar[jentry]=Particle_PT[part];
         }

//____________________________ W⁺ decay _________________________//

         if(Particle_PID[part]==5) {
            b[jentry].SetPx(Particle_Px[part]);
            b[jentry].SetPy(Particle_Py[part]);
            b[jentry].SetPz(Particle_Pz[part]);
            b[jentry].SetE(Particle_E[part]);
            etab[jentry]=Particle_Eta[part];
            pTb[jentry]=Particle_PT[part];
         }
         if(Particle_PID[part]==-11 || Particle_PID[part]==-13) {
            lbar[jentry].SetPx(Particle_Px[part]);
            lbar[jentry].SetPy(Particle_Py[part]);
            lbar[jentry].SetPz(Particle_Pz[part]);
            lbar[jentry].SetE(Particle_E[part]);
            etalbar[jentry]=Particle_Eta[part];
            pTlbar[jentry]=Particle_PT[part];
         }
         if(Particle_PID[part]==12 || Particle_PID[part]==14) {
            nu[jentry].SetPx(Particle_Px[part]);
            nu[jentry].SetPy(Particle_Py[part]);
            nu[jentry].SetPz(Particle_Pz[part]);
            nu[jentry].SetE(Particle_E[part]);
         }

//____________________________ W⁻ decay _________________________//

         if(Particle_PID[part]==-5) {
            bbar[jentry].SetPx(Particle_Px[part]);
            bbar[jentry].SetPy(Particle_Py[part]);
            bbar[jentry].SetPz(Particle_Pz[part]);
            bbar[jentry].SetE(Particle_E[part]);
            etabbar[jentry]=Particle_Eta[part];
            pTbbar[jentry]=Particle_PT[part];
         }
         if(Particle_PID[part]==11 || Particle_PID[part]==13) {
            l[jentry].SetPx(Particle_Px[part]);
            l[jentry].SetPy(Particle_Py[part]);
            l[jentry].SetPz(Particle_Pz[part]);
            l[jentry].SetE(Particle_E[part]);
            etal[jentry]=Particle_Eta[part];
            pTl[jentry]=Particle_PT[part];
         }
         if(Particle_PID[part]==-12 || Particle_PID[part]==-14) {
            nubar[jentry].SetPx(Particle_Px[part]);
            nubar[jentry].SetPy(Particle_Py[part]);
            nubar[jentry].SetPz(Particle_Pz[part]);
            nubar[jentry].SetE(Particle_E[part]);
         }
      // if (Cut(ientry) < 0) continue;
      // check advancement
      }
      if(jentry%100000==0)
         std::cout<<jentry<<std::endl;
   }
}
