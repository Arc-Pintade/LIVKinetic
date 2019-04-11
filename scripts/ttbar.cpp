#include "../include/kineticAnalyze.hpp"

#include <iostream>

using namespace std;

int main (){
    cout<<"Hello World !"<<endl;
    double wilson = 0.01;
    string XX[4];
    XX[0] = "XX";    XX[1] = "XY";    XX[2] = "XZ";    XX[3] = "YZ";

    KineticAnalyze k(wilson,wilson,wilson,wilson);
// 13 Tev is experiment 3
    for(int i=0; i<4; i++){
        k.fTime(i, 3);
        k.gTime(i,3);
        k.gTimeTT(i,3);
        k.compareFusAni(i,3);
        k.compareCMSD0(i);
        k.earthSignal(XX[i]);
    }

    k.amplEnergyComparaison(true);
//    k.amplEnergyComparaison(false);

    return 0;
}