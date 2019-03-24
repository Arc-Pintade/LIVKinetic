#include "../include/kineticAnalyze.hpp"

#include <iostream>

using namespace std;

int main (){

    cout<<"Hello World !"<<endl;
    double val = 0.01;

    KineticAnalyze k(val,val,val,val);
// 13 Tev is experiment 3
    for(int i=0; i<4; i++){
        k.fTime(i, 3);
        k.gTime(i,3);
        k.gTimeTT(i,3);
    }
    k.amplEnergy();

    return 0;
}
