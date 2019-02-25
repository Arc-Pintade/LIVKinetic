#include "../include/kineticAnalyze.hpp"

#include <iostream>

using namespace std;

int main (){

    cout<<"Hello World !"<<endl;

    KineticAnalyze k;
// 13 Tev is experiment 3
    for(int i=0; i<4; i++){
        k.fTime(i, 3);
        k.gTime(i,3);
        k.gTimeTT(i,3);
    }
    k.amplEnergy();

    return 0;
}
