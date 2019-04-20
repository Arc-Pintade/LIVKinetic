/////////////////////////////////////////////////////
//                                                 //
// Program Created by Arc-Pintade (Aurélien CARLE) //
// Thanks to Plut0n (Xavier Valcarce) for his help //
//                                                 //
/////////////////////////////////////////////////////

#ifndef const_h
#define const_h
#define const_cxx

#include <cmath>

    const int n13TeV = 5000000;
    const int n100TeV = 1000000;

//____________ MMSM ALL IN GeV ______________//

// Strong Interaction Coupling
    const double gS2 = 4*M_PI*0.1182;
    const double gS4 = gS2 * gS2;
// Weak Interaction Coupling
    const double gW2 = 4*M_PI*10e-5;
    const double gW4 = gW2 * gW2;
// Top Mass
    const double mt = 173.21;         //Uncertainties ~0.71
    const double mt2 = mt * mt;
// Top Width
    const double gammat = 1.41;       //Uncertainties ~ +0.19-0.15
    const double gammat2 = gammat * gammat;
// W Mass
    const double mW = 80.385;         //Uncertainties ~0.015  
    const double mW2 = mW * mW;
// W Width
    const double gammaW = 2.085;      //Uncertainties ~0.042
    const double gammaW2 = gammaW * gammaW;

//______________ Rotations ______________//

// Latitude
   const double latitudeCMS = (46.3099/180) * M_PI;    //Uncertainties ~0.0001
// Longitude
   const double longitudeCMS = (6.0766/180) * M_PI;    //Uncertainties ~0.0001
// Azimuth
   const double azimuthCMS = (101.2794/180) * M_PI;     //Uncertainties ~0.003
// Tilt
   const double tiltCMS = (0.20632/180) * M_PI;     //Uncertainties ~0.00003
// Universal Frequence
   const double omega = 7.29211515e-5;           // en rad/s

// Latitude
   const double latitudeATLAS = (46.2357/180) * M_PI;    //Uncertainties ~0.0001
// Longitude
   const double longitudeATLAS = (6.0553/180) * M_PI;    //Uncertainties ~0.0001
// Azimuth
   const double azimuthATLAS = (281.265/180) * M_PI;     //Uncertainties ~0.003

//______________ Rotation D0______________//

// Latitude D0
   const double colatitD0 = (49.8255/180) * M_PI;    //Uncertainties ~0.0001
   const double latitD0 = M_PI/2 - colatitD0;
// Azimuth D0 (Orientation of proton beam)
   const double azimD0 = (42.192/180) * M_PI;     //Uncertainties ~0.003
   const double azimD0Redifine = (222.192/180) * M_PI;

// Rotation Matrix Values
   const double s1D0 = sin(colatitD0);
   const double c1D0 = cos(colatitD0);
   const double s2D0 = sin(azimD0);
   const double c2D0 = cos(azimD0);

#endif
