#include <math.h>
#include "Propietats_termofisiques.h"
#include "Combustio.h"
#include <iostream>
const double massa_molarH2=0.2015939951e-02, massa_molarO2=0.3199880028e-1;//Kg/mol
const double T0=298;
double Tcc(double cabalH2,double cabalO2,double TinH2, double TinO2,double delta,double Qlost,double Tcc,double &v_H2O,double &v_H2exces, double v_O2exces){
    fluid H2, O2,H2O;
    double entalpia_f_O2=0, entalpia_f_H2=0, entalpia_f_H2O=-241820; //J/mol
    double cabal_molarH2=cabalH2/massa_molarH2, cabal_molarO2=cabalO2/massa_molarO2; //mol/s
    double cabal_molarO2_exces=0, cabal_molarH20=0, cabal_molarH2_exces=0;
    double fr=0.3;//Factor de relaxacio    
    //Comprovem si es fuel rich o oxidizer rich
    if(cabal_molarO2/2>cabal_molarH2){ //Oxidizer rich
        cabal_molarH20=cabal_molarH2;
        cabal_molarO2_exces=cabal_molarO2-cabal_molarH20/2;
    }
    else if(cabal_molarO2/2<cabal_molarH2) //Fuel rich
    {
        cabal_molarH20=cabal_molarO2/2;
        cabal_molarH2_exces=cabal_molarH2-cabal_molarH20;
    }
    else //Estequiometric
    {
        cabal_molarH20=cabal_molarH2;
    }
    //Trobem coeficient estequiometric
    
    double v_H2=1, v_O2=cabal_molarO2/cabal_molarH2;
    v_H2O=cabal_molarH20/cabal_molarH2;
    v_O2exces=cabal_molarO2_exces/cabal_molarH2;
    v_H2exces=cabal_molarH2_exces/cabal_molarH2; 
   
    //Entalpies dels reactius
    H2.Propietats_termofisiquesH2(T0,TinH2,1,1);//P i R valen 1 perque la densitat i viscositat no ens importen
    O2.Propietats_termofisiquesO2(T0,TinO2,1,1);
    double hH2=H2.cp*(TinH2-T0), hO2=O2.cp*(TinO2-T0);
    double error=1;
    //Process iteratiu per trobar Tcc
    double num=v_H2*hH2+ v_O2*hO2- v_H2O*entalpia_f_H2O- massa_molarH2*Qlost/cabalH2;
    while (error>delta)
    {   
        double Tcc_old=Tcc;
        H2O.Propietats_termofisiquesH2O(T0,Tcc,1,1);
        O2.Propietats_termofisiquesO2(T0,Tcc,1,1);
        H2.Propietats_termofisiquesH2(T0,Tcc,1,1);
        Tcc=T0+(num)/(v_H2O*H2O.cp+v_H2exces*H2.cp+v_O2exces*O2.cp);
        Tcc=Tcc_old+fr*(Tcc-Tcc_old);
        error=abs(Tcc_old-Tcc);
    }
   
    return(Tcc);
}