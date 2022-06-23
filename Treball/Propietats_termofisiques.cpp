#include <vector>
#include <math.h>
#include "Propietats_termofisiques.h"


void fluid::Calcul_Coeficients_anular(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa)
{
    Reynolds=(rho*v*D)/mu;
    Prandt=mu*cp/lambda;
    double A,B;
    //Churchill expression
    A=pow(2.457*log(1/(pow(7/Reynolds,0.9)+0.27*rugositat_relativa)),16);
    B=pow(37530/Reynolds,16);
    fregament=2*pow(pow(8/Reynolds,12)+1/pow(A+B,1.5),1.0/12.0);

    Nusselt=0.023*pow(Reynolds,0.8)*pow(Prandt,0.4); //Cal mirar quin nusselt, expressio amb Dh o expressio 16 de la seccio C3
    Alfa_i=lambda*Nusselt/D;
}
    

void fluid::Propietats_termofisiquesO2(double T, double P,double Rgas){
    double c0=0.3782456360e+01, c1=-0.2996734160e-02, c3=0.9847302010e-05, c4=-0.9681295090e-08, c5=0.3243728370e-11;
    densitat=P/(Rgas*T);
    viscositat=(5.27117e-8*pow(T,3)-9.24802e-5*pow(T,2)+9.43953e-2*T-7.79576e-1)*pow(10,-6);
   
    conductivitat=2.33100e-15*pow(T,4)+7.77e-14*pow(T,3)-3.5148e-8*pow(T,2)+1.02599e-4*T-8.9303e-4;
}
void fluid::Propietats_termofisiquesH2(double T0, double Tf, double P,double Rgas){
    double c0=0.2344331120e+01, c1=0.7980520750e-02, c2=-0.1947815100e-04, c3=0.2015720940e-07, c4=-0.7376117610e-11;   //cp 200 a 1000K
    double c01000=0.3337279200e+01,c11000=-0.4940247310e-04, c21000=0.4994567780e-06, c31000= -0.1795663940e-09, c41000=0.2002553760e-13; //cp de 1000 a 5000K
    double m0=-0.1614293964e+02, m1=0.1003491326e+01, m2=-0.5016044555e-01, m3=0.2330995224e-02; //Lambda i mu valen el mateix de 200 a 5000K
    double l0=-0.2277096638e+01, l1=-0.4674267764e+00,l2=0.1156734789e+00,l3=-0.2596025563e-02;
    if (T0==Tf and T0<1000){
        cp=c0+c1*T0+c2*pow(T0,2)+c3*pow(T0,3)+c4*pow(T0,4);
    }
    else if (T0==Tf and T0>=1000)
    {
        cp=c01000+c11000*T0+c21000*pow(T0,2)+c31000*pow(T0,3)+c41000*pow(T0,4);
    }
    else if((T0<1000 and Tf<1000) and T0!=Tf) //Integral de T0 a Tf/(Tf-T0)
    {
        double terme1=c0*T0+c1*pow(T0,2)/2+c2*pow(T0,3)/3+c3*pow(T0,4)/4+c4*pow(T0,5)/5;
        double terme2=c0*Tf+c1*pow(Tf,2)/2+c2*pow(Tf,3)/3+c3*pow(Tf,4)/4+c4*pow(Tf,5)/5;
        cp=(terme2-terme1)/(Tf-T0);
    }
    else if ((T0<1000 and Tf>1000) and T0!=Tf) //(Integral de T0 a 1000 + de 1000 a Tf)/(Tf-T0)
    {
        double terme1=c0*T0+c1*pow(T0,2)/2+c2*pow(T0,3)/3+c3*pow(T0,4)/4+c4*pow(T0,5)/5;
        double terme2=c0*1000+c1*pow(1000,2)/2+c2*pow(1000,3)/3+c3*pow(1000,4)/4+c4*pow(1000,5)/5;
        double terme3=c01000*1000+c11000*pow(1000,2)/2+c21000*pow(1000,3)/3+c31000*pow(1000,4)/4+c41000*pow(1000,5)/5;
        double terme4=c01000*Tf+c11000*pow(Tf,2)/2+c21000*pow(Tf,3)/3+c31000*pow(Tf,4)/4+c41000*pow(Tf,5)/5;
        cp=(terme2-terme1+terme4-terme3)/(Tf-T0);
    }
    else if ((T0>1000 and Tf>1000) and T0!=Tf)
    {
        double terme1=c01000*T0+c11000*pow(T0,2)/2+c21000*pow(T0,3)/3+c31000*pow(T0,4)/4+c41000*pow(T0,5)/5;
        double terme2=c01000*Tf+c11000*pow(Tf,2)/2+c21000*pow(Tf,3)/3+c31000*pow(Tf,4)/4+c41000*pow(Tf,5)/5;
        cp=(terme2-terme1)/(Tf-T0);
    }
    

    double T=(T0+Tf)/2;
    densitat=P/(Rgas*T);
    viscositat= m0+m1*T+m2*pow(T,2)+m3*pow(T,3);
    conductivitat=l0+l1*T+l2*pow(T,2)+l3*pow(T,3);
    
}
void fluid::Propietats_termofisiquesaire(double T, double P,double Rgas){
    densitat=P/(Rgas*T);
    viscositat=1.458e-6*pow(T,1.5)/(T+110.4); 
    cp=1034.09-2.849e-1*T+7.817e-4*pow(T,2)-4.971e-7*pow(T,3)+1.077e-10*pow(T,4);
    conductivitat=(2.648e-3*sqrt(T))/(1+(245.4/T)*pow(10,-12.0/T));
    beta=1/T;
}
void fluid::Calcul_coeficients_exterior(double cp, double lambda, double mu, double g, double beta, double rho, double Tw, double Tf, double X){
    Prandt=mu*cp/lambda;
    Grashof=g*beta*pow(rho,2)*abs(Tw-Tf)*pow(X,3)/pow(mu,2);
    Rayleigh=Grashof*Prandt;
    if (Rayleigh<10e+9)
    {
        Nusselt=0.47*pow(Rayleigh,1.0/4.0);
    }
    else{
        Nusselt=0.1*pow(Rayleigh,1.0/3.0);
    }
    Alfa_i=lambda*Nusselt/X;
    
}

double condMolibde(double T){
    return(-5.58681e-8*pow(T,3)+1.57652e-4*pow(T,2)-1.56785e-1*T+1.58929e+2);
}