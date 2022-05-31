#include <vector>
#include <math.h>
#include "Propietats_termofisiques.h"
void fluid::Calcul_Coeficients(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa)
{
    Reynolds=(rho*v*D)/mu;
    Prandt=mu*cp/lambda;
    Nusselt=0.023*pow(Reynolds,0.8)*pow(Prandt,0.4);
    Alfa_i=lambda*Nusselt/D;
    double A,B;
    //Churchill expression
    A=pow(2.457*log(1/(pow(7/Reynolds,0.9)+0.27*rugositat_relativa)),16);
    B=pow(37530/Reynolds,16);
    fregament=2*pow(pow(8/Reynolds,12)+1/pow(A+B,1.5),1.0/12.0);
}

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
    
    
   

void fluid::Propietats_termofisiquesmescla(double T, double P,double Rgas,double cpA,double cpB,double fmolA,double fmolB, double muA, double muB, double conductivitatA, double conductivitatB){
    densitat=P/(Rgas*T);
    //Suposem que la viscositat, i conductivitat son la mitjana, de moment
    viscositat=(muA+muB)/2;
    cp=fmolA*cpA+cpB*fmolB;
    conductivitat=(conductivitatA+conductivitatB)/2;
}

void fluid::Propietats_termofisiquesO2(double T, double P,double Rgas){
    densitat=P/(Rgas*T);
    viscositat=(5.27117e-8*pow(T,3)-9.24802e-5*pow(T,2)+9.43953e-2*T-7.79576e-1)*pow(10,-6);
    cp=(5.501166e-12*pow(T,4)-9.681119e-9*pow(T,3)+6.499907e-6*pow(T,2)-1.693904e-3*T+1.0597)*pow(10,3);
    conductivitat=2.33100e-15*pow(T,4)+7.77e-14*pow(T,3)-3.5148e-8*pow(T,2)+1.02599e-4*T-8.9303e-4;
}
void fluid::Propietats_termofisiquesH2(double T, double P,double Rgas){
    densitat=P/(Rgas*T);
    viscositat= (5.83673e-9*pow(T,3)-1.7564449e-5*pow(T,2)+0.030645686887*T+1.113450152679)*pow(10,-6);
    cp=(-1.468867e-11*pow(T,4)+5.0109515e-8*pow(T,3)-5.914804e-5*pow(T,2)+2.918448e-2*T+9.472783)*pow(10,3);
    conductivitat=3.69347e-15*pow(T,4)+6.18450e-11*pow(T,3)-3.08777e-7*pow(T,2)+6.8004e-4*T+2.8372e-3;
    
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