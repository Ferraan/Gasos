#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "solvers.h"
#include "Propietats_termofisiques.h"
using namespace std;
auto start = chrono::high_resolution_clock::now();


const int n = 1000; //Volums de control del fluid, n+1 nodes
const double delta = 1e-10;
const double pi = 2 * acos(0.0);
const double Runiversal=8.3144621;
const double massa_molarH2=2e-3, massa_molarO2=32e-3;//Kg/mol
const double g=9.81;

int main(){
    std::cout.precision(15);
    std::scientific;
    double L=0.5, Di=0.02, Do=0.023, Dm=0.033, DOut=0.035;
    double Tin1=1000, pin1=161e+5, Tin3=100, cabalin3=0.005, pin3=150e+5, Text=300, Pext=101325;
    double cabalin1H2=0.005, cabalin1O2=0.05, cabalin1Tot=cabalin1H2+cabalin1O2;
    double Ttub2_inic=800, Ttub4_inic=500;
    double rugositat2in=0.0002, rugositat2ext=0.0003, rugositat4in=0.0004, rugositat4ext=0.0002;
    
    //Calculs previs
    double Deltax=L/n, Dh=Dm-Do;
    double S1=pi*Di*Di/4, S3=pi*pow(Dm-Do,2)/4, Sl2int=pi*Di*Deltax, Sl2out=pi*Do*Deltax, Sl4in=pi*Dm*Deltax, Sl4out=pi*DOut*Deltax, Sx2=pi*(pow(Do/2,2)-pow(Di/2,2)),  Sx4=pi*(pow(DOut/2,2)-pow(Dm/2,2));
    double molsH2=cabalin1H2/massa_molarH2, molsO2=cabalin1O2/massa_molarO2;
    double fraccio_molarH2=molsH2/(molsH2+molsO2), fraccio_molarO2=molsO2/(molsH2+molsO2), massa_molar_cambra=fraccio_molarH2*massa_molarH2+fraccio_molarO2*massa_molarO2;
    double Rgas_cambra=Runiversal/massa_molar_cambra, Rhidrogen=Runiversal/massa_molarH2, Roxigen=Runiversal/massa_molarO2, Raire=287;
    double rhoin1=pin1/(Tin1*Rgas_cambra), rhoin3=pin3/(Tin3*Rhidrogen);
    double vin1=cabalin1Tot/(S1*rhoin1), vin3=cabalin3/(S3*rhoin3);
    std::vector<double> x1(n+1,0), v1(n+1,0), T1(n+1,0), p1(n+1,0), rho1(n+1,0), v3(n+1,0),T3(n+1,0), p3(n+1,0), rho3(n+1,0), T2(n,Ttub2_inic), T4(n,Ttub4_inic), alfa1(n,0), alfa3(n,0), alfa5(n,0);
    x1[0]=0;
    for (int i = 1; i < n+1; i++)
    {
        x1[i]=x1[i-1]+Deltax;
    }
    
   
    //Zona 1 i 3
    v1[0]=vin1; T1[0]=Tin1; p1[0]=pin1; rho1[0]=rhoin1;
    v3[0]=vin3; T3[0]=Tin3; p3[0]=pin3; rho3[0]=rhoin3;
    fluid H2cambra, O2cambra, H2ext, mescla_cambra, aire;
    double errorext=1;
    while(errorext>delta){
        for (int i=1; i<n+1; i++){
            v1[i]=v1[i-1]; p1[i]=p1[i-1]; T1[i]=T1[i-1]; rho1[i]=rho1[i-1];
            
            double error=1.0;
        
            while (error>delta){     
                double pold=p1[i], Told=T1[i], vold=v1[i], rhoold=rho1[i]; 
                double Ti=(T1[i]+T1[i-1])/2, Pi=(p1[i-1]+p1[i])/2,vi=(v1[i]+v1[i-1])/2,rhoi=(rho1[i]+rho1[i-1])/2;
                H2cambra.Propietats_termofisiquesH2(Ti,Pi,Rhidrogen); 
                O2cambra.Propietats_termofisiquesO2(Ti,Pi,Roxigen);
                mescla_cambra.Propietats_termofisiquesmescla(Ti,Pi,Rgas_cambra,H2cambra.cp,O2cambra.cp,fraccio_molarH2,fraccio_molarO2,H2cambra.viscositat,O2cambra.viscositat,H2cambra.conductivitat,H2cambra.conductivitat);
                mescla_cambra.Calcul_Coeficients(mescla_cambra.viscositat,mescla_cambra.cp,mescla_cambra.conductivitat,Di,vi,mescla_cambra.densitat,rugositat2in);
                alfa1[i-1]=mescla_cambra.Alfa_i;
                p1[i]=-cabalin1Tot*(v1[i]-v1[i-1])/S1+p1[i-1]-mescla_cambra.fregament*rhoi*pow(vi,2)/(2*S1)*pi*Di*Deltax;
                double cA=cabalin1Tot*(pow(v1[i],2)-pow(v1[i-1],2))/2.0, cB=cabalin1Tot*mescla_cambra.cp+mescla_cambra.Alfa_i*Sl2int/2; //Cal canviar el cp a la mitjana
                //cout<<cA<<"  "<<cB<<endl;
                T1[i]=(T1[i-1]*cabalin1Tot*mescla_cambra.cp-cA+Sl2int*mescla_cambra.Alfa_i*(T2[i]-T1[i-1]/2))/cB;
                rho1[i]=p1[i]/(T1[i]*Rgas_cambra);
                v1[i]=cabalin1Tot/(rho1[i]*S1);
                error=max(abs(p1[i]-pold), abs(Told-T1[i]));
                error=max(abs(vold-v1[i]),error);
                error=max(abs(rho1[i]-rhoold),error);
            }
            //Zona 3
            v3[i]=v3[i-1]; p3[i]=p3[i-1]; T3[i]=T3[i-1]; rho3[i]=rho3[i-1];
            error=1;
            while (error>delta){     
                double pold=p3[i], Told=T3[i], vold=v3[i], rhoold=rho3[i]; 
                double Ti=(T3[i]+T3[i-1])/2, Pi=(p3[i-1]+p3[i])/2,vi=(v3[i]+v3[i-1])/2,rhoi=(rho3[i]+rho3[i-1])/2;
                H2ext.Propietats_termofisiquesH2(Ti,Pi,Rhidrogen);    
                H2ext.Calcul_Coeficients(H2ext.viscositat,H2ext.cp,H2ext.conductivitat,Dh,vi,H2ext.densitat,rugositat2ext);
                alfa3[i-1]=H2ext.Alfa_i;
                p3[i]=-cabalin3*(v3[i]-v3[i-1])/S3+p3[i-1]-H2ext.fregament*rhoi*pow(vi,2)/(2*S3)*pi*Dh*Deltax;
                double cA=cabalin1Tot*(pow(v1[i],2)-pow(v1[i-1],2))/2.0, cB=cabalin3*H2ext.cp+H2ext.Alfa_i/2 *(Sl2out+Sl4in);
                T3[i]=(T3[i-1]*(cabalin3*H2ext.cp-H2ext.Alfa_i/2*(Sl2out+Sl4in))-cA+H2ext.Alfa_i*(T2[i-1]*Sl2out+T4[i-1]*Sl4in))/cB;
                rho3[i]=p3[i]/(T3[i]*Rhidrogen);
                v3[i]=cabalin3/(rho3[i]*S3);
                error=max(abs(p3[i]-pold), abs(Told-T3[i]));
                error=max(abs(vold-v3[i]),error);
                error=max(abs(rho3[i]-rhoold),error);
            }
            //Zona exterior
            double Tm=(T4[i]+Text)/2;
            aire.Propietats_termofisiquesaire(Tm,Pext,Raire);
            aire.Calcul_coeficients_exterior(aire.cp,aire.conductivitat,aire.viscositat,g,aire.beta,aire.densitat,T4[i],Text,DOut);
            alfa5[i-1]=aire.Alfa_i;
        }
        std::vector<double> aP(n), aE(n), aW(n), bP(n);
        //Tub 2 
        double T2old=T2[n/2];    
        aE[0]=-Deltax/(Deltax/(2*condMolibde(T2[0]))+Deltax/(2*condMolibde(T2[1])))*Sx2/Deltax;
        aW[0]=0; //Adiabatic
        bP[0]=-(T1[0]+T1[1])/2*alfa1[0]*Sl2int-alfa3[0]*(T3[0]+T3[1])/2*Sl2out;
        aP[0]=aW[0]+aE[0]-alfa1[0]*Sl2int-alfa3[0]*Sl2out;
        for (int i = 1; i < n-1; i++)
        {
            aE[i]=-Deltax/(Deltax/(2*condMolibde(T2[i]))+Deltax/(2*condMolibde(T2[i+1])))*Sx2/Deltax;
            aW[i]=Deltax/(Deltax/(2*condMolibde(T2[i]))+Deltax/(2*condMolibde(T2[i-1])))*Sx2/Deltax; //Adiabatic
            bP[i]=-(T1[i]+T1[i+1])/2*alfa1[i]*Sl2int-alfa3[i]*(T3[i]+T3[i+1])/2*Sl2out;
            aP[i]=aW[i]+aE[i]-alfa1[i]*Sl2int-alfa3[i]*Sl2out;
        }
        aE[n-1]=0;//Adiabatic
        aW[n-1]=Deltax/(Deltax/(2*condMolibde(T2[n-2]))+Deltax/(2*condMolibde(T2[n-1])))*Sx2/Deltax;
        bP[n-1]=-(T1[n-1]+T1[n])/2*alfa1[n-1]*Sl2int-alfa3[n-1]*(T3[n-1]+T3[n])/2*Sl2out;
        aP[n-1]=aW[n-1]+aE[n-1]-alfa1[n-1]*Sl2int-alfa3[n-1]*Sl2out;
        solverTDMA(T2,aP,aW,aE,bP,n);
        //Tub 4
        double T4old=T4[n/2];
        aE[0]=-Deltax/(Deltax/(2*condMolibde(T4[0]))+Deltax/(2*condMolibde(T4[1])))*Sx4/Deltax;
        aW[0]=0; //Adiabatic
        bP[0]=-(T3[0]+T3[1])/2*alfa3[0]*Sl4in-alfa5[0]*Text*Sl4out;
        aP[0]=aW[0]+aE[0]-alfa3[0]*Sl4in-alfa5[0]*Sl4out;
        for (int i = 1; i < n-1; i++)
        {
            aE[i]=-Deltax/(Deltax/(2*condMolibde(T4[i]))+Deltax/(2*condMolibde(T4[i+1])))*Sx4/Deltax;
            aW[i]=Deltax/(Deltax/(2*condMolibde(T4[i]))+Deltax/(2*condMolibde(T4[i-1])))*Sx4/Deltax; //Adiabatic
            bP[i]=-(T3[i]+T3[i+1])/2*alfa3[i]*Sl4in-alfa5[i]*Text*Sl4out;
            aP[i]=aW[i]+aE[i]-alfa3[i]*Sl4in-alfa5[i]*Sl4out;
        }
        aE[n-1]=0;//Adiabatic
        aW[n-1]=Deltax/(Deltax/(2*condMolibde(T4[n-2]))+Deltax/(2*condMolibde(T4[n-1])))*Sx4/Deltax;
        bP[n-1]=-(T3[n-1]+T3[n])/2*alfa3[n-1]*Sl4in-alfa5[n-1]*Text*Sl4out;
        aP[n-1]=aW[n-1]+aE[n-1]-alfa3[n-1]*Sl4in-alfa5[n-1]*Sl4out;
        solverTDMA(T4,aP,aW,aE,bP,n);
        errorext=max(abs(T4[n/2]-T4old),abs(T2[n/2]-T2old));
        
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout<<"Temps execucio (s)" <<static_cast<float>(duration.count())/1000000 << endl;
    
    
    ofstream fout;
    fout.open("Treball_sortida.csv");
    fout<<"i"<<","<<"x1[i]"<<","<<"T1[i]"<<","<<"P1[i]"<<","<<"v1[i]"<<","<<"rho1[i]"<<","<<"alfa1[i]"<<","<<"T3[i]"<<","<<"P3[i]"<<","<<"v3[i]"<<","<<"rho3[i]"<<","<<"alfa3[i]"<<","<<"T2[i]"<<","<<"T4[i]"<<endl; //alfa, T2,T4[n-1]=0 perque no hi ha nodes
    for (int i = 0; i < n+1; i++)
    {   
        fout<<setprecision(15)<<i<<","<<x1[i]<<","<<T1[i]<<","<<p1[i]<<","<<v1[i]<<","<<rho1[i]<<","<<alfa1[i]<<","<<T3[i]<<","<<p3[i]<<","<<v3[i]<<","<<rho3[i]<<","<<alfa3[i]<<","<<T2[i]<<","<<T4[i]<<endl;
        
        
    }
    
}



