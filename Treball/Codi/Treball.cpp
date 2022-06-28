#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <bits/stdc++.h>
#include "solvers.h"
#include "Propietats_termofisiques.h"
#include "Combustio.h"
using namespace std;
auto start = chrono::high_resolution_clock::now();


const int n =50; //Volums de control del fluid, n+1 nodes
const double delta = 1e-10;
const double pi = 2 * acos(0.0);
const double Runiversal=8.3144621;
const double massa_molarH2=0.2015939951e-2, massa_molarO2=0.3199880028e-1, massa_molarH2O=0.1801534009e-1;//Kg/mol
const double g=9.81;

int main(){
    std::cout.precision(15);
    std::scientific;
    double L=0.5, D1=0.2, D2=0.203, D3=0.25, D4=0.253;
    double pcc=100e+5, cabalin1H2=41.2, cabalin1O2=193.8, TinH2=200, TinO2=200; //La temperatura d'entrada no és la temperatura real, però si és inferior a 200K, el rang de cp ja no és vàlid
    double rhocambra=70.85;//Kg/m3 //Primer assumirem la del hidrogen líquid, després la modificarem amb la dels gasos després de la combustió i necessitarem calcular Tcc
    
    double Tin3=200, cabalin3=0.05, pin3=150e+5, Text=300, Pext=101325;
    double Tcomb=1000, Ttub2_inic=800, Ttub4_inic=500;
    double rugositat2ext=0.0003, rugositat4in=0.0004, rugositat4ext=0.0002;
    double Q_cambra=0, Tcc_inici=2500;
    //Calculs previs
    double Deltax=L/n, Dh=D3-D2;
    double S1=pi*D1*D1/4, S3=pi*pow(D3-D2,2)/4, Sl2int=pi*D1*Deltax, Sl2out=pi*D2*Deltax, Sl4in=pi*D3*Deltax, Sl4out=pi*D4*Deltax, Sx2=pi*(pow(D2/2,2)-pow(D1/2,2)),  Sx4=pi*(pow(D4/2,2)-pow(D3/2,2));
    double Rhidrogen=Runiversal/massa_molarH2, Roxigen=Runiversal/massa_molarO2, Raire=287;
    double rhoin3=pin3/(Tin3*Rhidrogen);
    double vin3=cabalin3/(S3*rhoin3);
    double v_cambra=cabalin1H2/(S1*rhocambra); //No sabem la velocitat dels gasos després de la combustió. Assumirem que és la mateixa a la que entra el fluid.
    std::vector<double> x1(n+1,0), v3(n+1,0),T3(n+1,0), p3(n+1,0), rho3(n+1,0), T2(n,Ttub2_inic), T4(n,Ttub4_inic), alfa3(n,0), alfa5(n,0), f3(n,0);
    double alfa1=100;
    x1[0]=0;
    for (int i = 1; i < n+1; i++)
    {
        x1[i]=x1[i-1]+Deltax;
    }
    double v_H2O,v_H2exces,v_O2exces;
    Tcomb=Tcc(cabalin1H2,cabalin1O2,TinH2,TinO2,delta,Q_cambra,Tcc_inici,v_H2O,v_H2exces,v_O2exces); //Calculem les fraccions molars
    double massa_molar_cambra=v_H2O*massa_molarH2O+v_H2exces*massa_molarH2+v_O2exces*massa_molarO2;
    double Rcambra=Runiversal/massa_molar_cambra;
    //Zona 1 i 3
    v3[0]=vin3; T3[0]=Tin3; p3[0]=pin3; rho3[0]=rhoin3;
    fluid H2ext, mescla_cambra, aire;  
    double errorext=1;
    while(errorext>delta){
        Tcomb=Tcc(cabalin1H2,cabalin1O2,TinH2,TinO2,delta,Q_cambra,Tcc_inici,v_H2O,v_H2exces,v_O2exces);
        mescla_cambra.Propietats_termofisiquescambra(Tcomb,pcc,Rcambra,v_H2O,v_H2exces,v_O2exces);
        mescla_cambra.Calcul_var_adim(mescla_cambra.viscositat,mescla_cambra.cp,mescla_cambra.conductivitat,D1,v_cambra,mescla_cambra.densitat,0);
        alfa1=mescla_cambra.Alfa_i;
        Q_cambra=0;
        for (int i = 0; i < n; i++)
        {    
            Q_cambra=alfa1*Sl2int*(Tcomb-T2[i])+Q_cambra;
        }
        
        for (int i=1; i<n+1; i++){           
            //Zona 3
            v3[i]=v3[i-1]; p3[i]=p3[i-1]; T3[i]=T3[i-1]; rho3[i]=rho3[i-1];
            double error=1.0;
            while (error>delta){     
                double pold=p3[i], Told=T3[i], vold=v3[i], rhoold=rho3[i]; 
                double Ti=(T3[i]+T3[i-1])/2, Pi=(p3[i-1]+p3[i])/2,vi=(v3[i]+v3[i-1])/2,rhoi=(rho3[i]+rho3[i-1])/2;
                H2ext.Propietats_termofisiquesH2(T3[i-1],T3[i],Pi,Rhidrogen);    
                H2ext.cp=H2ext.cp/massa_molarH2;
                H2ext.Calcul_var_adim(H2ext.viscositat,H2ext.cp,H2ext.conductivitat,Dh,vi,H2ext.densitat,rugositat2ext);
                alfa3[i-1]=H2ext.Alfa_i;
                p3[i]=-cabalin3*(v3[i]-v3[i-1])/S3+p3[i-1]-0.5*H2ext.fregament*rhoi*pow(vi,2)*(Sl2out+Sl4in)/S3;
                f3[i-1]=H2ext.fregament;
                double cA=cabalin3*(pow(v3[i],2)-pow(v3[i-1],2))/2.0, cB=cabalin3*H2ext.cp+H2ext.Alfa_i/2 *(Sl2out+Sl4in);
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
            aire.Calcul_coeficients_exterior(aire.cp,aire.conductivitat,aire.viscositat,g,aire.beta,aire.densitat,T4[i],Text,D4);
            alfa5[i-1]=aire.Alfa_i;
            
            
        }
        std::vector<double> aP(n), aE(n), aW(n), bP(n);
        //Tub 2 
        std::vector<double> T2old(n),T4old(n), D1f(n);
        T2old=T2;
        T4old=T4;
        aE[0]=Deltax/(Deltax/(2*condMolibde(T2[0]))+Deltax/(2*condMolibde(T2[1])))*Sx2/Deltax;
        aW[0]=0; //AD1abatic
        bP[0]=Tcomb*alfa1*Sl2int+alfa3[0]*(T3[0]+T3[1])/2*Sl2out;
        aP[0]=aW[0]+aE[0]+alfa1*Sl2int+alfa3[0]*Sl2out;
        for (int i = 1; i < n-1; i++)
        {
            aE[i]=Deltax/(Deltax/(2*condMolibde(T2[i]))+Deltax/(2*condMolibde(T2[i+1])))*Sx2/Deltax;
            aW[i]=Deltax/(Deltax/(2*condMolibde(T2[i]))+Deltax/(2*condMolibde(T2[i-1])))*Sx2/Deltax; 
            bP[i]=Tcomb*alfa1*Sl2int+alfa3[i]*(T3[i]+T3[i+1])/2*Sl2out;
            aP[i]=aW[i]+aE[i]+alfa1*Sl2int+alfa3[i]*Sl2out;
        }
        aE[n-1]=0;//AD1abatic
        aW[n-1]=Deltax/(Deltax/(2*condMolibde(T2[n-2]))+Deltax/(2*condMolibde(T2[n-1])))*Sx2/Deltax;
        bP[n-1]=Tcomb*alfa1*Sl2int+alfa3[n-1]*(T3[n-1]+T3[n])/2*Sl2out;
        aP[n-1]=aW[n-1]+aE[n-1]+alfa1*Sl2int+alfa3[n-1]*Sl2out;
        solverTDMA(T2,aP,aW,aE,bP,n);
        //Tub 4
        aE[0]=Deltax/(Deltax/(2*condMolibde(T4[0]))+Deltax/(2*condMolibde(T4[1])))*Sx4/Deltax;
        aW[0]=0; //AD1abatic
        bP[0]=(T3[0]+T3[1])/2*alfa3[0]*Sl4in+alfa5[0]*Text*Sl4out;
        aP[0]=aW[0]+aE[0]+alfa3[0]*Sl4in+alfa5[0]*Sl4out;
        for (int i = 1; i < n-1; i++)
        {
            aE[i]=Deltax/(Deltax/(2*condMolibde(T4[i]))+Deltax/(2*condMolibde(T4[i+1])))*Sx4/Deltax;
            aW[i]=Deltax/(Deltax/(2*condMolibde(T4[i]))+Deltax/(2*condMolibde(T4[i-1])))*Sx4/Deltax; 
            bP[i]=(T3[i]+T3[i+1])/2*alfa3[i]*Sl4in+alfa5[i]*Text*Sl4out;
            aP[i]=aW[i]+aE[i]+alfa3[i]*Sl4in+alfa5[i]*Sl4out;
        }
        aE[n-1]=0;//AD1abatic
        aW[n-1]=Deltax/(Deltax/(2*condMolibde(T4[n-2]))+Deltax/(2*condMolibde(T4[n-1])))*Sx4/Deltax;
        bP[n-1]=(T3[n-1]+T3[n])/2*alfa3[n-1]*Sl4in+alfa5[n-1]*Text*Sl4out;
        aP[n-1]=aW[n-1]+aE[n-1]+alfa3[n-1]*Sl4in+alfa5[n-1]*Sl4out;
        solverTDMA(T4,aP,aW,aE,bP,n);
        //Error  
        for (int i = 0; i < n-1; i++)
        {
            D1f[i]=max(abs(T2old[i]-T2[i]),abs(T4old[i]-T4[i]));
        }
        errorext=*max_element(D1f.begin(),D1f.end() );
    }
    //Validacions
    //Calor en els tubs=0
    double Q_tub2=0,Q_tub4=0;
    for (int i = 0; i < n; i++)
    {    
        Q_tub2=alfa1*Sl2int*(Tcomb-T2[i])-alfa3[i]*Sl2out*(-(T3[i]+T3[i+1])/2+T2[i])+Q_tub2;
        Q_tub4=alfa3[i]*Sl4in*((T3[i]+T3[i+1])/2-T4[i])-alfa5[i]*Sl4out*(-Text+T4[i])+Q_tub4;
    }
    double Q_34=0;
    for (int i = 0; i < n; i++)
    {    
        Q_34=alfa3[i]*Sl4in*((T3[i]+T3[i+1])/2-T4[i])+Q_34;
    }
    std::cout<<"Tub 2: "<<Q_tub2<<", Tub 4: "<<Q_tub4<<endl;
    //Conservacio massa tub
    double cons_massa3=rho3[0]*v3[0]-rho3[n]*v3[n];
    std::cout<<"Massa:"<<cons_massa3<<endl;
    //Momentum tub 3
    /* double mom_T=0;
    for (int i = 0; i < n-1; i++)
    {
        mom_T=p3[i]*S3-p3[i+1]*S3-0.5*f3[i]*(Sl2out+Sl4in)*(rho3[i]+rho3[i+1])/2*pow((v3[i]+v3[i+1])/2,2)-cabalin3*(v3[i+1]-v3[i])+mom_T;
    }
    cout<<"Moment: "<<mom_T<<endl; */
    //Momentum tub 3 balanç global
    double fregament_tot=0;
    for (int i = 0; i < n; i++)
    {
        fregament_tot=0.5*f3[i]*(Sl2out+Sl4in)*(rho3[i]+rho3[i+1])/2*pow((v3[i]+v3[i+1])/2,2)+fregament_tot;
    }
    double mom_2=p3[0]*S3-p3[n]*S3-fregament_tot-cabalin3*(v3[n]-v3[0]);
    cout<<"Moment: "<<mom_2<<endl;
    //Energia tub 3 volum de control a volum de control
    /* double energiatot=0;
    for (int i = 0; i < n; i++)
    {   
        H2ext.Propietats_termofisiquesH2(T3[i],T3[i+1],p3[i],Rhidrogen);     
        H2ext.cp=H2ext.cp/massa_molarH2;
        energiatot=alfa3[i]*Sl2out*(-(T3[i]+T3[i+1])/2+T2[i])-alfa3[i]*Sl4in*((T3[i]+T3[i+1])/2-T4[i])-cabalin3*H2ext.cp*(T3[i+1]-T3[i])-cabalin3*(pow(v3[i+1],2)-pow(v3[i],2))/2+energiatot;
    }
    cout<<"Energia intercanviada en el tub:"<<energiatot<<endl; */
    //Energia balanç global
    H2ext.Propietats_termofisiquesH2(T3[0],T3[n],p3[0]+p3[n]/2,Rhidrogen);
    H2ext.cp=H2ext.cp/massa_molarH2;
    double Energia=Q_cambra-Q_34-cabalin3*H2ext.cp*(T3[n]-T3[0])-cabalin3*((pow(v3[n],2)-pow(v3[0],2))/2);
    cout<<"Consercacio de l'energia: "<<Energia<<endl;
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    std::cout<<"Temps execucio (s)" <<static_cast<float>(duration.count())/1000000 << endl; 
    

    ofstream fout;
    fout.open("../Sortida/Treball_sortida.csv");
    fout<<"i"<<","<<"x1[i]"<<","<<"Tcc"<<","<<"Calor cambra"<<","<<"Calor extreta pel fluid"<<","<<"Calor intercanviada amb l'aire exterior"<<","<<"T3[i]"<<","<<"P3[i]"<<","<<"v3[i]"<<","<<"rho3[i]"<<","<<"alfa3[i]"<<","<<"T2[i]"<<","<<"T4[i]"<<endl; //alfa, T2,T4[n-1]=0 perque no hi ha nodes
    for (int i = 0; i < n+1; i++)
    {   
        fout<<setprecision(15)<<i<<","<<x1[i]<<","<<Tcomb<<","<<Q_cambra<<","<<abs(Q_34-Q_cambra)<<","<<Q_34<<","<<T3[i]<<","<<p3[i]<<","<<v3[i]<<","<<rho3[i]<<","<<alfa3[i]<<","<<T2[i]<<","<<T4[i]<<endl;
        
    }
    
}



