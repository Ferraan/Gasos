#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
using namespace std;
auto start = chrono::high_resolution_clock::now();
class fluid
{
    private:
    double Reynolds, Prandt, Nusselt;
    public:
    double Alfa_i, fregament, viscositat , densitat ,conductivitat ,cp;
    void Calcul_Coeficients(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    void Propietats_termofisiquesH2(double T, double P,double Rgas);
    void Propietats_termofisiquesO2(double T, double P,double Rgas);
    void Propietats_termofisiquesmescla(double T, double P,double Rgas,double cpA,double cpB,double fmolA,double fmolB, double muA, double muB, double conductivitatA, double conductivitatB);
    void Calcul_Coeficients_anular(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
};


const int n = 500; //Volums de control del fluid, n+1 nodes
const double delta = 1e-10;
const double pi = 2 * acos(0.0);
const double Runiversal=8.3144621;
const double massa_molarH2=2e-3, massa_molarO2=32e-3;//Kg/mol

int main(){
    std::cout.precision(15);
    std::scientific;
    double L=0.5, Di=0.02, Do=0.023, Dm=0.033, DOut=0.035;
    double Tin1=1000, pin1=161e+5, Tin3=100, cabalin3=0.005, pin3=150e+5, Text=300;
    double cabalin1H2=0.005, cabalin1O2=0.05, cabalin1Tot=cabalin1H2+cabalin1O2;
    double Ttub2_inic=800, Ttub4_inic=500;
    double rugositat2in=0.0002, rugositat2ext=0.0003, rugositat4in=0.0004, rugositat4ext=0.0002;
    
    //Calculs previs
    double Deltax=L/n, Dh=Dm-Do;
    double S1=pi*Di*Di/4, S3=pi*pow(Dm-Do,2)/4, Sl2int=pi*Di*Deltax, Sl2out=pi*Do*Deltax, Sl4in=pi*Dm*Deltax;
    double molsH2=cabalin1H2/massa_molarH2, molsO2=cabalin1O2/massa_molarO2;
    double fraccio_molarH2=molsH2/(molsH2+molsO2), fraccio_molarO2=molsO2/(molsH2+molsO2), massa_molar_cambra=fraccio_molarH2*massa_molarH2+fraccio_molarO2*massa_molarO2;
    double Rgas_cambra=Runiversal/massa_molar_cambra, Rhidrogen=Runiversal/massa_molarH2, Roxigen=Runiversal/massa_molarO2;
    double rhoin1=pin1/(Tin1*Rgas_cambra), rhoin3=pin3/(Tin3*Rhidrogen);
    double vin1=cabalin1Tot/(S1*rhoin1), vin3=cabalin3/(S3*rhoin3);
    std::vector<double> x1(n+1,0), v1(n+1,0), T1(n+1,0), p1(n+1,0), rho1(n+1,0), v3(n+1,0),T3(n+1,0), p3(n+1,0), rho3(n+1,0), T2(n,Ttub2_inic), T4(n,Ttub4_inic), alfa1(n,0), alfa3(n,0);
    x1[0]=0;
    for (int i = 1; i < n+1; i++)
    {
        x1[i]=x1[i-1]+Deltax;
    }
    
   
    //Zona 1 i 3
    v1[0]=vin1; T1[0]=Tin1; p1[0]=pin1; rho1[0]=rhoin1;
    v3[0]=vin3; T3[0]=Tin3; p3[0]=pin3; rho3[0]=rhoin3;
    fluid H2cambra, O2cambra, H2ext, mescla_cambra;
   
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
        
    }
    
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout<<"Temps execucio (s)" <<static_cast<float>(duration.count())/1000000 << endl;
    
    
    ofstream fout;
    fout.open("Treball.csv");
    fout<<"i"<<","<<"x1[i]"<<","<<"T1[i]"<<","<<"P1[i]"<<","<<"v1[i]"<<","<<"rho1[i]"<<","<<"alfa1[i]"<<","<<"T3[i]"<<","<<"P3[i]"<<","<<"v3[i]"<<","<<"rho3[i]"<<","<<"alfa3[i]"<<endl; //alfa, T2,T4[n]=0 perque no hi ha nodes
    for (int i = 0; i < n+1; i++)
    {   
        fout<<i<<","<<x1[i]<<","<<T1[i]<<","<<p1[i]<<","<<v1[i]<<","<<rho1[i]<<","<<alfa1[i]<<","<<T3[i]<<","<<p3[i]<<","<<v3[i]<<","<<rho3[i]<<","<<alfa3[i]<<endl;
        
        
    }
    
}



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