#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
struct Propietats_termo
{
    double viscositat,densitat,conductivitat,cp;
};

class Coeficients
{
    private:
    double Reynolds, Prandt, Nusselt;
    public:
    double Alfa_i, fregament;
    void Calcul_Coeficients(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    
};


Propietats_termo Propietats_termofisiques(double T, double P,double Rgas_cambra);

const int n = 1000;
const double delta = 1e-10;
const double pi = 2 * acos(0.0);
const double Runiversal=8.3144621;
const double massa_molarH2=1e-3, massa_molarO2=16e-3;//Kg/mol

int main(){
    std::cout.precision(15);
    std::scientific;
    double L=0.5, Di=0.02, Do=0.023, Dm=0.03, DOut=0.033;
    double Tin1=1000, cabalin1H2=41.2, cabalin1O2=250-cabalin1H2, cabalin1Tot=cabalin1H2+cabalin1O2, pin1=161e+5, Tin3=100, cabalin3=100, pin3=150e+5, Text=300;
    double Ttub2_inic=800, Ttub4_inic=500;
    double rugositat2in=0.0002, rugositat2ext=0.0003, rugositat4in=0.0004, rugositat4ext=0.0002;
    
    //Calculs previs
    double Deltax=L/n;
    double S1=pi*Di*Di/4, S3=pi*pow(Dm-Do,2)/4;
    double molsH2=cabalin1H2/massa_molarH2, molsO2=cabalin1O2/massa_molarO2;
    double fraccio_molarH2=molsH2/(molsH2+molsO2), fraccio_molarO2=molsO2/(molsH2+molsO2), massa_molar_cambra=fraccio_molarH2*massa_molarH2+fraccio_molarO2*massa_molarO2;
    double Rgas_cambra=Runiversal/massa_molar_cambra, Rhidrogen=Runiversal/massa_molarH2;
    double rhoin1=pin1/(Tin1*Rgas_cambra), rhoin3=pin3/(Tin3*Rhidrogen);
    double vin1=cabalin1Tot/(S1*rhoin1), vin3=cabalin3/(S3*rhoin3);
    std::vector<double> x1(n), v1(n), T1(n), p1(n), rho1(n), v3(n),T3(n), p3(n), rho3(n), T2(n-1,Ttub2_inic), T4(n-1,Ttub4_inic);
    x1[0]=0;
    for (int i = 1; i < n; i++)
    {
        x1[i]=x1[i-1]+Deltax;
    }
    
    //Zona 1
    v1[0]=vin1; T1[0]=Tin1; p1[0]=pin1; rho1[0]=rhoin1;
    Propietats_termo Propietats_H2, Propietats_H2_wall;
    Coeficients Coeficients_cambra;

    for (int i = 1; i < n-1; i++)
    {   
        v1[i]=v1[i-1]; p1[i]=p1[i-1]; T1[i]=T1[i-1]; rho1[i]=rho1[i-1];
        double error=1;
        while (error>delta)
        {
            
            
            double Ti=(T1[i]+T1[i-1])/2, Pi=(p1[i-1]+p1[i])/2,vi=(v1[i]+v1[i-1])/2,rhoi=(rho1[i]+rho1[i-1])/2;
            Propietats_H2=Propietats_termofisiques(Ti,Pi,Rgas_cambra); 
            // Propietats_H2_wall=Propietats_termofisiques(T2[i],Pi,Rgas_cambra); No necessitem si es flux turbulent de gasos
            Coeficients_cambra.Calcul_Coeficients(Propietats_H2.viscositat,Propietats_H2.cp,Propietats_H2.conductivitat,Di,vi,Propietats_H2.densitat,rugositat2in);
            p1[i]=-cabalin1Tot*(v1[i]-v1[i-1])+p1[i-1]*S1-Coeficients_cambra.fregament*rhoi*pow(vi,2)/2*pi*Di*Deltax;

        }
        
    }
    
    
}

Propietats_termo Propietats_termofisiques(double T, double P,double R){
    Propietats_termo Temp;
    Temp.densitat=P/(R*T);
    Temp.viscositat= (5.83673e-9*pow(T,3)-1.7564449e-5*pow(T,2)+0.030645686887*T+1.113450152679)*pow(10,-6);
    Temp.cp=(-1.468867e-11*pow(T,4)+5.0109515e-8*pow(T,3)-5.914804e-5*pow(T,2)+2.918448e-2*T+9.472783)*pow(10,3);
    Temp.conductivitat=3.69347e-15*pow(T,4)+6.18450e-11*pow(T,3)-3.08777e-7*pow(T,2)+6.8004e-4*T+2.8372e-3;
    return(Temp);
}

void Coeficients::Calcul_Coeficients(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa)
{
    Reynolds=(rho*v*D)/mu;
    Prandt=mu*cp/lambda;
    Nusselt=0.023*pow(Reynolds,0.8)*pow(Prandt,0.4);
    //Alfa_i=lambda*Nusselt/
    double A,B;
    //Churchill expression
    A=pow(2.457*log(1/(pow(7/Reynolds,0.9)+0.27*rugositat_relativa)),16);
    B=pow(37530/Reynolds,16);
    fregament=2*pow(pow(8/Reynolds,12)+1/pow(A+B,1.5),1.0/12.0);
}