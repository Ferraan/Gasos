#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>
using namespace std;

auto start = chrono::high_resolution_clock::now();

const int n = 100;
double re(double rP, double deltaR);
double rw(double rP, double deltaR);
void solverGS(double (&T)[n+3],double (&T_old)[n+3],double (&T_new)[n+3],double (&T_calc)[n+3], double aP[], double aW[], double aE[], double bP[],const int n,double delta,double fr);

int main()
{
    double pi = 2 * acos(0.0);
    double Rint = 2;
    double Rext = 10;
    double ef = 0.5;
    double lambda = 300;
    double Twall = 500;
    double Text = 25;
    double alphaext = 30;
    double alphaend = 20;
    double delta = 1e-10;
    double deltaR = (Rext - Rint) / n;
    double T[n + 3];
    double T_new[n + 3];
    double T_old[n + 3];
    double T_calc[n + 3];
    fill_n(T_new, n+3, 100);//T inicial
    fill_n(T_old, n+3, 100);
    fill_n(T, n+3, 100);
    double fr=1;
    double aP[n + 3];
    double aW[n + 3];
    double aE[n + 3];
    double bP[n + 3];
    double rP[n + 3];
    double AP[n+3];
    double Se[n+3];
    double Sw[n+3];
    aP[1] = 1;
    bP[1] = Twall;
    //Posicio nodes "reals"
    rP[0] = 0;
    T[0]=0;
    Se[0]=0;
    rP[1] = Rint;
    rP[2] = rP[1] + deltaR / 2;
    for (int i = 3; i < n+2; i++)
    {
        rP[i] = rP[i - 1] + deltaR;
    }
    rP[n + 2] = Rext;
    //Ap
    AP[0]=0;
    AP[1]=0;
    for (int i = 2; i < n+2; i++)
    {
        AP[i]=2*pi*(pow(re(rP[i],deltaR),2)-pow(rw(rP[i],deltaR),2));
    }
    AP[n+2]=0;
    
    Sw[0]=0;
    Sw[1]=0;
    Se[n+2]=0;
    for (int i = 1; i < n+3; i++)
    {
        Sw[i]=2*pi*(rP[i]-(rP[i]-rP[i-1])/2)*ef;
        Se[i]=2*pi*(rP[i]+(rP[i+1]-rP[i])/2)*ef;
    }
    aW[0]=0;
    aE[0]=0;
    bP[0]=0;
    aP[0]=0;
    aW[1]=0;
    aE[1]=0;
    bP[1]=Twall;
    aP[1]=1;
    for (int i = 2; i < n+2; i++)
    {
        aW[i]=(lambda*Sw[i])/(abs(rP[i-1]-rP[i]));
        aE[i]=(lambda*Se[i])/(rP[i]-rP[i-1]);
        aP[i]=aE[i]+aW[i]+alphaext*AP[i];
        bP[i]=alphaext*Text*AP[i];
    }
    aE[n+2]=0;
    aW[n+2]=(lambda)/(rP[n+2]-rP[n+1]);
    aP[n+2]=aW[n+2]+alphaend;
    bP[n+2]=alphaend*Text;

    solverGS(T,T_old,T_new,T_calc,aP,aW,aE,bP,n,delta,fr);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    cout<<"i"<<setw(17)<<"rP[i]"<<setw(17)<<"T[i]"<<setw(17)<<"Se[i]"<<setw(17)<<"Sw[i]"<<setw(17)<<"aE[i]"<<setw(17)<<"aW[i]"<<setw(17)<<"aP[i]"<<setw(17)<<"bP[i]"<<endl;
    for (int i = 0; i < n+3; i++)
    {   
        cout<<setw(17)<<left<<i<<setw(17)<<left<<rP[i]<<setw(17)<<left<<T[i]<<setw(17)<<left<<Se[i]<<setw(17)<<left<<Sw[i]<<setw(17)<<left<<aE[i]<<setw(17)<<left<<aW[i]<<setw(17)<<left<<aP[i]<<setw(17)<<left<<bP[i]<<endl;
        
    }
    cout<<"Temps execucio (s)" <<static_cast<float>(duration.count())/1000000 << endl;
}


double re(double rP, double deltaR){
        return(rP+deltaR/2);
}
double rw(double rP, double deltaR){
    return(rP-deltaR/2);
}
void solverGS(double (&T)[n+3],double (&T_old)[n+3],double (&T_new)[n+3],double (&T_calc)[n+3], double aP[], double aW[], double aE[], double bP[],const int n,double delta,double fr){
    float dif=0.5;
    while (dif>delta)
    {   
        dif=0;
        for (int i = 1; i < n+3; i++)
        {            
            T_calc[i]=(aE[i]*T_new[i+1]+aW[i]*T_new[i-1]+bP[i])/aP[i];
            
        }
        for (int i = 1; i < n+3; i++){
            T_new[i]=T_old[i]+fr*(T_calc[i]-T_old[i]);
            if(T_old[i]-T_new[i]>dif){
                dif=abs(T_old[i]-T_new[i]);
            }
        }
        copy(begin(T_new),end(T_new),begin(T_old));//Told=Tnew
    }   
    copy(begin(T_calc),end(T_calc),begin(T));
}