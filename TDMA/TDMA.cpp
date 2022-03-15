#include <iostream>
#include <cmath>
using namespace std;
double re(double rP, double deltaR);
double rw(double rP, double deltaR);
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
    const int n = 10;
    double deltaR = (Rext - Rint) / n;
    double T[n + 4];
    double aP[n + 2];
    double aW[n + 2];
    double aE[n + 2];
    double bP[n + 2];
    double rP[n + 2];
    double AP[n+2];
    double Se[n+2];
    double Sw[n+2];
    double P[n + 4];
    double R[n + 4];
    aP[1] = 1;
    aW[1] = 0;
    aE[1] = 0;
    bP[1] = Twall;
    //Posicio nodes "reals"
    rP[0] = Rint;
    rP[1] = rP[0] + deltaR / 2;
    for (int i = 2; i <= n; i++)
    {
        rP[i] = rP[i - 1] + deltaR;
    }
    rP[n + 1] = Rext;
    //Ap
    AP[0]=0;
    for (int i = 1; i < n+1; i++)
    {
        AP[i]=2*pi*(pow(re(rP[i],deltaR),2)-pow(rw(rP[i],deltaR),2));
    }
    AP[n+1]=0;
    
    Sw[0]=0;
    Se[n+1]=0;
    for (int i = 1; i < n+2; i++)
    {
        Sw[i]=2*pi*(rP[i]-deltaR/2)*ef;
        Se[i]=2*pi*(rP[i]+deltaR/2)*ef;
    }
    aW[0]=0;
    aE[0]=0;
    bP[0]=Twall;
    aP[0]=1;
    for (int i = 1; i < n+1; i++)
    {
        aW[i]=(lambda*Sw[i])/(deltaR);
        aE[i]=(lambda*Se[i])/(deltaR);
        aP[i]=aE[i]+aW[i]+alphaext*AP[i];
        bP[i]=alphaext*Text*AP[i];
    }
    aE[n+1]=0;
    aW[n+1]=lambda/deltaR;
    aP[n+1]=aW[n+1]+alphaend;
    bP[n+1]=alphaend*Text;

}


double re(double rP, double deltaR){
        return(rP+deltaR/2);
}
double rw(double rP, double deltaR){
    return(rP-deltaR/2);
}