#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
int main(){
    double pi = 2*acos(0.0);
    double SE=0;
    double SW=0;
    float Rint=2;
    float Rext=10;
    float ef=0.5;
    float lambda=300;
    float Twall=500;
    float Text=25;
    float alphaext=30;
    float alphaend=20;
    int n=10;
    float deltaR=(Rext-Rint)/n;
    vector<float> T;
    vector<float> aP;
    vector<float> aW;
    vector<float> aE;
    vector<float> bP;
    vector<float> rP;
    vector<float> A;
    vector<float> P;
    vector<float> R;
    for (int i=1; i<=n; i++){
        T.push_back(1);
        aW.push_back(1);
        aP.push_back(1);
        aE.push_back(1);
        bP.push_back(1);
        rP.push_back(1);
    }
    rP[0]=Rint;
    aP[0]=1;
    aW[0]=0;
    aE[0]=0;
    bP[0]=Twall;
    for (int i = 1; i < n; i++){
        rP[i]=rP[i-1]+deltaR;
    }
    rP.push_back(Rext);
    A.push_back(0);
    for (int i = 1; i < n; i++){
        A.push_back(2*pi*(pow(rP[i+1]-(deltaR/2),2)-pow(rP[i-1]+(deltaR/2),2)));
        SE=2*pi*(rP[i+1]-deltaR/2);
        SW=2*pi*(rP[i-1]+deltaR/2);
        aW[i]=lambda*SW/deltaR;
        aE[i]=lambda*SE/deltaR;
        aP[i]=aW[i]+aE[i]+alphaext*A[i];
        bP[i]=alphaext+Text*A[i];
    }
    aW.push_back(lambda*pow(rP[size(rP)-1],2)*pi/deltaR);
    bP.push_back(alphaend*Text+2*pi*Rext*ef+alphaext*Text*(pow(Rext,2)-pow(Rext-deltaR/2,2))*pi);
    aE.push_back(0);
    aP.push_back(aW[size(aW)-1]+alphaend*2*pi*Rext*ef+alphaext*(pow(Rext,2)-pow(Rext-deltaR/2,2))*pi);
    P.push_back(aE[0]/aP[0]);
    R.push_back(bP[0]/aP[0]);
    for (int i = 1; i <= n+1; i++)
    {
        P.push_back(aE[i]/(aP[i]-aW[i]*P[i-1]));
        R.push_back((bP[i]+aW[i]*R[i-1])/(aP[i]-aW[i]*P[i-1]));
    }
    T.push_back(0);
    for (int i = n+1; i > 0; i--)
    {
        T[i]=P[i]*T[i+1]+R[i];
        cout<<T[i]<<endl;
    }
    
    
}