#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
int main(){
    double pi = 2*acos(0.0);
    float Rint=2;
    float Rext=10;
    float ef=0.5;
    float cod_term=300;
    float Twall=500;
    float Text=25;
    float alphaext=30;
    float alphaend=20;
    int n=10;
    float deltaR=(Rext-Rint)/n;
    vector<float> T;
    vector<float> ap;
    vector<float> aw;
    vector<float> ae;
    vector<float> bp;
    vector<float> rp;
    vector<float> A;
    for (int i=1; i<=n; i++){
        T.push_back(1);
        aw.push_back(1);
        ap.push_back(1);
        ae.push_back(1);
        bp.push_back(1);
        rp.push_back(1);
    }
    rp[0]=Rint;
    ap[0]=1;
    aw[0]=0;
    ae[0]=0;
    bp[0]=Twall;
    for (int i = 1; i < n; i++){
        rp[i]=rp[i-1]+deltaR;
    }
    rp.push_back(Rext);
    for (int i = 1; i < n; i++){
        A.push_back(2*pi*(pow(rp[i+1],2)-pow(rp[i-1],2)));
        cout<<A[i-1]<<endl;
    }
    
    
}