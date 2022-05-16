#include "solvers.h"
#include <vector>
#include <cmath>

int solverTDMA(std::vector<double>& T,std::vector<double> aP,std::vector<double> aW, std::vector<double> aE, std::vector<double> bP,const int n){
    std::vector<double> P(n + 4);
    std::vector<double> R(n + 4);
    P[0]=0;
    R[0]=0;
    for (int i = 1; i < n+3; i++)
    {
        P[i]=aE[i]/(aP[i]-aW[i]*P[i-1]);
        R[i]=(bP[i]+aW[i]*R[i-1])/(aP[i]-aW[i]*P[i-1]);
    }
    for (int i = n+2; i >= 1; i--)
    {
        T[i]=P[i]*T[i+1]+R[i];
    }
    
    return(1);
}