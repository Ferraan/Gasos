#include "solvers.h"
#include <vector>
#include <cmath>
#include <iostream>
int solverTDMA(std::vector<double>& T,std::vector<double> aP,std::vector<double> aW, std::vector<double> aE, std::vector<double> bP,const int n){
    std::vector<double> P(n);
    std::vector<double> R(n);
    for (int i = 0; i < n; i++)
    {
        P[i]=aE[i]/(aP[i]-aW[i]*P[i-1]);
        R[i]=(bP[i]+aW[i]*R[i-1])/(aP[i]-aW[i]*P[i-1]);
    }
    for (int i = n-1; i >= 0; i--)
    {
        T[i]=P[i]*T[i+1]+R[i];
    }
    
    return(1);
}