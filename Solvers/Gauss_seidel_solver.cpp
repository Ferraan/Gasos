#include "solvers.h"
#include <vector>
#include <cmath>

int solverGS(std::vector<double>& T,std::vector<double> aP,std::vector<double> aW, std::vector<double> aE, std::vector<double> bP,const int n, double delta,double fr){
    float dif=0.5;
    double Told;
    
    while (dif>delta)
    {   
        dif=0;
        for (int i = 1; i < n+3; i++)
        {
            Told=T[i];
            T[i]=(aE[i]*T[i+1]+aW[i]*T[i-1]+bP[i])/aP[i];
            if (abs(Told-T[i])>dif)
            {
                dif=abs(Told-T[i]);
            }
            T[i]=Told+fr*(T[i]-Told);
            
            
        }    
        
    }  
    return(1);
}