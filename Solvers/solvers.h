#include <vector>
#include <cmath>

int solverGS(std::vector<double>& T,std::vector<double> aP,std::vector<double> aW, std::vector<double> aE, std::vector<double> bP,const int n,double delta,double fr);
int solverTDMA(std::vector<double>& T,std::vector<double> aP,std::vector<double> aW, std::vector<double> aE, std::vector<double> bP, const int n);