#include <math.h>
#include "Propietats_termofisiques.h"
#include "Combustio.h"
const double massa_molarH2=0.2015939951e+01, massa_molarO2=0.3199880028e+02;//Kg/kmol
double Tcc(double fraccio_molarO2, double fraccio_molarH2, double cabalH2,double cabalO2){
    double entalpia_f_O2=0, entalpia_f_H2=0, entalpia_f_H2O=-241.820; //KJ/Kmol
    double cabal_molarH2=cabalH2*massa_molarH2, cabal_molarO2=cabalO2*massa_molarO2; //kmol/s
    //Comprovem si es fuel rich o oxidizer rich
    if(cabal_molarO2/2>cabal_molarH2){ //Oxidizer rich
        
    }
    else if(cabal_molarO2/2<cabal_molarH2) //Fuel rich
    {
        /* code */
    }
    else //Estequiometric
    {

    }
    
}