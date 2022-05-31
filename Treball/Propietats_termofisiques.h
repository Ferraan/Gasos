class fluid
{
    private:
    double Reynolds, Prandt, Nusselt, Grashof, Rayleigh;
    public:
    double Alfa_i, fregament, viscositat , densitat ,conductivitat ,cp , beta;
    void Calcul_Coeficients(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    void Propietats_termofisiquesH2(double T, double P,double Rgas);
    void Propietats_termofisiquesO2(double T, double P,double Rgas);
    void Propietats_termofisiquesmescla(double T, double P,double Rgas,double cpA,double cpB,double fmolA,double fmolB, double muA, double muB, double conductivitatA, double conductivitatB);
    void Calcul_Coeficients_anular(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    void Propietats_termofisiquesaire(double T, double P,double Rgas);
    void Calcul_coeficients_exterior(double cp, double lambda, double mu, double g, double beta, double rho, double Tm, double Tf, double X);
};
double condMolibde(double T);