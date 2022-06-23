class fluid
{
    private:
    double Reynolds, Prandt, Nusselt, Grashof, Rayleigh;
    public:
    double Alfa_i, fregament, viscositat , densitat ,conductivitat ,cp , beta;
    void Calcul_Coeficients(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    void Propietats_termofisiquesH2(double T0, double Tf, double P,double Rgas);
    void Propietats_termofisiquesO2(double T, double P,double Rgas);
    void Calcul_Coeficients_anular(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    void Propietats_termofisiquesaire(double T, double P,double Rgas);
    void Calcul_coeficients_exterior(double cp, double lambda, double mu, double g, double beta, double rho, double Tm, double Tf, double X);
    double cp_c_avg_comb(double T0, double Tf);
};
double condMolibde(double T);