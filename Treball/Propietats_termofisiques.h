class fluid
{
    private:
    double Reynolds, Prandt, Nusselt, Grashof, Rayleigh;
    public:
    double Alfa_i, fregament, viscositat , densitat ,conductivitat ,cp , beta;
    void Calcul_var_adim(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    void Propietats_termofisiquesH2(double T0, double Tf, double P,double Rgas);
    void Propietats_termofisiquesO2(double T0,double Tf, double P,double Rgas);
    void Propietats_termofisiquesH2O(double T0,double Tf, double P,double Rgas);
    void Calcul_Coeficients_anular(double mu, double cp, double lambda, double D, double v, double rho, double rugositat_relativa);
    void Propietats_termofisiquesaire(double T, double P,double Rgas);
    void Calcul_coeficients_exterior(double cp, double lambda, double mu, double g, double beta, double rho, double Tm, double Tf, double X);
    void Propietats_termofisiquescambra(double T0, double P,double Rgas,double v_H2O,double v_H2exces, double v_O2exces);
};
double condMolibde(double T);