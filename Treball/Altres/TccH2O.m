clear;
h_fch4=-74850;
W_ch4=0.01604;
h_fCO2=-393520;
h_fH2O=-241820;
T_fin=300;
TAin=300;
T0=298;

vH2O=1;
v_H2=1;
vO2=0.5;

a0h2o=0.4198640560e+01 ;
a1h2o=-0.2036434100e-02 ;
a2h2o=0.6520402110e-05 ;
a3h2o=-0.5487970620e-08 ;
a4h2o=0.1771978170e-11;
b0h2o=0.3033992490e+01;
b1h2o=0.2176918040e-02 ;
b2h2o=-0.1640725180e-06 ;
b3h2o=-0.9704198700e-10 ;
b4h2o=0.1682009920e-13;

a0o2=0.3782456360e+01 ;
a1o2=-0.2996734160e-02 ;
a2o2=0.9847302010e-05 ;
a3o2=-0.9681295090e-08 ;
a4o2=0.3243728370e-11;
b0o2=0.3282537840e+01 ;
b1o2=0.1483087540e-02 ;
b2o2=-0.7579666690e-06 ;
b3o2=0.2094705550e-09 ;
b4o2=-0.2167177940e-13;

a0n2=0.3298677000e+01;
a1n2=0.1408240400e-02 ;
a2n2=-0.3963222000e-05 ;
a3n2=0.5641515000e-08 ;
a4n2=-0.2444854000e-11;
b0n2=0.2926640000e+01;
b1n2=0.1487976800e-02 ;
b2n2=-0.5684760000e-06 ;
b3n2=0.1009703800e-09 ;
b4n2=-0.6753351000e-14;

a0ch4=0.5149876130e+01 ;
a1ch4=-0.1367097880e-01 ;
a2ch4=0.4918005990e-04 ;
a3ch4=-0.4847430260e-07; 
a4ch4=0.1666939560e-10;
b0ch4=0.7485149500e-01;
b1ch4=0.1339094670e-01 ;
b2ch4=-0.5732858090e-05 ;
b3ch4=0.1222925350e-08 ;
b4ch4=-0.1018152300e-12;

a0co2=0.2356773520e+01;
a1co2=0.8984596770e-02 ;
a2co2=-0.7123562690e-05 ;
a3co2=0.2459190220e-08 ;
a4co2=-0.1436995480e-12;
b0co2=0.3857460290e+01 ;
b1co2=0.4414370260e-02 ;
b2co2=-0.2214814040e-05 ;
b3co2=0.5234901880e-09 ;
b4co2=-0.4720841640e-13;

a0h2=0.2344331120e+01 ;
a1h2=0.7980520750e-02 ;
a2h2=-0.1947815100e-04 ;
a3h2=0.2015720940e-07 ;
a4h2=-0.7376117610e-11;
b0h2=0.3337279200e+01;
b1h2=-0.4940247310e-04 ;
b2h2=0.4994567780e-06 ;
b3h2=-0.1795663940e-09 ;
b4h2=0.2002553760e-13;


terme2=cp(T0,T_fin,a0h2,a1h2,a2h2,a3h2,a4h2,b0h2,b1h2,b2h2,b3h2,b4h2)*(T_fin-T0)+vO2*cp(T0,TAin,a0o2,a1o2,a2o2,a3o2,a4o2,b0o2,b1o2,b2o2,b3o2,b4o2)*(TAin-T0);
Tcc=3000;
error=1;

terme1=(-(vH2O*h_fH2O))
H2in=cp(T0,TAin,a0h2,a1h2,a2h2,a3h2,a4h2,b0h2,b1h2,b2h2,b3h2,b4h2)*(TAin-T0)
O2in=vO2*cp(T0,TAin,a0o2,a1o2,a2o2,a3o2,a4o2,b0o2,b1o2,b2o2,b3o2,b4o2)*(TAin-T0)
H2Oout=vH2O*cp(298,Tcc,a0h2o,a1h2o,a2h2o,a3h2o,a4h2o,b0h2o,b1h2o,b2h2o,b3h2o,b4h2o)*(Tcc-T0)
while (error>10^-10)
    Tccold=Tcc;
    Tcc=T0+(terme1+terme2)/(+vH2O*cp(T0,Tcc,a0h2o,a1h2o,a2h2o,a3h2o,a4h2o,b0h2o,b1h2o,b2h2o,b3h2o,b4h2o));
    Tcc=Tccold+0.3*(Tcc-Tccold);
    error=abs(Tcc-Tccold);
end
% terme1=(h_fch4-(v_co2*h_fCO2+vH2O*h_fH2O))
% Co2out=v_co2*cp(298,Tcc,a0co2,a1co2,a2co2,a3co2,a4co2,b0co2,b1co2,b2co2,b3co2,b4co2)*(Tcc-T0)
% H2Oout=vH2O*cp(298,Tcc,a0h2o,a1h2o,a2h2o,a3h2o,a4h2o,b0h2o,b1h2o,b2h2o,b3h2o,b4h2o)*(Tcc-T0)
% n2out=v_N2*cp(T0,Tcc,a0n2,a1n2,a2n2,a3n2,a4n2,b0n2,b1n2,b2n2,b3n2,b4n2)*(Tcc-T0)
% ch4in=cp(T0,T_fin,a0ch4,a1ch4,a2ch4,a3ch4,a4ch4,b0ch4,b1ch4,b2ch4,b3ch4,b4ch4)*(T_fin-T0)
% O2in=vO2*cp(T0,TAin,a0o2,a1o2,a2o2,a3o2,a4o2,b0o2,b1o2,b2o2,b3o2,b4o2)*(TAin-T0)
% n2in=v_N2*cp(T0,TAin,a0n2,a1n2,a2n2,a3n2,a4n2,b0n2,b1n2,b2n2,b3n2,b4n2)*(TAin-T0)
% sum=-terme1-ch4in-O2in-n2in+n2out+Co2out+H2Oout
