clear
a0h2=0.2344331120e+01 ;
a1h2=0.7980520750e-02 ;
a2h2=-0.1947815100e-04 ;
a3h2=0.2015720940e-07 ;
a4h2=-0.7376117610e-11;
a0h21=0.2344331120e+01 ;
a1h21=0.7980520750e-02 ;
a2h21=-0.1947815100e-04 ;
a3h21=0.2015720940e-07 ;
a4h21=-0.7376117610e-11;
a0h2o=0.4198640560e+01 ;
a1h2o=-0.2036434100e-02 ;
a2h2o=0.6520402110e-05 ;
a3h2o=-0.5487970620e-08 ;
a4h2o=0.1771978170e-11;
a0o2=0.3782456360e+01 ;
a1o2=-0.2996734160e-02 ;
a2o2=0.9847302010e-05 ;
a3o2=-0.9681295090e-08 ;
a4o2=0.3243728370e-11;
cabalin1H2=41.2;
cabalin1O2=208.8;
entalpia_f_H2O=-241820;
massa_molarH2=0.2015939951e-02;
massa_molarO2=0.3199880028e-1;
T0=298;
TinH2=500;
TinO2=500;
tic

cabal_molarH2=cabalin1H2/massa_molarH2;
cabal_molarO2=cabalin1O2/massa_molarO2;
cabal_molarH20=cabal_molarO2/2;
cabal_molarH2_excess=cabal_molarH2-cabal_molarH20;
v_H2=1; 
v_O2=0.5;%cabal_molarO2/cabal_molarH2; 
v_H20=1;%;cabal_molarH20/cabal_molarH2;
v_H2excess=cabal_molarH2_excess/cabal_molarH2; 
cph2=cp(T0,TinH2,a0h2,a1h2,a2h2,a3h2,a4h2);
cpo2=cp(T0,TinO2,a0o2,a1o2,a2o2,a3o2,a4o2);
hh2r=cph2*(TinH2-T0);
ho2r=cpo2*(TinO2-T0);
Tcc=3000;
error=1;
num=v_H2*hh2r+v_O2*ho2r-v_H20*entalpia_f_H2O;
while (error>10^-10)
  Tccold=Tcc;
  Tcc=(T0+num/(+v_H20*cp(T0,Tcc,a0h2o,a1h2o,a2h2o,a3h2o,a4h2o)));%v_H2excess*cp(T0,Tcc,a0h2,a1h2,a2h2,a3h2,a4h2)
  Tcc=Tccold+0.3*(Tcc-Tccold);
  error=abs(Tcc-Tccold);  
  disp(Tcc);
end
toc