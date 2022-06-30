clear;
path='Treball_sortida_n=2048_delta-10_geometria.csv';
Sortida=readmatrix(path);
x=Sortida(:,2);
t=Sortida(1,14);
x2=Sortida(1:end-1,2);
Tcc=Sortida(1,3)
Q_cambra=Sortida(1,4)
Q_fluid=Sortida(1,5)
%nom=join(['Sortida n=',n(i),' \delta=10^{-',delta(i),'}  Temps execució=',num2str(t),' s']);
Q_aire=Sortida(1,6)
T3=Sortida(:,7);
P3=Sortida(:,8);
v3=Sortida(:,9);
rho3=Sortida(:,10);
T2=Sortida(1:end-1,12);
T4=Sortida(1:end-1,13);
T2max=max(T2)
T4max=max(T4)
f=figure(Name='Grafics');
f.Position(3:4)=[1538,929];
sgtitle('Propietats fluid tub 3 canviant la geometria',Interpreter='latex')
subplot(2,2,1)
plot(x,v3)
ylabel('Velocitat m/s')
xlabel('x (m)')
title('Velocitat 3-x')
subplot(2,2,2)
plot(x,rho3)
xlabel('x (m)')
ylabel('Densitat (kg/m^3)')
title('Densitat 3-x')
subplot(2,2,3)
plot(x,T3)
xlabel('x (m)')
ylabel('T (K)')
title('Temperatura 3-x')
subplot(2,2,4)
plot(x,P3)
xlabel('x (m)')
ylabel('Pressió (Pa)')
title('Pressió 3-x')
f2=figure(Name='Temperatures i calors');
f2.Position(3:4)=[1538,929];
sgtitle('Temperatures dels tubs')
subplot(2,1,1)
plot(x2,T2)
xlabel('x (m)')
ylabel('T (K)')
title('Temperatura 2-x')
subplot(2,1,2)
plot(x2,T4)
xlabel('x (m)')
ylabel('T (K)')
title('Temperatura 4-x')