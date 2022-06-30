clear;
n=["4" "16" "32" "64" "512" "1024" "1024" "2048"];
delta=["10" "8" "8" "8" "8" "8" "10" "10"];

for i=1:8
    path=replace(join(['Treball_sortida_n=',n(i),'_delta-',delta(i),'.csv'])," ","");
    
    Sortida=readmatrix(path);
    x=Sortida(:,2);
    t=Sortida(1,14);
    x2=Sortida(1:end-1,2);
    Tcc=Sortida(1,3);
    Q_cambra=Sortida(1,4);
    Q_fluid=Sortida(1,5);
    nom=join(['Sortida n=',n(i),' \delta=10^{-',delta(i),'}  Temps execució=',num2str(t),' s']);
    Q_aire=Sortida(1,6);
    T3=Sortida(:,7);
    P3=Sortida(:,8);
    v3=Sortida(:,9);
    rho3=Sortida(:,10);
    T2=Sortida(1:end-1,12);
    T4=Sortida(1:end-1,13);
    f=figure('visible','off',Name='Grafics');
    f.Position(3:4)=[1538,929];
    sgtitle(nom)
    subplot(2,2,1)
    yyaxis left
    scatter([1,2,3],[Q_cambra,Q_fluid,Q_aire],25,"blue",'filled')
    ylabel('$$\dot{Q} (W)$$',Interpreter='latex')
    yyaxis right
    scatter([3.95],[Tcc],50,"red","filled")
    ylabel('Temperatura cambra combustió (K)')
    title('Calors extrets i temperatura cambra combustió ')
    subplot(2,2,2)
    plot(x2,T2)
    xlabel('x (m)')
    ylabel('T (K)')
    title('Temperatura 2-x')
    subplot(2,2,3)
    plot(x,T3)
    xlabel('x (m)')
    ylabel('T (K)')
    title('Temperatura 3-x')
    subplot(2,2,4)
    plot(x2,T4)
    xlabel('x (m)')
    ylabel('T (K)')
    title('Temperatura 4-x')
    saveas(f,replace(join(['n_',n(i),'_delta-',delta(i)])," ",""),'epsc')
end