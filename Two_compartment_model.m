%Two-compartment_model.m
%Edward Tekwa Aug 24, 18
%plot two-compartment, chemostat-resource, type-II heterotroph equilibria

%environmental characteristics:
syms T; %temperature (C)

%resource characteristics:
syms F; %chemostatic growth (dilution rate) (0.0075)
syms B0max; %maximum resource biomass (5)

%heterotroph characteristics:
syms s; %size (g)
syms lambda; %consumption efficiency (0.4)
v=1.52*s^(-0.264); %maximum search rate
D=1.34E9*exp((-7310)/(T+273)+s^0.13)*s^(-0.29); %metabolic rate
Pmax=1.01E9*exp((-7310)/(T+273))*s^(-0.239); %maximum growth rate
%Pmax=1.01E9*exp((-7310)/(12.5+273))*s^(-0.239); %maximum growth rate

%solutions:
B0=D*(1+D/Pmax)/(v*lambda); %equilibrium resource biomass
Bh=F*(lambda*B0max/D-(D/Pmax+1)/v); %equilibrium heterotroph biomass
Phhalf=lambda*F*B0max/2-D;
Phmax=lambda*F*B0max-D;

%parameterize:
F=0.0075;
B0max=5;
lambda=0.4;


scrsz = get(0,'ScreenSize');
figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.5 scrsz(4)/2.5]);
subplot(1,2,1)
yyaxis left
s=(457*0.01)^1.32 %optimal consumer size of 0.01g resource
B0num=eval(B0);
Bhnum=eval(Bh);
Brationum=eval(Bh/B0);
Phhalfnum=eval(Phhalf);
Phmaxnum=eval(Phmax);
fplot(B0num,[0 30],'b','linewidth',2);
hold on
fplot(Bhnum,[0 30],'r','linewidth',2);
%fplot(Brationum,[0 30],'k','linewidth',2);
ylim([0 2.5])
title ''
ylabel 'equilibrium biomass [g/m^3]'
xlabel 'temperature [\circC]'

yyaxis right
fplot(Phmaxnum,[0 30],'k','linewidth',2);

subplot(1,2,2)
s=10*(457*0.01)^1.32 %10x optimal consumer size of 0.01g resource
B0num=eval(B0);
Bhnum=eval(Bh);
Brationum=eval(Bh/B0);
fplot(B0num,[0 30],'b','linewidth',2);
hold on
fplot(Bhnum,[0 30],'r','linewidth',2);
%fplot(Brationum,[0 30],'k','linewidth',2);
ylim([0 2.5])
title ''
ylabel 'equilibrium biomass [g/m^3]'
xlabel 'temperature [\circC]'