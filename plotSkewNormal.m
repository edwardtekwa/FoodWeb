set(0,'defaultaxeslinewidth',1)
set(0,'DefaultAxesFontSize',20)
scrsz = get(0,'ScreenSize');
figs(1)=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/10 scrsz(4)/8]);

v_t=[];
Z_opt=14;
Temps=10:0.01:18;
for T=Temps
    v_t=[v_t skewThEnv(1,T,Z_opt)];
end
plot(Temps,v_t,'k','LineWidth',2);
xlabel '^oC'
ylabel 'v'
xlim([11 16])