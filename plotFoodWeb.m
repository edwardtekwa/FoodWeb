%plotFoodWeb.m
%Ed Tekwa Dec 20, 2021

%load a single simulation run. Use following line for example in paper, or replace name with file generated from
%Make_warming_MultPtsstats_Parallel.m:
load('Foodweb_numSpecies200_dT3_basalSize0.01_meanD3_pInedible0_fIII.mat')

plot_demog_spatial_end; %plot time series and spatial distributions

set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',20)
scrsz = get(0,'ScreenSize');
figs(1)=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/9 scrsz(4)/2.2]);

timePt=[1,200];
patches=[8 9 10 11];

for t=timePt
    [gainB gainZ dB dZ v TE PB TLik TLi TLk TLall Cij(:,:,:,t)] = sub_demog_diagnostics(t*365,Bw_yrs(:,:,t),Zw_yrs(:,:,t),repmat(P.T,1,1) + (t*365-1)*P.dT,P);
end  

for pl=1:length(patches)
    subplot(length(patches),1,pl)
    ptch=patches(pl)
    hold on
    for t=timePt
        if t==1
            sColor=[0 0 1];
            cColor=[0.9 0.9 1];
        else
            sColor=[1 0 0];
            cColor=[1 0.9 0.9];
        end
        Zw_t=Zw_yrs(ptch,1,t);
        Bw_t=Bw_yrs(ptch,:,t);
        Bw_t(Bw_t==0)=NaN;
        SearchTopt=[P.T(ptch)+P.dT*365*t P.z];
        ZBw_t=[Zw_t Bw_t];
        for sp1=2:size(Cij,1) %predator
            for sp2=1:size(Cij,2) %prey
                if Cij(sp1,sp2,ptch,t)>0.00001 && Bw_t(sp1)>0 && ZBw_t(sp2)>0
                    %plot([P.T(ptch)+P.dT*365*t P.z],P.S,'LineWidth',Cij(sp1,sp2+1,ptch,t));
                    plot([SearchTopt(sp1+1),SearchTopt(sp2)],[P.S(sp1+1),P.S(sp2)],'Color',cColor,'LineWidth',log10(100000*Cij(sp1,sp2,ptch,t)));
                end
            end
        end
        sc=scatter(SearchTopt,P.S,ZBw_t*200,sColor,'filled');
        alpha(sc,0.3)

        
        if pl==length(patches) && t==1
            text(P.T(ptch)+P.dT*365*t-0.55,-2,'yr 800');
        elseif pl==length(patches) && t==200
            text(P.T(ptch)+P.dT*365*t-0.75,-2,'yr 1000');
        end
    end
    
    if pl==length(patches)
        sr=scatter(11,-1,1*1000,'k'); %reference
        text(10.4,-1,'1g/m^3');
    end
    
    ylim([-4 7])
    xlim([9.5 18])
    if pl==length(patches)
        xlabel('optimal search temperature [^oC]')
    end
    %ylabel('log_{10}(body size [g])')
    title(['patch ' num2str(ptch)],'Fontsize',11)
end