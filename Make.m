clear all; %close all;


%%%% Testing realization ofsize-shifts
global P v % P has parameters, v are current thermal envelopes
%if matlabpool('size') ~= 0; matlabpool close; end; matlabpool local 2

FoodWebFile=sprintf('Food web %s.mat', string(datetime));


%%-- Load ensemble trait values
%load ./Data/Data_traits
make_traits;

%%-- Make storage
P = make_parameters(TR,1); % get parameters of model
!rm -f ./Data/Data_ensembles.nc
nccreate('./Data/Data_ensembles.nc','B','Dimensions',...
    {'TR',size(TR,2),'nx',P.nx,'sp',P.n,'t',P.Tend+1});
nccreate('./Data/Data_ensembles.nc','Z','Dimensions',...
    {'TR',size(TR,2),'nx',P.nx,'t',P.Tend+1});


%%-- Run ensemble demographics
%for i = 200; size(TR,1);

%create matrices to store ensemble data
EnsembleGbar=[];
EnsembleGBbins=[];
EnsembleprodBbins=[];
avgB=[];
avgP=[];
Ps=P;

for i = 1:size(TR,2);
    %i=1;
    disp([num2str(i) ', basalSize=' num2str(TR{i}.s.m0) ', meanD=' num2str(HP.sdm(i))])
    
    %% make parameters for a given set of traits
    [P B Z T] = make_parameters(TR,i);
    TE = zeros(P.n, P.nx, P.Tend);
    PB = zeros(P.n, P.nx, P.Tend);
    prodZ=zeros(P.nx, 1, P.Tend);
    prodB=zeros(P.nx, P.n, P.Tend);
    
    %% Iterate model forward in time (days at the moment)
    YearStartT=1;
    for t = 1:P.Tend;
        %    for t = 1:5000;
        
        % organize
        B1      = []; B2 = [];
        B1      = B(:,:,t);Z1 = Z(:,:,t);
        
        % fix
        B1(B1<eps) = 0;% eps;
        Z1(Z1<eps) = 0;% eps;
        
        % Shift thermal gradient
        T1      = P.T;% + t.*P.dT; %<<< add this when time is right
        if t<100000 %initial period
            T1w      = P.T;
        elseif t>100000 && t<=200000 %warming period (2 degree increase total)
            T1w      = P.T + (t-100000)*P.dT;
        else %new stable period
            T1w      = P.T +2;
        end
        %T1      = P.T + 2*sin(2*pi*((365-t.*P.dt)/365)); %seasonal cycle
        %T1      = P.T + 2*sin(2*pi*((365-t.*P.dt/10)/365)); %decadal cycle
        
        
        % demographics
        B2        = sub_move(B1); % move
        [gainB gainZ dB dZ v TE(:,:,t) PB(:,:,t)] = sub_demog(1,B2,Z1,T1); % grow/die
        
        % time step
        Z(:,:,t+1) = Z1 + (dZ .* P.dt);
        B(:,:,t+1) = B2 + (dB .* P.dt);
        prodZ(:,:,t) = gainZ;
        prodB(:,:,t) = gainB;
        
        %fix again?
        %B(B(:,:,t+1)<eps) = eps;
        
        %disp(num2str(t))
    end
    
    %     %% Save ensemble member here
    %     ncwrite('./Data/Data_ensembles.nc','B',...
    %             reshape(B,[1 P.nx P.n P.Tend+1]),[i 1 1 1]);
    %     ncwrite('./Data/Data_ensembles.nc','Z',...
    %             reshape(Z,[1 P.nx P.Tend+1]),[i 1 1]);
    
    
    
    %%-- Data analysis
    %B(B<eps) = eps;
    %Z(Z<eps) = eps;
    %% Simple plotting to file (see Fig/biomass_plots.pdf)
    %plot_demog_spatial(25, Z, B, P, v) % plots ts from patch 25
    %plot_demog_spatial_B(1,Z, B, P) % plots ts from patch 1 (coldest)
    set(0,'DefaultFigureVisible', 'off'); %use this to suppress plot display but still runs through the stats
    [EnsembleGbar(i,:),EnsembleGBbins(i,:),EnsembleprodBbins(i,:),avgB(:,:,i),avgP(:,:,i)]=plot_demog_spatial_B(P.nx,Z, B, P,HP,prodZ,prodB,i);
    %plot_demog_spatial_spectra(P.nx,Z, B(:,:,ceil(4*t/5):t), P,HP,prodZ,prodB,i) %spectra plots only
    Ps(i)=P; %store parameters
    
    if min(min(min(B)))<log10(eps)
       disp('negative biomass error'); 
    end
    %save and overwrite after every new simulation in the ensemble
    save(FoodWebFile,'EnsembleGbar','EnsembleGBbins','EnsembleprodBbins','avgB','avgP','Ps');
end
%compute and plot ensemble stats
% scrsz = get(0,'ScreenSize');
% fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);
% subplot(2,3,1)
% b=bar(nanmean(EnsembleGbar)-floor(log10(eps)),'FaceColor',BiomassColor);
%     hold on
%     b=bar(Y2bins-floor(log10(eps)),'FaceColor',ProdColor);
%     errorbar([1:1+numBins]-0.04,Ybins-floor(log10(eps)),Ybins-YbinsTimeLo,YbinsTimeHi-Ybins,'.','Color',BiomassErrColor,'LineWidth',2);
%     errorbar([1:1+numBins]+0.04,Y2bins-floor(log10(eps)),Y2bins-Y2binsTimeLo,Y2binsTimeHi-Y2bins,'.','Color',ProdErrColor,'LineWidth',2);
%    
