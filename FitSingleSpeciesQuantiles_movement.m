
%Edward Wong Oct 22, 17
%collect data from Make_warming_MultPtsstats_Parallel0 simulations

%clear
set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',14)

numFit=5; %number of steps in quantile ranges to try
maxQuant=0.9;
minQuant=0.5;
QuantRange=[minQuant:(maxQuant-minQuant)/(numFit-1):maxQuant];
%choose temperature change scenario to display: (position of tempChange =[-2 2 4 6])
TempScenario=3; %position is without 0 change position [+2 +4 +6]
recordYrs=20; %use these number of years at the end of time series to compute mean statistics
%Cases={'basalSize0.01_meanD-Inf_stdD0';'basalSize0.01_meanD0_stdD0';'basalSize0.01_meanD3_stdD0';'basalSize0.01_meanD6_stdD0';'basalSize0.01_meanD9_stdD0';'basalSize0.01_meanD12_stdD0'};
Cases={'basalSize0.01_meanD-Inf';'basalSize0.01_meanD0';'basalSize0.01_meanD3';'basalSize0.01_meanD6';'basalSize0.01_meanD9';'basalSize0.01_meanD12'};
%Cases={'basalSize0.001_meanD-Inf_stdD0';'basalSize0.001_meanD0_stdD0';'basalSize0.001_meanD3_stdD0';'basalSize0.001_meanD6_stdD0';'basalSize0.001_meanD9_stdD0';'basalSize0.001_meanD12_stdD0'};
numCases=length(Cases);
moveRates=[-Inf,0,3,6,9,12];
%moveRates=[-Inf,0,3,6];
caseIts=[]; %record number of iterations for each case



'Select files containing simulation results in first folder'

Files={};
Paths={};

reply1 = 'a';

while (reply1 == 'a')
    [FileName,PathName]=uigetfile('.mat','Simulation data','MultiSelect','on')
    newFiles=cellstr(FileName);
    newPaths=cellstr(PathName);
    Files=[Files, newFiles];
    for f=1:length(newFiles)
        Paths=[Paths, newPaths];
    end
    reply1 = input('Finished? (enter to continue, a to add files, q to quit): ','s');
    if isempty(reply1)
        reply1 = 'continue';
    end
end

numFiles=length(Files);

if (reply1 ~='q') %continue
    for CaseNumber=1:numCases %first, check for minimum number of iterations for each case
        iteration=0;
        Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to movement rate treatment of CaseNumber
        for filePos=1:numFiles
            if(Positions(filePos)==1)
                iteration=iteration+1;
                caseIts(CaseNumber)=iteration; %update number of iterations for each case
            end
        end
    end
    numIt=min(caseIts); %minimum number of iterations for all cases
    
    for CaseNumber=1:numCases
        iteration=0;
        Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to landscape type CaseNumber
        for filePos=1:numFiles
            if(Positions(filePos)==1 && iteration<numIt)
                iteration=iteration+1;
                load([Paths{filePos} Files{filePos}]);
                
                %add data to aggregate matrices
                %parameters: (record the case parameters at indices in each variable that follows)
                findPos1=findstr(Cases{CaseNumber},'basalSize')+length('basalSize');
                findPos2=findstr(Cases{CaseNumber},'_meanD')-1;
                findPos3=findstr(Cases{CaseNumber},'_meanD')+length('_meanD');
                findPos4=findstr(Cases{CaseNumber},'_pInedible')-1;
                %findPos4=findstr(Cases{CaseNumber},'_stdD')-1;
                %findPos5=findstr(Cases{CaseNumber},'_stdD')+length('_stdD');
                
                %Redo single-species LV hindcast fits by finding the best biomass
                %and production quantiles to match model projections to
                %forecast:
                for i=1:length(QuantRange)
                    for j=1:length(QuantRange)
                    [r4,a4,z4,K4,flag4,raR24,r_T4,K_T4,K_T_ratio4,r_T_ratio4]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(4,:)); %fit growth model to all patches at once
                    end
                end
            end
        end
    end
end
