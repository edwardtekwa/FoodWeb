%Edward Wong Oct 22, 17
%collect data from Make_warming_endstats simulations

clear
Cases={'basalSize0.01_meanD-Inf_stdD0';'basalSize0.01_meanD0_stdD0';'basalSize0.01_meanD1_stdD0';'basalSize0.01_meanD10_stdD0'};
numCases=length(Cases);

%parameters: (record the case parameters at indices in each variable that follows)
meanDispersals=[];
stdDispersals=[];
basalSizes=[];
paramIndices=[];
warmingIndices=[];
%variables for no warming/warming (col1, col2 unless otherwise commented) cases:
%global (for heterotrophs)
totBiomass=[]; %total biomass
totProd=[]; %total productivity
totRich=[]; %species richness
totBeta=[]; %beta diversity
totTrophicLevel=[]; %mean global trophic level
maxTrophicLevel=[]; %maximum global non-extinct trophic species
%by space
varBiomass=[];
varProd=[];
varTrophicLevel=[];
%by body mass
allBody_Biomass=[]; %record all non-extinct species. col1: body mass, col2: biomass for no warming cases, col3: biomass for warming cases.
allBody_Prod=[]; %record all non-extinct species. col1: body mass, col2: productivity for no warming cases, col3: productivity for warming cases.



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
    for CaseNumber=1:numCases
        iteration=0;
        Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to landscape type CaseNumber
        for filePos=1:numFiles
            if(Positions(filePos)==1)
                iteration=iteration+1;
                load([Paths{filePos} Files{filePos}]);
                
                %add data to aggregate matrices
                %parameters: (record the case parameters at indices in each variable that follows)
                findPos1=findstr(Cases{CaseNumber},'basalSize')+length('basalSize');
                findPos2=findstr(Cases{CaseNumber},'_meanD')-1;
                findPos3=findstr(Cases{CaseNumber},'_meanD')+length('_meanD');
                findPos4=findstr(Cases{CaseNumber},'_stdD')-1;
                findPos5=findstr(Cases{CaseNumber},'_stdD')+length('_stdD');
                basalSizes=[basalSizes;str2num(Cases{CaseNumber}(findPos1:findPos2))];
                meanDispersals=[meanDispersals;str2num(Cases{CaseNumber}(findPos3:findPos4))];
                stdDispersals=[stdDispersals;str2num(Cases{CaseNumber}(findPos5:end))];
                paramIndices=[paramIndices;CaseNumber];
                warmingIndices=[warmingIndices;[0 1]];
        
                %variables for no warming/warming (col1, col2 unless otherwise commented) cases:
                %global (for heterotrophs)
                totBiomass=[totBiomass; [sum(B(:)) sum(Bw(:))]]; %total biomass
                totProd=[totProd; [sum(gainB(:)) sum(gainBw(:))]]; %total productivity
                totRich=[totRich; [sum(sum(B)>eps) sum(sum(Bw)>eps)]]; %species richness
                totBeta=[totBeta; [sum(sum(B)>eps)/nanmean(sum(B>eps,2)) sum(sum(Bw)>eps)/nanmean(sum(Bw>eps,2))]]; %beta diversity
                totTrophicLevel=[totTrophicLevel; [TLall TLallw]]; %mean global trophic level
                maxTrophicLevel=[maxTrophicLevel; [max(TLi) max(TLiw)]]; %maximum non-extinct global trophic species
            end
        end
    end
end

scrsz = get(0,'ScreenSize');
figs(1)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
subplot(3,2,1)
%gscatter([paramIndices;paramIndices],totBiomass(:),warmingIndices(:),'br','o')
%violin(reshape(totBiomass,20,8));
biomFit=fitlm(paramIndices,totBiomass(:,1),'quadratic');
hold on
biomWFit=fitlm(paramIndices,totBiomass(:,2),'quadratic');
biomPlot=plot(biomFit);
set(biomPlot,'color','b','linewidth',2)
biomWPlot=plot(biomWFit);
set(biomWPlot,'color','r','linewidth',2)
title ''
xlabel 'movement rate'
ylabel 'biomass'
legend off
subplot(3,2,2)
hold on
%gscatter([paramIndices;paramIndices],totProd(:),warmingIndices(:),'br','o')
prodFit=fitlm(paramIndices,totProd(:,1),'quadratic');
prodWFit=fitlm(paramIndices,totProd(:,2),'quadratic');
prodPlot=plot(prodFit);
set(prodPlot,'color','b','linewidth',2)
prodWPlot=plot(prodWFit);
set(prodWPlot,'color','r','linewidth',2)
title ''
xlabel 'movement rate'
ylabel 'productivity'
legend off
subplot(3,2,3)
%gscatter([paramIndices;paramIndices],totRich(:),warmingIndices(:),'br','o')
hold on
richFit=fitlm(paramIndices,totRich(:,1),'quadratic');
richWFit=fitlm(paramIndices,totRich(:,2),'quadratic');
richPlot=plot(richFit);
set(richPlot,'color','b','linewidth',2)
richWPlot=plot(richWFit);
set(richWPlot,'color','r','linewidth',2)
xlabel 'movement rate'
ylabel 'richness'
title ''
legend off
subplot(3,2,4)
%gscatter([paramIndices;paramIndices],totBeta(:),warmingIndices(:),'br','o')
hold on
betaFit=fitlm(paramIndices,totBeta(:,1),'quadratic');
betaWFit=fitlm(paramIndices,totBeta(:,2),'quadratic');
betaPlot=plot(betaFit);
set(betaPlot,'color','b','linewidth',2)
betaWPlot=plot(betaWFit);
set(betaWPlot,'color','r','linewidth',2)
xlabel 'movement rate'
ylabel 'beta diversity'
title ''
legend off
subplot(3,2,5)
%gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
hold on
trophFit=fitlm(paramIndices,totTrophicLevel(:,1),'quadratic');
trophWFit=fitlm(paramIndices,totTrophicLevel(:,2),'quadratic');
trophPlot=plot(trophFit);
set(trophPlot,'color','b','linewidth',2)
trophWPlot=plot(trophWFit);
set(trophWPlot,'color','r','linewidth',2)
xlabel 'movement rate'
ylabel 'mean trophic level'
title ''
legend off
subplot(3,2,6)
%gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
hold on
trophFit=fitlm(paramIndices,maxTrophicLevel(:,1),'quadratic');
trophWFit=fitlm(paramIndices,maxTrophicLevel(:,2),'quadratic');
trophPlot=plot(trophFit);
set(trophPlot,'color','b','linewidth',2)
trophWPlot=plot(trophWFit);
set(trophWPlot,'color','r','linewidth',2)
xlabel 'movement rate'
ylabel 'max trophic level'
title ''
legend off

figs(2)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
subplot(3,2,1)
violin(reshape(log2(totBiomass(:,2)./totBiomass(:,1)),20,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
hold on
refl=refline(0,0);
set(refl,'color','k')
xlabel 'movement rate'
ylabel 'log ratio biomass'
subplot(3,2,2)
violin(reshape((totProd(:,2)-totProd(:,1)),20,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
hold on
refl=refline(0,0);
set(refl,'color','k')
xlabel 'movement rate'
ylabel 'difference productivity'
subplot(3,2,3)
violin(reshape(log2(totRich(:,2)./totRich(:,1)),20,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
hold on
refl=refline(0,0);
set(refl,'color','k')
xlabel 'movement rate'
ylabel 'log ratio richness'
subplot(3,2,4)
violin(reshape(log2(totBeta(:,2)./totBeta(:,1)),20,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
hold on
refl=refline(0,0);
set(refl,'color','k')
xlabel 'movement rate'
ylabel 'log ratio beta richness'
subplot(3,2,5)
violin(reshape(log2(totTrophicLevel(:,2)./totTrophicLevel(:,1)),20,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
hold on
refl=refline(0,0);
set(refl,'color','k')
xlabel 'movement rate'
ylabel 'log ratio trophic level'
subplot(3,2,6)
violin(reshape(log2(maxTrophicLevel(:,2)./maxTrophicLevel(:,1)),20,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
hold on
refl=refline(0,0);
set(refl,'color','k')
xlabel 'movement rate'
ylabel 'log ratio max trophic level'
