function []=plot_demog_spatial_Bavg(B_all,BS_all,Prod_all,Z_all,Shift_all,P,tempChanges,numIt,moveRate)

%Edward Tekwa May 4, 2018
%plot spatial species distributions averaged over simulation replicates

maxRich=40; %6 or 40
maxBiomass=9; %1 or 8
maxShift=10;
minShift=-10;
scrsz = get(0,'ScreenSize');
figs(1)=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.25 scrsz(4)]);

BiomassColor=[0.7 0.7 1];
ProdColor=[1 0.7 0.7];
BiomassErrColor='b';
ProdErrColor='r';
CM=colormap(jet(128)); % set colormap
caxis([0 6]); % set colormap range
C = caxis;
left_color = [0 0 1];
right_color = [0 0 0];
set(figs(1),'defaultAxesColorOrder',[left_color; right_color]);

numPatch=size(B_all,1);
%scale matrices:
%B_all=B_all/numIt;
%Prod_all=Prod_all/numIt;

midCase=ceil(length(tempChanges)/2);
for Case=1:length(tempChanges)
    subplot(length(tempChanges),3,Case*3-2)
    %pick out data for current temperature 
    B_case=reshape(B_all(:,:,(Case-1)*numIt+1:Case*numIt),numPatch,[]);
    numNonZeroIt=sum(nansum(nansum(B_all(:,:,(Case-1)*numIt+1:Case*numIt),1),2)>eps)
    B_case=B_case/numNonZeroIt; %scale to non-extinct communities
    B_case(isnan(B_case))=0;
    Prod_case=reshape(Prod_all(:,:,(Case-1)*numIt+1:Case*numIt),numPatch,[])/numNonZeroIt;
    Prod_case(isnan(Prod_case))=0;
    BS_case=reshape(BS_all((Case-1)*numIt+1:Case*numIt,:)',1,[]);
    %sum sorted and scaled biomasses and productions of species in 128 body size
    %bins
    binnedB=zeros(numPatch,128);
    binnedProd=zeros(numPatch,128);
    for i = 1:length(B_case)
        Ci = fix((log10(BS_case(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
        binnedB(:,Ci)=binnedB(:,Ci)+B_case(:,i);
        binnedProd(:,Ci)=binnedProd(:,Ci)+Prod_case(:,i);
    end
    
    yyaxis left
    SpeciesDistr=area(binnedB,'LineStyle','none');
    hold on
    for i = 1:length(binnedB)
        SpeciesDistr(i).FaceColor=CM(i,:);
    end
    set(gca,'fontsize',15);
    box on;
    %xlabel(['Patch (Basal size=' num2str(P.s.m0) 'g, mean(D)=' num2str(HP.sdm(ei))  ', std(D)=' num2str(HP.sdv(ei)) ', std(optimal T)=' num2str(HP.tv(ei)) ')'],'fontsize',16)
    %xlabel(['Patch (Basal size=' num2str(P.s.m0) 'g'],'fontsize',16)
    if Case==1
        title(['connectivity=' moveRate ', n=' num2str(numIt)],'fontsize',16)
        %         at = linspace(1, 6, 6);
        %         b = colorbar('YTick',[0:1/5:1],'YTickLabel',at);
        %         set(get(b,'ylabel'), 'String', 'body size (log_{10}g', 'fontsize', 16)
    end
    if Case==2
        ylabel '^*no warming'
    end
    if Case==midCase
        ylabel('biomass (gm^{-3})','fontsize',16)
    end
    if Case==length(tempChanges)
        xlabel(['patch end temperature (^{\circ}C)'],'fontsize',16)
    end
    
    %plot basal biomass
    boundedline([1:numPatch], nanmean(Z_all(:,(Case-1)*numIt+1:Case*numIt),2)',[nanstd(Z_all(:,(Case-1)*numIt+1:Case*numIt)')*1.96./((sum(~isnan(Z_all(:,(Case-1)*numIt+1:Case*numIt)'))).^0.5)]','b','alpha');
    
    
    xlim([0.8 numPatch+0.2]);
    xticks([1:2:numPatch]);
    xticklabels(P.T(1:2:end)+tempChanges(Case));
    ylim([0 maxBiomass])
    
    yyaxis right
    B_space=B_all(:,:,(Case-1)*numIt+1:Case*numIt);
    B_space_richness=reshape(sum(B_space>eps,2),numPatch,[]);
    AlphaRichnessDistr=sort(B_space_richness,2); %sort each row separately from low to high richness per spatial location
    meanAlphaRichness=mean(AlphaRichnessDistr,2);
    %     loAlphaRichness=AlphaRichnessDistr(:,:,ceil(lengthWindow*0.025));
    %     hiAlphaRichness=AlphaRichnessDistr(:,:,ceil(lengthWindow*0.975));
    %     bar(meanAlphaRichness,'FaceColor','none','LineWidth',2,'EdgeColor',[0 0.45 0.74]);
    boundedline([1:numPatch], meanAlphaRichness',[nanstd(AlphaRichnessDistr')*1.96./((sum(~isnan(AlphaRichnessDistr'))).^0.5)]','k','alpha');
    %     hold on
    %     er=errorbar([1:numpatch],meanAlphaRichness,meanAlphaRichness-loAlphaRichness,hiAlphaRichness-meanAlphaRichness,'.');
    %     set(er,'LineWidth',2);
    if Case==midCase
        ylabel('local richness (\alpha)')
    end
    ylim([0 maxRich])
    
    
    subplot(length(tempChanges),3,Case*3-1) %plot histograms of # species at different body sizes
%    yyaxis left
    histogram(log10(BS_case(sum(B_case)>eps)),12,'FaceColor',[0.5 0.5 0.5]); %include only species that are not extinct
    if Case==midCase
        ylabel('frequency','fontsize',16)
    end
%     yyaxis right
%     ksdensity(log10(BS_case(sum(B_case)>eps)));
%     if Case==midCase
%         ylabel('probability density','fontsize',16)
%     end
    if Case==length(tempChanges)
        xlabel(['log_{10}(body size)'],'fontsize',16)
    end
    xlim([1 6])
    
    subplot(length(tempChanges),3,Case*3) %plot relationship between body size and biomass
    %BS_case=reshape(BS_all((Case-1)*numIt+1:Case*numIt,:)',1,[]);
    Shift_case=reshape(Shift_all(1,:,(Case-1)*numIt+1:Case*numIt),1,[]);
    BS_Shift=fitlm(log10(BS_case),Shift_case)
    numVals=100;
    xval = min(log10(BS_case)):(max(log10(BS_case))-min(log10(BS_case)))/(numVals-1):max(log10(BS_case));
    [ynew, ci] = predict(BS_Shift, xval');
    b = [ynew-ci(:,1) ci(:,2)-ynew];
    hold on
    boundedline(xval,ynew',b,'k')
    scatter(log10(BS_case),Shift_case,'k.')
    xlim([1 6])
    ylim([minShift maxShift])
    
    if Case==midCase
        ylabel('species shift','fontsize',16)
    end
    if Case==length(tempChanges)
        xlabel(['log_{10}(body size)'],'fontsize',16)
    end
end
