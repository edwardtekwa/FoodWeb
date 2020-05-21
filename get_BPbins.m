function [Ybins,Y2bins] = get_BPbins(B, P, prodB,numYrs)
%Edward Tekwa 5/1/18
%Take in B and P for all species in one simulation and return biomass sums
%binned by body size

X = log10(P.s.mi);
%construct histogram bins:
%binWidth=(ceil(max(X))-floor(min(X)))/numBins;
binWidth=0.25; %unit in log10 scale
numBins=(ceil(max(X))-floor(min(X)))/binWidth; %number of bins for heterotrophs (total number of bins is 1+numBins) 7
avgWindow=size(B,3)-numYrs+1:size(B,3); %average over last numYrs years
lengthWindow=length(avgWindow);

Ybins=zeros(1,numBins+1); %histogram bins for biomass
Y2bins=zeros(1,numBins+1); %histogram bins for productivity

for bin=2:numBins+1 %other bins are heterotrophs
    lowerX=floor(min(X))+(bin-2)*binWidth;
    upperX=floor(min(X))+(bin-1)*binWidth;
    %     Ybins(bin)=log10(sum(Yall(find(X>=lowerX & X<upperX))));
    %     Y2bins(bin)=log10(sum(Y2all(find(X>=lowerX & X<upperX))));
%     sortedTimeBiomass=sort(log10(mean(sum(B(:,find(X>=lowerX & X<upperX),avgWindow),2))));
%     sortedTimeProd=sort(log10(mean(sum(prodB(:,find(X>=lowerX & X<upperX),avgWindow),2))));
    sortedTimeBiomass=sort((mean(nansum(B(:,find(X>=lowerX & X<upperX),avgWindow),2))))./size(B,1);
    sortedTimeProd=sort((mean(nansum(prodB(:,find(X>=lowerX & X<upperX),avgWindow),2))))./size(B,1);
    Ybins(bin)=mean(sortedTimeBiomass);
    Y2bins(bin)=mean(sortedTimeProd);
end