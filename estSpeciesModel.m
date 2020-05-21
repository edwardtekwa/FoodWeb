function [r,a,K,flag,raR2,r_T,betas,K_T_ratio,r_T_ratio]=estSpeciesModel(Btrans,dBtrans,gainBtrans,P) %r, a dimensions: species, patch (temperature)
%Edward Tekwa Nov 17, 17
%estimate per-mass intrinsic growth rate and competition terms
global zs
options=optimoptions(@lsqlin,'OptimalityTolerance',1e-12,'StepTolerance',1e-12,'ConstraintTolerance',1e-12,'display','off');

r=NaN(size(Btrans,1),size(Btrans,2));
a=NaN(2,size(Btrans,2)); %a=NaN(size(Btrans,1),size(Btrans,2)); %first row: self competition, second row: interaction with all other species in patch
raR2=NaN(1,size(Btrans,2)); %model R^2 values for each species
numraPts=NaN(size(Btrans,1),size(Btrans,2));
K=NaN(size(Btrans,1),size(Btrans,2)); %maximum carrying capacity
flag=NaN(1,size(Btrans,2));
Bmax=NaN(size(Btrans,1),size(Btrans,2)); %biomass
Bmean=NaN(size(Btrans,1),size(Btrans,2));
dBmax=NaN(size(dBtrans,1),size(dBtrans,2)); %change
dBmean=NaN(size(dBtrans,1),size(dBtrans,2));
gainBmax=NaN(size(dBtrans,1),size(dBtrans,2)); %productivity
gainBmean=NaN(size(dBtrans,1),size(dBtrans,2));
r_T=NaN(size(Btrans,1),size(Btrans,2));
K_T=NaN(size(Btrans,1),size(Btrans,2));
betas=NaN(1,size(Btrans,2));

dBperM=dBtrans./Btrans;
for species=1:size(Btrans,2)
    %X3=NaN(size(Btrans,1)*size(Btrans,3),size(Btrans,1)+1);
    X3=NaN(size(Btrans,1)*size(Btrans,3),size(Btrans,1)+2);
    Y3=NaN(size(Btrans,1)*size(Btrans,3),1); %all dB data for each species
    nPatches=0;
    if sum(Btrans(:,species,end))>eps %check if species is non-extinct at the end of the transcient period
        for patch=1:size(Btrans,1)
            if Btrans(patch,species,end)>eps && max(gainBtrans(patch,species,:))>-Inf
                nPatches=nPatches+1;
                Y=reshape(dBperM(patch,species,:),[],1);
                X1=reshape(Btrans(patch,species,:),[],1); %self biomass within patch (for estimating a_ii)
                X0=reshape(sum(Btrans(patch,[1:species-1,species+1:end],:),2),[],1); %total other biomass within patch (for estimating a_ij)
                Bmax(patch,species)=max(X1);
                Bmean(patch,species)=mean(X1);
                dBmax(patch,species)=max(dBtrans(patch,species,:));
                dBmean(patch,species)=mean(dBtrans(patch,species,:));
                gainBmax(patch,species)=max(gainBtrans(patch,species,:));
                gainBmean(patch,species)=mean(gainBtrans(patch,species,:));
                %X2=[ones(size(X1)), X1]; %separate intercepts (r) and biomass-dependence (a) at each location for each species
                %separate intercepts (r) for each location and one biomass-dependence (a)
                X3((patch-1)*length(X1)+1:patch*length(X1),patch)=ones(size(X1));
                X3((patch-1)*length(X1)+1:patch*length(X1),end-1)=X1;
                X3((patch-1)*length(X1)+1:patch*length(X1),end)=X0; %add interaction with others
                Y3((patch-1)*length(X1)+1:patch*length(X1))=Y;
                %lm=fitlm(X1,Y, 'y ~ x1');
                %[lm2,resnorm,residual,exitflag,output,lambda]=lsqlin(X2,Y,[-1 -max(X1)],0,[],[],[-Inf;-Inf],[Inf;0],[],options); %r/a<max(Btrans) or x1/x2<=max(Btrans)
                %r(patch,species)=lm.Coefficients.Estimate(1);
                %a(patch,species)=lm.Coefficients.Estimate(2);
                %r(patch,species)=lm2(1);
                %a(patch,species)=lm2(2);
                %flag(patch,species)=exitflag;
                %raR2(patch,species)=lm.Rsquared.Ordinary;
                %numraPts(patch,species)=lm.NumObservations;
            end
        end
    end
    X4=X3(~isnan(X3(:,end)),:); %delete NaN rows
    X4=X4(:,nansum(X3)>0); %delete NaN columns;
    X4(isnan(X4))=0; %replace NaNs in the design matrix with zero;
    Y4=Y3(~isnan(X3(:,end)),:); %delete NaN rows
    %     %construct lower and upper constraints
    %     lower=-Inf*ones(size(X4,2),1);
    %     upper=[Inf*ones(size(X4,2)-1,1);0];
    %construct Bmax constraints
    %nPatches=size(X4,2)-1;
    
    if ~isempty(Y4)
        %IneqA=[-diag(ones(nPatches,1)) -Bmax((nansum(X3(:,1:end-2))>0)',species) zeros(nPatches,1); diag(ones(nPatches,1)).*Bmax((nansum(X3(:,1:end-2))>0)',species) zeros(nPatches,2); [zeros(1,nPatches) 1 0]]; %constrains by rows: r/-a<=Bmax for each patch, r*Bmax<=gainBmax, a<=0
        %Ineqb=[zeros(nPatches,1);gainBmax((nansum(X3(:,1:end-2))>0)',species);0]; %last number is how large intraspecific competition must be (-a_ii)
        IneqA=[diag(ones(nPatches,1)).*Bmax((nansum(X3(:,1:end-2))>0)',species) zeros(nPatches,2); [zeros(1,nPatches) 1 0]]; %constrains by rows: r/-a<=Bmax for each patch, r*Bmax<=gainBmax, a<=0
        Ineqb=[gainBmax((nansum(X3(:,1:end-2))>0)',species);0]; %last number is how large intraspecific competition must be (-a_ii)
        [lm4,resnorm,residual,exitflag,output,lambda]=lsqlin(X4,Y4,IneqA,Ineqb,[],[],[],[],[],options); %r/a<max(Btrans) or x1/x2<=max(Btrans) %r/a<max(Btrans) or x1/x2<=max(Btrans)
        %[lm4,resnorm,residual,exitflag,output,lambda]=lsqlin(X4,Y4,[],[],[],[],[],[],[],options); %r/a<max(Btrans) or x1/x2<=max(Btrans) %r/a<max(Btrans) or x1/x2<=max(Btrans)
        EqA=[zeros(1,nPatches) 0 1]; %set a_ij
        EqB=0; %set a_ij=0
        [lm4LG,resnormLG,residualLG,exitflagLG,outputLG,lambdaLG]=lsqlin(X4,Y4,IneqA,Ineqb,EqA,EqB,[],[],[],options); %no interspecific interaction (LG=logistic growth model)
        r(nansum(X3(:,1:end-2))>0,species)=lm4(1:end-2);
        a(:,species)=lm4(end-1:end);
        raR2(species)=1-resnorm/(var(Y4)*(length(Y4)-1));
        flag(1,species)=exitflag;
        %estimate left-skewed normal thermal envelope of species s:
        zs=P.z(species); %get current species optimal temperature
        beta=nlinfit(P.T(nansum(X3(:,1:end-2))>0),lm4(1:end-2),@skewThEnv,1e-5);
        r_T(:,species)=skewThEnv(beta,P.T);
        betas(1,species)=beta;
    end
end

K=r./-a(1,:);
K_T=r_T./-a(1,:); %carrying capacities based on temperature-dependent r and
sumBmean=nansum(Bmean(:));
sumBmax=nansum(Bmax(:));
sumK=nansum(K(:));
sumK_T=nansum(K_T(:));
K_T_ratio=sumK_T/sumBmax; %predicted maximum total biomass over observed maximum total biomass
r_T_ratio=nansum(r_T(:))/nansum(gainBmax(:));%predicted maximum total productivity over observed maximum total productivity

% scrsz = get(0,'ScreenSize');
% figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/5 scrsz(4)]);
% subplot(3,1,1)
% bar(max(r));
% xlabel 'species (ordered by increasing body size)'
% ylabel 'maximum intrinsic growth rate (r)'
% xlim([0,100])
% title([{'Lotka-Volterra parameter estimates from'}; {'food-web simulation'}])
% subplot(3,1,2)
% bar((-a(1,:)).^.1);
% xlabel 'species (ordered by increasing body size)'
% ylabel 'self competition (-a_i_i)^{1/10}'
% xlim([0,100])
% subplot(3,1,3)
% bar(a(2,:));
% xlabel 'species (ordered by increasing body size)'
% ylabel 'competition with others (a_i_o)'
% xlim([0,100])

%temp-independent a