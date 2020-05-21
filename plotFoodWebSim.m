%plotFoodWebSim
%Edward Tekwa May 15, 18
%Diagnostics: plot food web simulation outcomes versus model projections
%Instruction: load simulation .mat file first

figure('Color', [1 1 1]);

subplot(2,2,1);
bar(sum(mean(Btrans,3)),'k');
xlabel 'species'
ylabel 'transcient mean biomass'
ylim([0 15])
yyaxis right
refline(0,sum(sum(Btrans(:,:,end))>11*eps))
ylabel 'richness'
ylim([0 sum(sum(Btrans(:,:,end))>11*eps)+1])

subplot(2,2,2);
bar(sum(K_T1),'b');
xlabel 'species'
ylabel 'single-species transcient mean biomass'
ylim([0 15])
yyaxis right
refline(0,sum(sum(K_T1)>11*eps))
ylabel 'richness'
ylim([0 sum(sum(Btrans(:,:,end))>11*eps)+1])

subplot(2,2,3);
bar(sum(mean(B_yrs,3)),'r');
xlabel 'species'
ylabel 'warming mean biomass'
ylim([0 15])
yyaxis right
refline(0,sum(sum(B_yrs(:,:,end))>11*eps))
ylabel 'richness'
ylim([0 sum(sum(Btrans(:,:,end))>11*eps)+1])

subplot(2,2,4);
bar(sum(mean(BLV1_yrs,3)),'m');
xlabel 'species'
ylabel 'single-species warming mean biomass'
ylim([0 15])
yyaxis right
refline(0,sum(sum(BLV1_yrs(:,:,end))>11*eps))
ylabel 'richness'
ylim([0 sum(sum(Btrans(:,:,end))>11*eps)+1])

numYrs=20;
% plot_demog_spatial_Ball(B_yrs, P, gainB_yrs,0,numYrs);
% plot_demog_spatial_Ball(BLV1_yrs, P, gainBLV1_yrs,0,numYrs);
% plot_demog_spatial_Ball(Btrans, P, gainBtrans,0,numYrs);