%Add_warming_endstats_movement_script.m
%Edward Tekwa
%script to run analyses and save as .mat files

for pa=1:1
    Paths={'Lambda 02';'GA1'; 'EA 069'; 'PPMR 208'};
    pEdibleCases={'specialist';'generalist'};
    TempScenario=1; %position is without 0 change position [+2 +4 +6] (3 or 1)
    Path=Paths{pa}; %'GA1'; %'EA 069'; %'PPMR 208'; %''; %'Specialist food web'; %'Generalist food web';
    %first, use these lines for specialist food webs:
    TimeData=[Path ' ' pEdibleCases{1}]; %name of file
    Cases={'basalSize0.01_meanD-Inf_pInedible0.5';'basalSize0.01_meanD0_pInedible0.5';'basalSize0.01_meanD3_pInedible0.5';'basalSize0.01_meanD6_pInedible0.5';'basalSize0.01_meanD9_pInedible0.5';'basalSize0.01_meanD12_pInedible0.5'};
    AntiCases={'basalSize0.01_meanD-Inf_pInedible1';'basalSize0.01_meanD0_pInedible1';'basalSize0.01_meanD3_pInedible1';'basalSize0.01_meanD6_pInedible1';'basalSize0.01_meanD9_pInedible1';'basalSize0.01_meanD12_pInedible1'};
    Add_warming_endstats_movement_fixedDir_new
    
    %second, use these lines for generalist food webs:
    TimeData=[Path ' ' pEdibleCases{2}]; %name of file
    %Cases={'basalSize0.01_meanD-Inf';'basalSize0.01_meanD0';'basalSize0.01_meanD3';'basalSize0.01_meanD6';'basalSize0.01_meanD9';'basalSize0.01_meanD12'};
    Cases={'basalSize0.01_meanD-Inf_pInedible0';'basalSize0.01_meanD0_pInedible0';'basalSize0.01_meanD3_pInedible0';'basalSize0.01_meanD6_pInedible0';'basalSize0.01_meanD9_pInedible0';'basalSize0.01_meanD12_pInedible0'};
    AntiCases={'basalSize0.01_meanD-Inf_pInedible0.5';'basalSize0.01_meanD0_pInedible0.5';'basalSize0.01_meanD3_pInedible0.5';'basalSize0.01_meanD6_pInedible0.5';'basalSize0.01_meanD9_pInedible0.5';'basalSize0.01_meanD12_pInedible0.5'};
    Add_warming_endstats_movement_fixedDir_new
    clear
end