%shiftMap.m
%Edward Tekwa June 28, 2019
%heat map of centroid shifts

%load directory of a simulation set, then extract using:
%nanmean(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[])) %extract specialist centroid shift means
%nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])) %extract generalist centroid shift means

BaseSpecialistShifts=[-0.0208165884265031,-0.855171938213938,-1.18358275616791,-1.61088019961662,-1.38446529263982,-0.206001596088462];
BaseGeneralistShifts=[-0.00446395823358861,-0.725794623695584,-0.891381778706218,-1.20019370160724,-1.12140611431660,-0.106229182792871];
PPMR208SpecialistShifts=[8.88178419700125e-17,-0.902689428802254,-1.27666877586124,-1.57973867331982,-1.56431715375861,-0.366226318550947];
PPMR208GeneralistShifts=[-0.499999820690379,-0.787154748453894,-0.951500081163275,-1.15780023514091,-1.46195278464926,-0.225698924233901];
GA1SpecialistShifts=[-0.144017926106498,-0.751070945925928,-1.13073668906497,-1.68411281962896,-1.18858632938823,-0.182939953823094];
GA1GeneralistShifts=[-0.0792679823951551,-0.754127338196848,-0.804235729615487,-1.20861021714902,-1.10018784959955,-0.0753113108806235];
EA069SpecialistShifts=[-0.179603211333503,-0.255961868265714,-0.318160311054744,-0.473793117461974,-0.657244953936487,-0.00270428889219053];
EA069GeneralistShifts=[-0.266986792008677,-0.335326003126187,-0.366617642197923,-0.576793009690616,-0.654372491989711,-0.00534194786082062];
Lambda02SpecialistShifts=[-0.0876189001550679,-1.19044356427958,-1.30043421279162,-1.85869585689099,-1.76831391312072,-0.267846700854546];
Lambda02GeneralistShifts=[-0.109539490061866,-0.801256558393785,-0.918307399477871,-1.24577872551754,-1.41731733784171,-0.156386779059745];

%descending order of specialist shift: Lambda, GA1, Base, PPMR, EA
ShiftTable=-(100/3)*[Lambda02GeneralistShifts;Lambda02SpecialistShifts;GA1GeneralistShifts;GA1SpecialistShifts;BaseGeneralistShifts;BaseSpecialistShifts;PPMR208GeneralistShifts;PPMR208SpecialistShifts;EA069GeneralistShifts;EA069SpecialistShifts];
yTreatments={'Generalist: low consumption efficiency';'Specialist: low consumption efficiency';'Generalist: high search rate';'Specialist: high search rate';'Generalist: base';'Specialist: base';'Generalist: small PPMR';'Specialist: small PPMR';'Generalist: low metabolism';'Specialist: low metabolism'};
%yTreatments={'G: \lambda-';'S: \lambda-';'G: v+';'S: v+';'G: base';'S: base';'G: PPMR-';'S: PPMR-';'G: D-';'S: D-'};
yOrder=[1 3 5 7 9 2 4 6 8 10];
%dispersalLabels={'0' '1' '10$$^3$$' '10$$^6$$' '10$$^9$$' '10$$^{12}$$'};
dispersalLabels={'0' '1' '10^3' '10^6' '10^9' '10^{12}'};


set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',14)
set(0,'defaulttextinterpreter','latex');  
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
scrsz = get(0,'ScreenSize');

shiftfig2=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.8 scrsz(4)/3]); %3.5
set(shiftfig2,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
h=heatmap(dispersalLabels,yTreatments(yOrder),ShiftTable(yOrder,:));
h.Colormap = flipud(jet);
h.FontSize = 14;
h.CellLabelColor = 'none';
h.GridVisible = 'off';
