function [doPlot,doPlotEach,doDispR2,pos,FontSz,FontSzLabel,CapSz,...
    FontWeightLabel,NConditions,CindCond,...
    CindCondInBetween,CtCond,Cx2] = ...
    SetPlotParametersAllCondAllSPeed2(TreadmillCondition,NTime)

doPlot = true;%false;%
doPlotEach =  1;%false;%true;%
doDispR2 = false;% true;

%Plot parameters
pos = 1000*[0.01,    0,    1,    0.9];% pos = 1000*[0.01,    0,    2,    1];
FontSz = 6;% FontSz = 10;
FontSzLabel =12;
CapSz = 3;
FontWeightLabel ='bold';



NConditions = 6;
% if strcmp(TreadmillCondition,'TIED')
%     NConditions = 3;
% else
%     NConditions = 6;
% end
%NConditionsPlus = NConditions + 1;

%H time frames for condition/combination
CindCond(1:NConditions) = {1:NTime};
CindCond(NConditions+1) = {1:NConditions*NTime};%comCondComSpeed
CindCond(NConditions+2) = {1:3*NTime};%IntactComSpeed
CindCond(NConditions+3) = {1:3*NTime};%SpinalComSpeed
CindCond(NConditions+4) = {1:2*NTime};%ConCondSpeed0.4
CindCond(NConditions+5) = {1:2*NTime};%ConCondSpeed0.7
CindCond(NConditions+6) = {1:2*NTime};%ConCondSpeed1.0

CindCondInBetween = cell(1,NConditions*2);
CtCond = cell(1,NConditions*2);
Cx2 = cell(1,NConditions*2);
for iCond = 1:NConditions
    CindCondInBetween(iCond) = {[1:NTime,1+(2*NConditions-1)*NTime:2*NConditions*NTime]};
    CtCond{iCond} = 1+NTime*(iCond-1):iCond*NTime;
end
CindCondInBetween(NConditions+1) = {1:2*NConditions*NTime};%comCondComSpeed
CindCondInBetween(NConditions+2) = {[1:3*NTime,1+(2*NConditions-3)*NTime:2*NConditions*NTime]};%IntactComSpeed
CindCondInBetween(NConditions+3) = {[1:3*NTime,1+(2*NConditions-3)*NTime:2*NConditions*NTime]};%SpinalComSpeed
CindCondInBetween(NConditions+4) = {[1:2*NTime,1+(2*NConditions-2)*NTime:2*NConditions*NTime]};%ConCondSpeed0.4
CindCondInBetween(NConditions+5) = {[1:2*NTime,1+(2*NConditions-2)*NTime:2*NConditions*NTime]};%ConCondSpeed0.7
CindCondInBetween(NConditions+6) = {[1:2*NTime,1+(2*NConditions-2)*NTime:2*NConditions*NTime]};%ConCondSpeed1.0

CtCond{NConditions+1} = 1:NConditions*NTime;
CtCond{NConditions+2} = 1:3*NTime;
CtCond{NConditions+3} = 1:3*NTime;
CtCond{NConditions+4} = 1:2*NTime;
CtCond{NConditions+5} = 1:2*NTime;
CtCond{NConditions+6} = 1:2*NTime;

% for NConditions = 2;
%CindCondInBetween(1:NConditions) = {[1:NTime,1+(2*NConditions-1)*NTime:2*NConditions*NTime]};

% set for NConditions = 3, will change for split
%CtCond = {1:NTime,1+NTime:2*NTime,1+2*NTime:3*NTime,1:3*NTime};%time for plot H conditions

for iCond = 1:NConditions+6
    Cx2{iCond} = [CtCond{iCond},fliplr(CtCond{iCond})];
end
