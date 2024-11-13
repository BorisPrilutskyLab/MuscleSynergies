clear variables;
DataFolder = '..\DataHindTied\';

N_syn = 5;
NShuffles = 50;
CondName = 'Tables_AllCondAllSpeeds';
StatTablesFile = strcat(DataFolder,CondName,'_syn-',num2str(N_syn),'_sh',num2str(NShuffles),'.mat');
load(StatTablesFile,'TH','TW','TH_WoAll','TW_WoAll')

CondNameR = 'AllCondAllSpeedsR'; %renamed
NameStr = strcat(CondNameR,'_syn',num2str(N_syn),'_sh',num2str(NShuffles),'.txt');

% testing
TH0=TH;TW0=TW;TH_WoAll0=TH_WoAll;TW_WoAll0=TW_WoAll;
%rename muscles
Muscles = unique(TW.Muscle);
% Plots order
MusclesNew ={'FDL','PLO','SOL','LG','MG','PLA','VL','BFA'...
    ,'GLU','CF','BFP','ST','TA','SRT','IP'}; %,'FL'
NumbersAdd = [1,2,3,4,5,6,7,8,91,92,93,94,95,96,97];

TW1 = TW0;
TW_WoAll1=TW_WoAll0;
for iMus = 1:length(MusclesNew)
    Muscle = MusclesNew{iMus};
    MuscleNewName = {strcat(num2str(NumbersAdd(iMus)),Muscle)};
    TW_Mus = TW0(ismember(TW0.Muscle,Muscle),:);
    TW1(ismember(TW1.Muscle,Muscle),:).Muscle = repmat(MuscleNewName,height(TW_Mus),1);
    TW_WoAll_Mus = TW_WoAll0(ismember(TW_WoAll0.Muscle,Muscle),:);
    TW_WoAll1(ismember(TW_WoAll0.Muscle,Muscle),:).Muscle = repmat(MuscleNewName,height(TW_WoAll_Mus),1);
end


TH = renamevars(TH0,{'SideSpeed','ContraSideSpeed'},{'IpsiSpeed','ContraSpeed'});
TW = renamevars(TW1,{'SideSpeed','ContraSideSpeed'},{'IpsiSpeed','ContraSpeed'});
TH_WoAll = renamevars(TH_WoAll0,{'SideSpeed','ContraSideSpeed'},{'IpsiSpeed','ContraSpeed'});
TW_WoAll = renamevars(TW_WoAll1,{'SideSpeed','ContraSideSpeed'},{'IpsiSpeed','ContraSpeed'});

writetable(TH,strcat('../Data/TH_',NameStr))
writetable(TW,strcat('../Data/TW_',NameStr))

writetable(TH_WoAll,strcat('../Data/TH_WoAll_',NameStr))
writetable(TW_WoAll,strcat('../Data/TW_WoAll_',NameStr))

