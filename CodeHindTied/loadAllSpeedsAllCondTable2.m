function [TClean,CleanMuscles,NCleanMuscles,N_Allmuscles,SideSpeeds,...
    NSideSpeeds,MusMap,MusclesNew] = loadAllSpeedsAllCondTable2(DataFolder,TreadmillCondition)


TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB.mat');
load(TableSpeedsFilename,'TCleanTCond15R')

TClean = TCleanTCond15R;
CleanMuscles = unique(TClean.Muscle); %all muscles in table

NCleanMuscles = cellfun(@(x) sum(ismember(TClean.Muscle,x)), CleanMuscles);

N_Allmuscles = length(CleanMuscles);
SideSpeeds = unique(TClean.SideSpeed);%Speed of treadmill of side of limb
NSideSpeeds = cellfun(@(x) sum(ismember(TClean.SideSpeed,x)), SideSpeeds);

% Plots order
MusclesNew ={'FDL','PLO','SOL','LG','MG','PLA','VL','BFA'...
    ,'GLU','CF','BFP','ST','TA','SRT','IP'}; %,'FL'
% MusclesNew ={'FDL','FHL','PLO','SOL','LG','MG','PLA','VL','VM','BFA'...
%     ,'GLU','CF','BFP','ST','TA','SRT','FL','IP'};
% ,'SM','GR','VLL'
MusclesOriginal = CleanMuscles;%Original order
MusMap = zeros(1,N_Allmuscles);
for i = 1:N_Allmuscles
    ID = find(ismember(MusclesOriginal,MusclesNew(i)));
    if ~isempty(ID)
        MusMap(i) = ID;
    end
end