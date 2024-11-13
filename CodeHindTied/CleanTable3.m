clear variables
DataFolder = '..\DataHindTied\';

Conditions = {'INTACT','SPINAL'};
TreadmillConditions = {'SPLIT','TIED'};

Condition = 'INTACT';%'SPINAL';%
TreadmillCondition = 'TIED';%'SPLIT';%
switch TreadmillCondition
    case 'SPLIT'
        SpeedSets = {'0.4-0.5','0.4-0.7','0.4-1.0','0.5-0.4','0.7-0.4','1.0-0.4'};%
    case 'TIED'
        SpeedSets = {'0.4','0.7','1.0'};
end

TCleanSpeeds = [];

for iSpeedSet = 1:length(SpeedSets)
    SpeedSet = SpeedSets{iSpeedSet};

    TableFilename = strcat(DataFolder,Condition,'_',TreadmillCondition,'_',SpeedSet,'_','TablesVicB3.mat');
    load(TableFilename,'TWithOutliers');

    ListOfOuliersFilename = strcat(DataFolder,Condition,'_',TreadmillCondition,'_',SpeedSet,'_','OutliersList.mat');
    load(ListOfOuliersFilename,'TOutList');
    TOutList = unique(TOutList);

    TAddVars = {'Cat','MSide','Muscle','Conn','SideSpeed','BeltsSpeeds','CycleN'};
    TWithOutliersRedused = TWithOutliers(:,TAddVars);
    [C,ia,ib] = intersect(TWithOutliersRedused,TOutList);

    TClean = TWithOutliers;
    TClean(ia,:).Outlier = ones(length(ia),1);

    TCleanSpeeds = [TCleanSpeeds;TClean];

    save(TableFilename,'TClean','-append')
end

TableSpeedsFilename = strcat(DataFolder,Condition,'_',TreadmillCondition,'_','TablesVicB3.mat');

save(TableSpeedsFilename,'TCleanSpeeds')