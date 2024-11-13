clear variables
DataFolder = '..\DataHindTied\';

TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';

TBurstsFilename = strcat(DataFolder,TreadmillCondition,'_','TBursts3.txt');
TBurst = readtable(TBurstsFilename);

OutBurst = {'IP.Ext','TA.Ext','BFA.Flex','MG.Flex','VL.Flex','SOL.Flex'};
TBurstClean = TBurst;
OutMusFunc = squeeze(split(OutBurst,'.'));
for iOut = 1:length(OutMusFunc)
    TBurstClean(ismember(TBurstClean.Muscle,OutMusFunc(iOut,1))&ismember(TBurstClean.Function,OutMusFunc(iOut,2)),:)= [];
end

TBurstsFilename = strcat(DataFolder,TreadmillCondition,'_','TBurstClean3.txt');
writetable(TBurstClean,TBurstsFilename)
