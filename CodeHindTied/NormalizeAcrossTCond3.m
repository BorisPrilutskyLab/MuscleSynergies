clear variables
DataFolder = '..\DataHindTied\';

Conditions = {'INTACT','SPINAL'};
TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';%'SPLIT';%

TCleanTCond = [];

for iCond = 1:length(Conditions)
    Condition = Conditions{iCond};

    TableSpeedsFilename = strcat(DataFolder,Condition,'_',TreadmillCondition,'_','TablesVicB3.mat');
    load(TableSpeedsFilename,'TCleanSpeeds')

    TClean = TCleanSpeeds(~TCleanSpeeds.Outlier,:);
    TClean(ismember(TClean.Muscle,{'VLL','GR','VM','FHL','VLM','SRL','RTL','RTM','SM'}),:) = [];
    TClean.Outlier = repmat({Condition},height(TClean),1);
    TClean = renamevars(TClean,'Outlier','Condition');
    TCleanTCond = [TCleanTCond;TClean];
end

TCleanTCond2 = addvars(TCleanTCond,TCleanTCond.EMG_Norm,'NewVariableNames','EMG_Norm2'); 

Muscles = unique(TCleanTCond.Muscle);
for iMus = 1:length(Muscles)
    Muscle = Muscles{iMus};
    TMus = TCleanTCond(ismember(TCleanTCond.Muscle,Muscle),:);
    Cats = unique(TMus.Cat);
    for iCat = 1:length(Cats)
        Cat = Cats{iCat};
        TCat = TMus(ismember(TMus.Cat,Cat),:);
        Sides = unique(TCat.MSide);
        for iSide =1:length(Sides)
            Side = Sides{iSide};
            TSide = TCat(ismember(TCat.MSide,Side),:);
            Conns = unique(TSide.Conn);
            for iConn = 1:length(Conns)
                Conn = Conns{iConn};
                TConn = TSide(ismember(TSide.Conn,Conn),:);
                EMGmax = max(cellfun(@max, TConn.EMG));
                EMGNorm = cellfun(@(x) x/EMGmax, TConn.EMG,UniformOutput=false);
                TCleanTCond2(ismember(TCleanTCond2.Muscle,Muscle)&...
                    ismember(TCleanTCond2.Cat,Cat)&ismember(TCleanTCond2.MSide,Side)&...
                    ismember(TCleanTCond2.Conn,Conn),:).EMG_Norm2 = EMGNorm;
            end 
        end
    end
end

TableTCondFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB3.mat');
save(TableTCondFilename,'TCleanTCond2')
