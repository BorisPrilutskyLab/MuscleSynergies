clear variables
DataFolder = '..\DataHindTied\';

NTime= 101;
Conditions = {'INTACT','SPINAL'};
TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';
TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB3.mat');
load(TableSpeedsFilename,'TCleanTCond16R')

% Norm2 across all speeds/conditions
TCleanTCond3 = addvars(TCleanTCond16R,TCleanTCond16R.FBurstAveNormMax,...
    TCleanTCond16R.FBurstAveNormPeakEnv, TCleanTCond16R.EBurstAveNormMax,...
    TCleanTCond16R.EBurstAveNormPeakEnv, TCleanTCond16R.EBurstAveNormMax,...
    TCleanTCond16R.EBurstAveNormPeakEnv,...
    TCleanTCond16R.EBurstAvePrePC,TCleanTCond16R.EBurstAvePostPC,...
     'NewVariableNames',{'FBurstAveNormMax2',...
    'FBurstAveNormPeakEnv2','EBurstAveNormMax2','EBurstAveNormPeakEnv2',...
    'AllBurstAveNormMax2','AllBurstAveNormPeakEnv2',...
    'EBurstAvePrePCNormMax2','EBurstAvePostPCNormMax2'});

Muscles = unique(TCleanTCond16R.Muscle);
for iMus = 1:length(Muscles)
    Muscle = Muscles{iMus};
    TMus = TCleanTCond16R(ismember(TCleanTCond16R.Muscle,Muscle),:);
    Cats = unique(TMus.Cat);
    for iCat = 1:length(Cats)
        Cat = Cats{iCat};
        TCat = TMus(ismember(TMus.Cat,Cat),:);
        Sides = unique(TCat.MSide);
        for iSide = 1:length(Sides)
            Side = Sides{iSide};
            TSide = TCat(ismember(TCat.MSide,Side),:);
            Conns = unique(TSide.Conn);
            for iConn = 1:length(Conns)
                Conn = Conns{iConn};
                TConn = TSide(ismember(TSide.Conn,Conn),:);

                maxEnv = max(cellfun(@max, TConn.EMG));
                if ~isempty(maxEnv)
                    FBurstAveEmptyID = cellfun(@isempty, TConn.FBurstAve);
                    if ~all(FBurstAveEmptyID)
                        FBurstAvemaxC = cellfun(@max, TConn.FBurstAve(~FBurstAveEmptyID),UniformOutput=false);
                        FBurstAvemax = max(cellfun(@max,FBurstAvemaxC));
                        FBurstAveNorm2Max = cellfun(@(x) x/FBurstAvemax, TConn.FBurstAve,UniformOutput=false);

                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).FBurstAveNormMax2 = FBurstAveNorm2Max;

                        FBurstAveNorm2Env = cellfun(@(x) x/maxEnv, TConn.FBurstAve,UniformOutput=false);

                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).FBurstAveNormPeakEnv2 = FBurstAveNorm2Env;
                    end

                    EBurstAveEmptyID = cellfun(@isempty, TConn.EBurstAve);
                    if ~all(EBurstAveEmptyID)

                        EBurstAvemaxC = cellfun(@max, TConn.EBurstAve(~EBurstAveEmptyID),UniformOutput=false);
                        EBurstAvemax = max(cellfun(@max,EBurstAvemaxC));
                        
                        EBurstAveNorm2Max = cellfun(@(x) x/EBurstAvemax, TConn.EBurstAve,UniformOutput=false);
                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).EBurstAveNormMax2 = EBurstAveNorm2Max;

                        EBurstAvePrePCNormMax2 = cellfun(@(x) x/EBurstAvemax, TConn.EBurstAvePrePC,UniformOutput=false);
                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).EBurstAvePrePCNormMax2 = EBurstAvePrePCNormMax2;
                        
                        EBurstAvePostPCNormMax2 = cellfun(@(x) x/EBurstAvemax, TConn.EBurstAvePostPC,UniformOutput=false);
                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).EBurstAvePostPCNormMax2 = EBurstAvePostPCNormMax2;

                        EBurstAveNorm2Env = cellfun(@(x) x/maxEnv, TConn.EBurstAve,UniformOutput=false);
                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).EBurstAveNormPeakEnv2 = EBurstAveNorm2Env;
                    end
                    %All
                    AllBurstAveEmptyID = cellfun(@isempty, TConn.All_BurstAve);
                    if ~all(AllBurstAveEmptyID)

                        AllBurstAvemaxC = cellfun(@max, TConn.All_BurstAve(~AllBurstAveEmptyID),UniformOutput=false);
                        AllBurstAvemax = max(cellfun(@max,AllBurstAvemaxC));
                        AllBurstAveNorm2Max = cellfun(@(x) x/AllBurstAvemax, TConn.All_BurstAve,UniformOutput=false);
                                                
                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).AllBurstAveNormMax2 = AllBurstAveNorm2Max;

                        AllBurstAveNorm2Env = cellfun(@(x) x/maxEnv, TConn.All_BurstAve,UniformOutput=false);

                        TCleanTCond3(ismember(TCleanTCond3.Muscle,Muscle)&...
                            ismember(TCleanTCond3.Cat,Cat)&ismember(TCleanTCond3.MSide,Side)&...
                            ismember(TCleanTCond3.Conn,Conn),:).AllBurstAveNormPeakEnv2 = AllBurstAveNorm2Env;
                    end

                end
            end
        end
    end
end

save(TableSpeedsFilename,'TCleanTCond3','-append');

TBurstVars = {'Cat','Muscle','Condition','IpsiSpeed','ContraSpeed','Function',...
    'CycleTime_sec','StanceTime_sec','SwingTime_sec','SwingTimeNorm','DF','BurstOnNorm','BurstOffNorm',...
    'BurstAveNormMax2','BurstAveNormPeakEnv2',...
    'BurstAvePrePCNormMax2','BurstAvePostPCNormMax2',...
    'BurstPrePCtime_sec','BurstPostPCtime_sec',...
    'BurstPrePCtimeNorm','BurstPostPCtimeNorm'};
%};
%    'All_BurstAvePrePCNormMax2','All_BurstAvePostPCNormMax2',...

TBurst = cell2table(cell(0,length(TBurstVars)),'VariableNames',TBurstVars);
for iSample = 1:height(TCleanTCond3)
    TSample = TCleanTCond3(iSample,:);

    if ~isempty(TSample.FBurstOn{1})
        for iB = 1:length(TSample.FBurstOn{1})
            Cat = TSample.Cat;
            Muscle = TSample.Muscle;
            Condition = TSample.Condition;
            IpsiSpeed = TSample.SideSpeed;
            ContraSpeed = TSample.ContraSideSpeed;
            CycleTime = TSample.CycleTime;
            StanceTime = TSample.StanceTime;
            DF = TSample.DF;
            BurstOnNorm = TSample.FBurstOnNorm{1}(iB);
            BurstOffNorm = TSample.FBurstOffNorm{1}(iB);
            BurstAveNormMax2 = TSample.FBurstAveNormMax2{1}(iB);
            BurstAveNormPeakEnv2 = TSample.FBurstAveNormPeakEnv2{1}(iB);
            Function = {'Flex'};

            if isempty(TSample.EBurstAvePrePCNormMax2{1})
                BurstAvePrePCNormMax2 = nan;
                BurstAvePostPCNormMax2 = nan;
                BurstPrePCtime_sec = nan;
                BurstPostPCtime_sec = nan;
                BurstPrePCtimeNorm = nan;
                BurstPostPCtimeNorm = nan;

            else
                BurstAvePrePCNormMax2 = TSample.EBurstAvePrePCNormMax2{1}(iB);
                BurstAvePostPCNormMax2 = TSample.EBurstAvePostPCNormMax2{1}(iB);
                BurstPrePCtime_sec = TSample.EBurstPrePCtime_sec{1}(iB);
                BurstPostPCtime_sec = TSample.EBurstPostPCtime_sec{1}(iB);
                BurstPrePCtimeNorm = TSample.EBurstPrePCtimeNorm{1}(iB);
                BurstPostPCtimeNorm = TSample.EBurstPostPCtimeNorm{1}(iB);
            end
            PCNorm = TSample.PCNorm;
            PC_sec = TSample.PC_sec;

            TF = cell2table([Cat,Muscle,Condition,IpsiSpeed,ContraSpeed,...
                Function,CycleTime,StanceTime,PC_sec,PCNorm,DF,BurstOnNorm,BurstOffNorm,...
                BurstAveNormMax2,BurstAveNormPeakEnv2,...
                BurstAvePrePCNormMax2,BurstAvePostPCNormMax2,...
                BurstPrePCtime_sec,BurstPostPCtime_sec,...
                BurstPrePCtimeNorm,BurstPostPCtimeNorm,...
                ],'VariableNames',TBurstVars);
            TBurst = [TBurst;TF];
        end
    end
    if ~isempty(TSample.EBurstOn{1})
        for iB = 1:length(TSample.EBurstOn{1})
            Cat = TSample.Cat;
            Muscle = TSample.Muscle;
            Condition = TSample.Condition;
            IpsiSpeed = TSample.SideSpeed;
            ContraSpeed = TSample.ContraSideSpeed;
            CycleTime = TSample.CycleTime;
            StanceTime = TSample.StanceTime;
            DF = TSample.DF;
            BurstOnNorm = TSample.EBurstOnNorm{1}(iB);
            BurstOffNorm = TSample.EBurstOffNorm{1}(iB);
            BurstAveNormMax2 = TSample.EBurstAveNormMax2{1}(iB);
            BurstAveNormPeakEnv2 = TSample.EBurstAveNormPeakEnv2{1}(iB);
            Function = {'Ext'};

            BurstAvePrePCNormMax2 = TSample.EBurstAvePrePCNormMax2{1}(iB);
            BurstAvePostPCNormMax2 = TSample.EBurstAvePostPCNormMax2{1}(iB);
            BurstPrePCtime_sec = TSample.EBurstPrePCtime_sec{1}(iB);
            BurstPostPCtime_sec = TSample.EBurstPostPCtime_sec{1}(iB);
            BurstPrePCtimeNorm = TSample.EBurstPrePCtimeNorm{1}(iB);
            BurstPostPCtimeNorm = TSample.EBurstPostPCtimeNorm{1}(iB);
            PCNorm = TSample.PCNorm;
            PC_sec = TSample.PC_sec;

            TE = cell2table([Cat,Muscle,Condition,IpsiSpeed,ContraSpeed,...
                Function,CycleTime,StanceTime,PC_sec,PCNorm,DF,BurstOnNorm,BurstOffNorm,...
                BurstAveNormMax2,BurstAveNormPeakEnv2,...
                BurstAvePrePCNormMax2,BurstAvePostPCNormMax2,...
                BurstPrePCtime_sec,BurstPostPCtime_sec,...
                BurstPrePCtimeNorm,BurstPostPCtimeNorm,...
                ],'VariableNames',TBurstVars);
            TBurst = [TBurst;TE];
        end
    end
end

TBurstsFilename = strcat(DataFolder,TreadmillCondition,'_','TBursts3.txt');
writetable(TBurst,TBurstsFilename)


