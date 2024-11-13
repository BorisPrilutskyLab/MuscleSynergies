clear variables
DataFolder = '..\DataHindTied\';

Min2BurstFlexOn = -.25;

Max2BurstExtOff = .7;
Min2BursExtOn = 0;

MaxFlexOn = 0.0; 
MinFlexOn = -.25;
MaxFlexOffAfterPC = .1;

MinExtOn = .15;
MaxExtOff = 1.0; %
MinExOffAfterPC = .1;

TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';
TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB3.mat');
load(TableSpeedsFilename,'TCleanTCond2')

TCleanTCond15 = TCleanTCond2(~ismember(TCleanTCond2.Muscle,'FL'),:);
TwoJT = TCleanTCond15(ismember(TCleanTCond15.Function,'NA'),:);
TCleanTCond15(ismember(TCleanTCond15.Function,'NA'),:).Function = repmat({'TwoJoint'},height(TwoJT),1);
TwoJointMuscles = unique(TCleanTCond15(ismember(TCleanTCond15.Function,'TwoJoint'),:).Muscle);
FlexExtMuscles = unique(TCleanTCond15(~ismember(TCleanTCond15.Function,'TwoJoint'),:).Muscle);


TCleanTCond16 = addvars(TCleanTCond15,...
    TCleanTCond15.All_BurstAvePrePC,TCleanTCond15.All_BurstAvePostPC,...
    TCleanTCond15.All_BurstAvePrePC,TCleanTCond15.All_BurstAvePrePC,...
    TCleanTCond15.All_BurstAvePrePC,TCleanTCond15.All_BurstAvePrePC,...
    TCleanTCond15.StanceTime,TCleanTCond15.StanceTime,...
    'NewVariableNames',{'EBurstAvePrePC','EBurstAvePostPC',...
    'EBurstPrePCtime_sec','EBurstPostPCtime_sec',...
    'EBurstPrePCtimeNorm','EBurstPostPCtimeNorm','PCNorm','PC_sec'});

TCleanTCond16R = TCleanTCond16;

for iMus = 1:length(TwoJointMuscles)
    Muscle = TwoJointMuscles{iMus};
    TMus = TCleanTCond16R(ismember(TCleanTCond16R.Muscle,Muscle),:);
    TMusR = TMus;
    for iSample = 1:height(TMus)
        TSample = TMus(iSample,:);
        SwTimeNorm = (TSample.CycleTime - TSample.StanceTime)/TSample.CycleTime;
        SwTime = (TSample.CycleTime - TSample.StanceTime);

        TMusR(iSample,:).PCNorm = SwTimeNorm;
        TMusR(iSample,:).PC_sec = SwTime;


        AllBurstOn = TSample.All_BurstOn;
        AllBurstOff = TSample.All_BurstOff;
        AllBurstOnNorm = TSample.AllBurstOnNorm;
        AllBurstOffNorm = TSample.AllBurstOffNorm;
        AllBurstAve = TSample.All_BurstAve;
        AllBurstAveNormMax = TSample.AllBurstAveNormMax;
        AllBurstAveNormPeakEnv = TSample.AllBurstAveNormPeakEnv;

        All_BurstAvePrePC = TSample.All_BurstAvePrePC;
        All_BurstAvePostPC = TSample.All_BurstAvePostPC;

        TMusR(iSample,:).FBurstOn = cellfun(@(x,y,z) z(x>Min2BurstFlexOn&x<SwTimeNorm*0.5&y<SwTimeNorm),AllBurstOnNorm,AllBurstOffNorm,AllBurstOn,'UniformOutput', false);
        TMusR(iSample,:).FBurstOff = cellfun(@(x,y,z) z(x>Min2BurstFlexOn&x<SwTimeNorm*0.5&y<SwTimeNorm),AllBurstOnNorm,AllBurstOffNorm,AllBurstOff,'UniformOutput', false);
        TMusR(iSample,:).FBurstOnNorm = cellfun(@(x,y,z) z(x>Min2BurstFlexOn&x<SwTimeNorm*0.5&y<SwTimeNorm),AllBurstOnNorm,AllBurstOffNorm,AllBurstOnNorm,'UniformOutput', false);
        TMusR(iSample,:).FBurstOffNorm = cellfun(@(x,y,z) z(x>Min2BurstFlexOn&x<SwTimeNorm*0.5&y<SwTimeNorm),AllBurstOnNorm,AllBurstOffNorm,AllBurstOffNorm,'UniformOutput', false);
        TMusR(iSample,:).FBurstAve = cellfun(@(x,y,z) z(x>Min2BurstFlexOn&x<SwTimeNorm*0.5&y<SwTimeNorm),AllBurstOnNorm,AllBurstOffNorm,AllBurstAve,'UniformOutput', false);
        TMusR(iSample,:).FBurstAveNormMax = cellfun(@(x,y,z) z(x>Min2BurstFlexOn&x<SwTimeNorm*0.5&y<SwTimeNorm),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormMax,'UniformOutput', false);
        TMusR(iSample,:).FBurstAveNormPeakEnv = cellfun(@(x,y,z) z(x>Min2BurstFlexOn&x<SwTimeNorm*0.5&y<SwTimeNorm),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormPeakEnv,'UniformOutput', false);

        TMusR(iSample,:).EBurstOn = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOn,'UniformOutput', false);
        TMusR(iSample,:).EBurstOff = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOff,'UniformOutput', false);
        TMusR(iSample,:).EBurstOnNorm = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOnNorm,'UniformOutput', false);
        TMusR(iSample,:).EBurstOffNorm = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOffNorm,'UniformOutput', false);
        TMusR(iSample,:).EBurstAve = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstAve,'UniformOutput', false);
        TMusR(iSample,:).EBurstAveNormMax = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormMax,'UniformOutput', false);
        TMusR(iSample,:).EBurstAveNormPeakEnv = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormPeakEnv,'UniformOutput', false);

        TMusR(iSample,:).EBurstAvePrePC = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,All_BurstAvePrePC,'UniformOutput', false);
        TMusR(iSample,:).EBurstAvePostPC = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,All_BurstAvePostPC,'UniformOutput', false);
 
        %remove nan and timing for nan
        EBurstAvePrePC_isnan = cellfun(@isnan, TMusR(iSample,:).EBurstAvePrePC,'UniformOutput', false);
        EBurstAvePrePC_length = cellfun(@length, EBurstAvePrePC_isnan);

        EBurstOnNorm = cellfun(@(x,y) x(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,'UniformOutput', false);
        EBurstOffNorm = cellfun(@(x,y) y(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,'UniformOutput', false);
         
        EBurstOn = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOn,'UniformOutput', false);
        EBurstOff = cellfun(@(x,y,z) z(x>Min2BursExtOn&y>SwTimeNorm&y<Max2BurstExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOff,'UniformOutput', false);

        TMusR(iSample,:).EBurstPrePCtimeNorm =  cellfun(@(x,y) x - y,num2cell(SwTimeNorm),EBurstOnNorm,'UniformOutput', false);
        TMusR(iSample,:).EBurstPostPCtimeNorm =  cellfun(@(x,y) - x + y,num2cell(SwTimeNorm),EBurstOffNorm,'UniformOutput', false);
        
        TMusR(iSample,:).EBurstPrePCtime_sec = cellfun(@(x,y) x - y,num2cell(SwTime),EBurstOn,'UniformOutput', false);
        TMusR(iSample,:).EBurstPostPCtime_sec =  cellfun(@(x,y) - x + y,num2cell(SwTime),EBurstOff,'UniformOutput', false);
    end
    TCleanTCond16R(ismember(TCleanTCond16R.Muscle,Muscle),:) = TMusR;
end

for iMus = 1:length(FlexExtMuscles)
    Muscle = FlexExtMuscles{iMus};
    TMus = TCleanTCond16R(ismember(TCleanTCond16R.Muscle,Muscle),:);
    TMusR = TMus;
    for iSample = 1:height(TMus)
        TSample = TMus(iSample,:);
        SwTimeNorm = (TSample.CycleTime - TSample.StanceTime)/TSample.CycleTime;
        SwTime = (TSample.CycleTime - TSample.StanceTime);

        TMusR(iSample,:).PCNorm = SwTimeNorm;
        TMusR(iSample,:).PC_sec = SwTime;
              
        AllBurstOn = TSample.All_BurstOn;
        AllBurstOff = TSample.All_BurstOff;
        AllBurstOnNorm = TSample.AllBurstOnNorm;
        AllBurstOffNorm = TSample.AllBurstOffNorm;
        AllBurstAve = TSample.All_BurstAve;
        AllBurstAveNormMax = TSample.AllBurstAveNormMax;
        AllBurstAveNormPeakEnv = TSample.AllBurstAveNormPeakEnv;
         
        All_BurstAvePrePC = TSample.All_BurstAvePrePC;
        All_BurstAvePostPC = TSample.All_BurstAvePostPC;


        TMusR(iSample,:).FBurstOn = cellfun(@(x,y,z) z(x>MinFlexOn&x<MaxFlexOn&y<SwTimeNorm+MaxFlexOffAfterPC),AllBurstOnNorm,AllBurstOffNorm,AllBurstOn,'UniformOutput', false);
        TMusR(iSample,:).FBurstOff = cellfun(@(x,y,z) z(x>MinFlexOn&x<MaxFlexOn&y<SwTimeNorm+MaxFlexOffAfterPC),AllBurstOnNorm,AllBurstOffNorm,AllBurstOff,'UniformOutput', false);
        TMusR(iSample,:).FBurstOnNorm = cellfun(@(x,y,z) z(x>MinFlexOn&x<MaxFlexOn&y<SwTimeNorm+MaxFlexOffAfterPC),AllBurstOnNorm,AllBurstOffNorm,AllBurstOnNorm,'UniformOutput', false);
        TMusR(iSample,:).FBurstOffNorm = cellfun(@(x,y,z) z(x>MinFlexOn&x<MaxFlexOn&y<SwTimeNorm+MaxFlexOffAfterPC),AllBurstOnNorm,AllBurstOffNorm,AllBurstOffNorm,'UniformOutput', false);
        TMusR(iSample,:).FBurstAve = cellfun(@(x,y,z) z(x>MinFlexOn&x<MaxFlexOn&y<SwTimeNorm+MaxFlexOffAfterPC),AllBurstOnNorm,AllBurstOffNorm,AllBurstAve,'UniformOutput', false);
        TMusR(iSample,:).FBurstAveNormMax = cellfun(@(x,y,z) z(x>MinFlexOn&x<MaxFlexOn&y<SwTimeNorm+MaxFlexOffAfterPC),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormMax,'UniformOutput', false);
        TMusR(iSample,:).FBurstAveNormPeakEnv = cellfun(@(x,y,z) z(x>MinFlexOn&x<MaxFlexOn&y<SwTimeNorm+MaxFlexOffAfterPC),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormPeakEnv,'UniformOutput', false);

        TMusR(iSample,:).EBurstOn = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOn,'UniformOutput', false);
        TMusR(iSample,:).EBurstOff = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOff,'UniformOutput', false);
        TMusR(iSample,:).EBurstOnNorm = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOnNorm,'UniformOutput', false);
        TMusR(iSample,:).EBurstOffNorm = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOffNorm,'UniformOutput', false);
        TMusR(iSample,:).EBurstAve = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstAve,'UniformOutput', false);
        TMusR(iSample,:).EBurstAveNormMax = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormMax,'UniformOutput', false);
        TMusR(iSample,:).EBurstAveNormPeakEnv = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstAveNormPeakEnv,'UniformOutput', false);

        TMusR(iSample,:).EBurstAvePrePC = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,All_BurstAvePrePC,'UniformOutput', false);
        TMusR(iSample,:).EBurstAvePostPC = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,All_BurstAvePostPC,'UniformOutput', false);

        EBurstOnNorm = cellfun(@(x,y) x(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,'UniformOutput', false);
        EBurstOffNorm = cellfun(@(x,y) y(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,'UniformOutput', false);
         
        EBurstOn = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOn,'UniformOutput', false);
        EBurstOff = cellfun(@(x,y,z) z(x>MinExtOn&y>SwTimeNorm+MinExOffAfterPC&y<MaxExtOff),AllBurstOnNorm,AllBurstOffNorm,AllBurstOff,'UniformOutput', false);

        TMusR(iSample,:).EBurstPrePCtimeNorm =  cellfun(@(x,y) x - y,num2cell(SwTimeNorm),EBurstOnNorm,'UniformOutput', false);
        TMusR(iSample,:).EBurstPostPCtimeNorm =  cellfun(@(x,y) - x + y,num2cell(SwTimeNorm),EBurstOffNorm,'UniformOutput', false);
        
        TMusR(iSample,:).EBurstPrePCtime_sec =cellfun(@(x,y) x - y,num2cell(SwTime),EBurstOn,'UniformOutput', false);
        TMusR(iSample,:).EBurstPostPCtime_sec =  cellfun(@(x,y) - x + y,num2cell(SwTime),EBurstOff,'UniformOutput', false);



     end
    TCleanTCond16R(ismember(TCleanTCond16R.Muscle,Muscle),:) = TMusR;
end

save(TableSpeedsFilename,'TCleanTCond16','TCleanTCond16R','-append')
