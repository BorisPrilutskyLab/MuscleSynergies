clear variables
DataFolder = '..\DataHindTied\';

Nreplicates = 20;
NTime= 101;
NShuffles = 50;


Conditions = {'INTACT', 'SPINAL'};
TreadmillCondition = 'TIED';
SideSpeeds = {'0.4','0.7','1.0'};
TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB.mat');
load(TableSpeedsFilename,'TCleanTCond15R')
NCondMuscles = length(unique(TCleanTCond15R.Muscle));


for N_syn = 1:14
    V_allConds_allSpeeds = [];
    TCondSpeeds = [];
    for iSideSpeed = 1:length(SideSpeeds)
        V_AllCond_OneSpeeds = [];
        TAllCondOneSpeeds = [];
        SideSpeed = SideSpeeds{iSideSpeed};
        ContraSideSpeed = SideSpeed;
        for iCond =1:length(Conditions)
            Condition = Conditions{iCond};

            f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',num2str(N_syn),'_SideSpeed-',SideSpeed,...
                '_ContraSideSpeed-',ContraSideSpeed,'_NShuffles-',...
                num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
            load(f_name)
            % loaded: 'R2','Rmse','V_all','W_all','H_all','CleanMuscles'
            % 'F_all','R2EachMus','RmseEachMus','CondMuscles','TCond'
            %   CondMuscles'

            % all conditions/ all speeds
            V_allConds_allSpeeds = cat(3,V_allConds_allSpeeds,V_all);
            TCondSpeeds = [TCondSpeeds;TCond];

            % all conditions/ one speeds
            V_AllCond_OneSpeeds = cat(3,V_AllCond_OneSpeeds,V_all);
            TAllCondOneSpeeds = [TAllCondOneSpeeds;TCond];
        end
        % all conditions/ one speeds
        FileNameStr = strcat(DataFolder,TreadmillCondition,'CondAllSpeed',SideSpeed);
        GetSynShuffles(V_AllCond_OneSpeeds,NTime,Conditions,SideSpeeds(iSideSpeed),NCondMuscles,N_syn,...
            NShuffles,Nreplicates,FileNameStr,TAllCondOneSpeeds)
    end

    % all conditions/ all speeds combined
    FileNameStr = strcat(DataFolder,TreadmillCondition,'CondAllSpeedAll');
    GetSynShuffles(V_allConds_allSpeeds,NTime,Conditions,SideSpeeds,NCondMuscles,N_syn,...
        NShuffles,Nreplicates,FileNameStr,TCondSpeeds)

    for iCond =1:length(Conditions)
        Condition = Conditions{iCond};
        V_OneCond_AllSpeeds = [];
        TOneCondAllSpeeds = [];
        for iSideSpeed = 1:length(SideSpeeds)
            SideSpeed = SideSpeeds{iSideSpeed};
            ContraSideSpeed = SideSpeed;
            f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',num2str(N_syn),'_SideSpeed-',SideSpeed,...
                '_ContraSideSpeed-',ContraSideSpeed,'_NShuffles-',...
                num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
            load(f_name)

            % all conditions/ one speeds
            V_OneCond_AllSpeeds = cat(3,V_OneCond_AllSpeeds,V_all);
            TOneCondAllSpeeds = [TOneCondAllSpeeds;TCond];
        end
        FileNameStr = strcat(DataFolder,TreadmillCondition,'Cond',Condition,'SpeedAll');
        GetSynShuffles(V_OneCond_AllSpeeds,NTime,Conditions(iCond),SideSpeeds,NCondMuscles,N_syn,...
            NShuffles,Nreplicates,FileNameStr,TOneCondAllSpeeds)
    end
end