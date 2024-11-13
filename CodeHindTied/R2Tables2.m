clear variables
DataFolder = '..\DataHindTied\';

Nreplicates = 20;
NTime= 101;
NShuffles = 50;
Conditions = {'INTACT','SPINAL'};
TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';
TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB.mat');

load(TableSpeedsFilename,'TCleanTCond15R')
MaxNsyn = 14;
TVars = {'Condition','Speed','Syn','R2'};
TR2AllSyns = cell2table(cell(0,length(TVars)),'VariableNames',TVars);
for N_syn = 1:MaxNsyn
    TR2 = cell2table(cell(0,length(TVars)),'VariableNames',TVars);
    for iCond = 1:length(Conditions)
        Condition = Conditions{iCond};%'SPINAL';%'INTACT';%

        TClean = TCleanTCond15R(ismember(TCleanTCond15R.Condition,Condition),:);

        SideSpeeds = unique(TClean.SideSpeed);%Speed of treadmill of side of limb

        for iSideSpeed = 1:length(SideSpeeds)
            SideSpeed = SideSpeeds{iSideSpeed};

            % possible speeds of opposite belt for choosen this side speed
            ContraSideSpeeds = unique(TClean(ismember(TClean.SideSpeed,SideSpeed),:).ContraSideSpeed);

            for iContraSideSpeed = 1:length(ContraSideSpeeds)
                ContraSideSpeed = ContraSideSpeeds{iContraSideSpeed};
                % one condition is unique combination of SideSpeed and ContraSideSpeed
                TCond = TClean(ismember(TClean.SideSpeed,SideSpeed)&...
                    ismember(TClean.ContraSideSpeed,ContraSideSpeed),:);

                CondMuscles = unique(TCond.Muscle);
                NCondMuscles = length(CondMuscles);

                f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',num2str(N_syn),'_SideSpeed-',SideSpeed,...
                    '_ContraSideSpeed-',ContraSideSpeed,'_NShuffles-',...
                    num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
                load(f_name,'R2')
                NSamples = length(R2);
                ConditonV = repmat({Condition},NSamples,1);
                SpeedV = repmat({SideSpeed},NSamples,1);
                SynV = repmat(num2cell(N_syn),NSamples,1);
                T1 = cell2table([ConditonV,SpeedV,SynV,num2cell(R2)],'VariableNames',TVars);
                TR2 = [TR2;T1];
            end
        end
    end
    f_name = strcat(DataFolder,'CondAllSpeedsAll_',TreadmillCondition,'_Syn-',num2str(N_syn),...
        '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
    load(f_name,'R2')
    NSamples = length(R2);
    ConditonV = repmat({'ComCond'},NSamples,1);
    SpeedV = repmat({'ComSpeed'},NSamples,1);
    SynV = repmat(num2cell(N_syn),NSamples,1);
    T1 = cell2table([ConditonV,SpeedV,SynV,num2cell(R2)],'VariableNames',TVars);
    TR2 = [TR2;T1];

    f_nameT =  strcat(DataFolder,'R2_CondAllSpeedsAll_',TreadmillCondition,'_Syn-',num2str(N_syn),...
        '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.txt');
    writetable(TR2,f_nameT);
    TR2AllSyns = [TR2AllSyns;TR2];
end
f_nameTAll =  strcat(DataFolder,'R2_CondAllSpeedsAll_',TreadmillCondition,'_AllSyn',...
    '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.txt');
writetable(TR2AllSyns,f_nameTAll);

