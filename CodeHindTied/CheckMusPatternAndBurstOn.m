clear variables
DataFolder = '..\DataHindTied\';

NTime= 101;
Conditions = {'INTACT','SPINAL'};
TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';
TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB.mat');
load(TableSpeedsFilename,'TCleanTCond15R')

Speeds2Check = {'0.4','0.7','1.0'};

Muscles2Check = {'BFA','CF','FDL','GLU','LG','MG','PLA','PLO','SOL','BFP'};
t = 0:100;
for iMus = 1:length(Muscles2Check)
    Muscle = Muscles2Check{iMus};
        
    f = figure('Name',strcat(Muscle,'.Pattern'));
    f.Position = [282         195        1172         701];
    for iSpeed = 1:length(Speeds2Check)
        Speed2Check = Speeds2Check{iSpeed};

        TBurstSpeedMus = TCleanTCond15R(ismember(TCleanTCond15R.Muscle,Muscle)&ismember(TCleanTCond15R.SideSpeed,Speed2Check),:);

        for iCond = 1:length(Conditions)
            Condition = Conditions{iCond};
            subplot(length(Speeds2Check),length(Conditions),iCond + (iSpeed-1)*length(Conditions))
            TBurstSpeedMusCond = TBurstSpeedMus(ismember(TBurstSpeedMus.Condition,Condition),:);

            plot(t,cell2mat(TBurstSpeedMusCond.EMG_Norm'))
            hold on
            PC = 100*(TBurstSpeedMusCond.CycleTime - TBurstSpeedMusCond.StanceTime)./TBurstSpeedMusCond.CycleTime;
            BurstsOn = cell2mat(TBurstSpeedMusCond.EBurstOnNorm)*100;
            xline(BurstsOn)
            xline(mean(PC),'r--','LineWidth',5)
            title(strcat(Muscle,'.',Condition,'.',Speed2Check))
            text(80,0.8,strcat('B_O_n = ',num2str(mean(BurstsOn),2)));
            hold off
            xlim([0,100])
        end
    end
    ax = gcf;
    FigName = strcat('../Figures/Mus_',Muscle,'.pdf');%,'_',Condition,'_',Speed2Check
    exportgraphics(ax,FigName,'BackgroundColor','none','ContentType','vector')

end
