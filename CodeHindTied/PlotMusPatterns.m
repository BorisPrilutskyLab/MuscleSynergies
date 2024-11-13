clear variables
close all
DataFolder = '..\DataHindTied\';
Conditions = {'INTACT','SPINAL'};
TreadmillCondition = 'TIED';

[TClean,CleanMuscles,NCleanMuscles,N_Allmuscles,SideSpeeds,NSideSpeeds,...
    MusMap,MusclesNew] = loadAllSpeedsAllCondTable2(DataFolder,TreadmillCondition);

grayColor = [.7 .7 .7];
BeltSpeeds = unique(TClean.SideSpeed);
t = linspace(0,100,101);

StanceOnsetFrames = nan(length(BeltSpeeds),length(Conditions));
for iSpeed = 1:length(BeltSpeeds)
    BeltSpeed = BeltSpeeds(iSpeed);
    for iCond = 1:length(Conditions)
        Condition = Conditions(iCond);
        TCase = TClean(ismember(TClean.SideSpeed,BeltSpeed)&...
            ismember(TClean.Condition,Condition),:);
        StanceOnsetFrames(iSpeed,iCond) = round((1 - mean(TCase.DF))*100);
    end
end


figure('Name','Muscles')
for iMus = 1:length(MusclesNew)
    Muscle = MusclesNew(iMus);
    for iSpeed = 1:length(BeltSpeeds)
        BeltSpeed = BeltSpeeds(iSpeed);
        for iCond = 1:length(Conditions)
            Condition = Conditions(iCond);
            TMusCase = TClean(ismember(TClean.Muscle,Muscle)&...
                ismember(TClean.SideSpeed,BeltSpeed)&...
                ismember(TClean.Condition,Condition),:);
            CMusEMG = TMusCase.EMG_Norm2;

            %mean across muscles
            StanceOnsetFrame = StanceOnsetFrames(iSpeed,iCond);

            DMusEMG = cell2mat(CMusEMG')';
            mMusEMG = mean(DMusEMG);
            %sdMusEMG = std(DMusEMG);
            subplot(length(MusclesNew),length(BeltSpeeds)*length(Conditions),...
                (iMus-1)*length(BeltSpeeds)*length(Conditions) + iCond + (iSpeed-1)*length(Conditions))
            plot(t,DMusEMG','color',grayColor,'LineWidth',.5)
            hold on
            plot(t,mMusEMG,'color','k','LineWidth',1)
            xline(StanceOnsetFrame)
            yline(0)
            hold off
            text(101,.5,Muscle{1,1})
            if iMus==1
                title(strcat(Condition{1,1},'.',BeltSpeed{1,1}))
            end
            if iMus~=length(MusclesNew)
                set(gca,'XTick',[]);
            end
            set(gca,'YTick',[]);
            ylim([0,1])
            box off
        end
    end
end

FigName3 = strcat('../Figures/MusPatterns_4.pdf');
f3 = gcf;
exportgraphics(f3,FigName3,'ContentType','vector')

FigName3e = strcat('../Figures/MusPatterns_4.eps');
exportgraphics(f3,FigName3e,'ContentType','vector')
