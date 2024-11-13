clear variables
DataFolder = 'D:\Dropbox (GaTech)\FrigonLab\\ForeLimbEMGsFrigon\Data\';
%DataFolder = 'E:\DropBox(Gatech)\Dropbox (GaTech)\FrigonLab\ForeLimbEMGsFrigon\Data\';

Nreplicates = 20;
NTime= 101;
NShuffles = 50;
Conditions = {'INTACT','SPINAL'};
TreadmillConditions = {'SPLIT','TIED'};
UsedConfLevelHigh = 0.95;
UsedConfLevelLow = 0.90;
UsedConfLevelLow2 = 0.85;

TreadmillCondition = 'TIED';
TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB.mat');
load(TableSpeedsFilename,'TCleanTCond15R')
MaxNsyn = 14;
mR2 = nan(length(Conditions),3,MaxNsyn);
sdR2 = nan(length(Conditions),3,MaxNsyn);
mR2Com = nan(1,MaxNsyn);
sdR2Com = nan(1,MaxNsyn);
mR2ComSpeed = nan(length(Conditions),MaxNsyn);
sdR2ComSpeed = nan(length(Conditions),MaxNsyn);
mR2ComCond = nan(3,MaxNsyn);
sdR2ComCond = nan(3,MaxNsyn);

SideSpeeds = unique(TCleanTCond15R.SideSpeed);
for N_syn = 1:MaxNsyn
    for iCond = 1:length(Conditions)
        Condition = Conditions{iCond};%'SPINAL';%'INTACT';%
        TClean = TCleanTCond15R(ismember(TCleanTCond15R.Condition,Condition),:);
        %SideSpeeds = unique(TClean.SideSpeed);%Speed of treadmill of side of limb

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
                mR2(iCond,iSideSpeed,N_syn) = mean(R2(:));
                sdR2(iCond,iSideSpeed,N_syn) = std(R2(:));
            end
        end

        FileNameStr = strcat(DataFolder,TreadmillCondition,'Cond',Condition,'SpeedAll');
        f_name = strcat(FileNameStr,'_Syn-',num2str(N_syn),...
            '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
        load(f_name,'R2')
        mR2ComSpeed(iCond,N_syn) = mean(R2(:));
        sdR2ComSpeed(iCond,N_syn) = std(R2(:));

    end
    f_name = strcat(DataFolder,'CondAllSpeedsAll_',TreadmillCondition,'_Syn-',num2str(N_syn),...
        '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
    load(f_name,'R2')
    mR2Com(N_syn) = mean(R2(:));
    sdR2Com(N_syn) = std(R2(:));

    for iSideSpeed = 1:length(SideSpeeds)
        SideSpeed = SideSpeeds{iSideSpeed};
        ContraSideSpeed = SideSpeed;
        TCond = TCleanTCond15R(ismember(TCleanTCond15R.SideSpeed,SideSpeed)&...
            ismember(TCleanTCond15R.ContraSideSpeed,ContraSideSpeed),:);

        CondMuscles = unique(TCond.Muscle);
        NCondMuscles = length(CondMuscles);
        FileNameStr = strcat(DataFolder,TreadmillCondition,'CondAllSpeed',SideSpeed);
        f_name = strcat(FileNameStr,'_Syn-',num2str(N_syn),...
            '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');

        load(f_name,'R2')
        mR2ComCond(iSideSpeed,N_syn) = mean(R2(:));
        sdR2ComCond(iSideSpeed,N_syn) = std(R2(:));
    end
end

%% plot
UsedColors = {'b','g','r'};
ColorNames = {'Blue','Green','Red'};
x0 = (1:MaxNsyn)-0.5;
UsedStyles = {'-','--'};
figure('Name','SynergyNumbers')
%plot(mR2Com,'color','k','LineStyle', ':')
errorbar(x0+0.5,mR2Com,sdR2Com,'color','k','LineStyle', ':')
hold on
yline(UsedConfLevelHigh)
text(.5,UsedConfLevelHigh + 0.02,num2str(UsedConfLevelHigh,2))
yline(UsedConfLevelLow,'--')
text(0.5,UsedConfLevelLow + 0.02,num2str(UsedConfLevelLow,2))
for iCond = 1:length(Conditions)
    Condition = Conditions{iCond};%'SPINAL';%'INTACT';%
    TClean = TCleanTCond15R(ismember(TCleanTCond15R.Condition,Condition),:);
    SideSpeeds = unique(TClean.SideSpeed);%Speed of treadmill of side of limb

    for iSideSpeed = 1:length(SideSpeeds)
        SideSpeed = SideSpeeds{iSideSpeed};
        %            plot(squeeze(mR2(iCond,iSideSpeed,:)),'color',UsedColors{iSideSpeed},'LineStyle', UsedStyles{iCond})
        x = x0 + .5*(iCond-1) + 0.1*(iSideSpeed-1);
        errorbar(x,squeeze(mR2(iCond,iSideSpeed,:)),squeeze(sdR2(iCond,iSideSpeed,:)),'color',UsedColors{iSideSpeed},'LineStyle', UsedStyles{iCond})
    end
    % x = x0 + 1.5 + .1*(iCond-1);
    % errorbar(x,squeeze(mR2ComSpeed(iCond,:)),squeeze(sdR2ComSpeed(iCond,:)),'color','k','LineStyle', UsedStyles{iCond})
end
hold off
%%
   R2FigName = strcat('../Figures/Fig_sh50_v1R2m15_1-14.pdf');
    print(R2FigName,'-dpdf','-bestfit','-r0','-painters')%);
%%
UsedColors = {'b','g','r'};
ColorNames = {'Blue','Green','Red'};
x0 = (1:MaxNsyn)-0.5;
UsedStyles = {'-','--'};
figure('Name','SynergyNumbersComb')
%plot(mR2Com,'color','k','LineStyle', ':')
errorbar(x0+0.5,mR2Com,sdR2Com,'color','k','LineStyle', ':')
hold on
yline(UsedConfLevelHigh)
text(.5,UsedConfLevelHigh + 0.02,num2str(UsedConfLevelHigh,2))
yline(UsedConfLevelLow,'--')
text(0.5,UsedConfLevelLow + 0.02,num2str(UsedConfLevelLow,2))
for iCond = 1:length(Conditions)
    Condition = Conditions{iCond};%'SPINAL';%'INTACT';%
    x = x0 + .1*(iCond-1);
    errorbar(x,squeeze(mR2ComSpeed(iCond,:)),squeeze(sdR2ComSpeed(iCond,:)),'color',UsedColors{iCond},'LineStyle', '--')
end
for iSideSpeed = 1:length(SideSpeeds)
    SideSpeed = SideSpeeds{iSideSpeed};
    x = x0 + .1*length(Conditions) + 0.1*(iSideSpeed-1);
    errorbar(x,squeeze(mR2ComCond(iSideSpeed,:)),squeeze(sdR2ComCond(iSideSpeed,:)),'color',UsedColors{iSideSpeed},'LineStyle', '-')
end
hold off

%%
CapSz = 1;
lightBlue = [[0.6 0.8 1]]; 
darkGreen = [0.4667,    0.6745,    0.1882];
lightRed = [0.9686    0.6078    0.5529];
grey = [0.6510    0.6510    0.6510];
ExraColors = [lightBlue;darkGreen;lightRed];

figure('Name','SynergyNumbersBar')
mR2_r = reshape(permute(mR2,[3,1,2]),MaxNsyn,[]);
sdR2_r = reshape(permute(sdR2,[3,1,2]),MaxNsyn,[]);

mR2_rCom = mR2_r;% [mR2_r,mR2Com'];
sdR2_rCom = sdR2_r;% [sdR2_r,sdR2Com'];

[numgroups,numbars] = size(mR2_rCom);%N_mus, N_dir
groupwidth = min(0.8, numbars/(numbars+1.5));


bsh= bar(mR2_rCom);
hold on
yline(UsedConfLevelHigh)
text(.5,UsedConfLevelHigh + 0.02,num2str(UsedConfLevelHigh,2))
yline(UsedConfLevelLow,'--')
text(0.5,UsedConfLevelLow + 0.02,num2str(UsedConfLevelLow,2))
yline(UsedConfLevelLow2,'.')
text(0.5,UsedConfLevelLow2 + 0.02,num2str(UsedConfLevelLow2,2))

for iSideSpeed = 1:length(SideSpeeds)
    set(bsh(1+(iSideSpeed-1)*2:iSideSpeed*2),'FaceColor',UsedColors{iSideSpeed});
end
set(bsh(2),'FaceColor',lightBlue);
set(bsh(3),'FaceColor',darkGreen);
set(bsh(6),'FaceColor',lightRed);
%set(bsh(7),'FaceColor',grey);

set(bsh,'BarWidth',1);
for j = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x,mR2_rCom(:,j),sdR2_rCom(:,j), 'k', 'linestyle', 'none','CapSize',CapSz)
end

hold off
set(gca,'box','off','TickDir','out')
axis([0.5 MaxNsyn+0.5 0 1.2]);
ylabel(strcat('mean R2'))%,'fontweight',FontWeightLabel, 'FontSize', FontSzLabel
set(gca,'XTick',linspace(1,MaxNsyn,MaxNsyn));%,'XTickLabel',MusclesNew, 'FontSize', FontSz
title('intact(dark)/spinal(light); blue- 0.4, green-0.7, red - 1.0') %; grey - common
%%
   R2FigName = strcat('../Figures/Fig_SynNumberInd_sh50_R2m15_1-14.pdf');
    print(R2FigName,'-dpdf','-bestfit','-r0','-painters')%);
%%
figure('Name','SynergyNumbersCombBar')
mR2ComCond_r = reshape(mR2ComCond',MaxNsyn,[]);
sdR2ComCond_r = reshape(sdR2ComCond',MaxNsyn,[]);

mR2ComSpeed_r = reshape(mR2ComSpeed',MaxNsyn,[]);
sdR2ComSpeed_r = reshape(sdR2ComSpeed',MaxNsyn,[]);

mR2_rCom = [mR2ComCond_r,mR2ComSpeed_r,mR2Com'];
sdR2_rCom = [sdR2ComCond_r,sdR2ComSpeed_r,sdR2Com'];

[numgroups,numbars] = size(mR2_rCom);%N_mus, N_dir
groupwidth = min(0.8, numbars/(numbars+1.5));


bsh= bar(mR2_rCom);
hold on
yline(UsedConfLevelHigh)
text(.5,UsedConfLevelHigh + 0.02,num2str(UsedConfLevelHigh,2))
yline(UsedConfLevelLow,'--')
text(0.5,UsedConfLevelLow + 0.02,num2str(UsedConfLevelLow,2))
yline(UsedConfLevelLow2,'.')
text(0.5,UsedConfLevelLow2 + 0.02,num2str(UsedConfLevelLow2,2))

for iSideSpeed = 1:length(SideSpeeds)
    set(bsh(iSideSpeed),'FaceColor',UsedColors{iSideSpeed});
end
set(bsh(4),'FaceColor',lightBlue);
set(bsh(5),'FaceColor',darkGreen);
%set(bsh(6),'FaceColor',lightRed);
set(bsh(6),'FaceColor',grey);

set(bsh,'BarWidth',1);
for j = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x,mR2_rCom(:,j),sdR2_rCom(:,j), 'k', 'linestyle', 'none','CapSize',CapSz)
end

hold off

set(gca,'box','off','TickDir','out')
axis([0.5 MaxNsyn+0.5 0 1.2]);
ylabel(strcat('mean R2'))%,'fontweight',FontWeightLabel, 'FontSize', FontSzLabel
set(gca,'XTick',linspace(1,MaxNsyn,MaxNsyn));%,'XTickLabel',MusclesNew, 'FontSize', FontSz
title({'Combined Conditions at speeds: blue- 0.4, green-0.7, red - 1.0;';...
    'Combined velocities at: Intact -light blue,SPINAL- dark green';'grey - common'})
%%
   R2FigName = strcat('../Figures/Fig_SynNumberBarCombined_sh50_R2m15_1-14.pdf');
    print(R2FigName,'-dpdf','-bestfit','-r0','-painters')%);