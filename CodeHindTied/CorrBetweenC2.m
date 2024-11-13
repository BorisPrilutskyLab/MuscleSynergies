clear variables
close all
DataFolder = '..\DataHindTied\';

Conditions = {'INTACT','SPINAL'};
TreadmillCondition = 'TIED';
N_syn =5;
NShuffles =50;
SortedSynFile = strcat('../Data/Sorted_syn-',num2str(N_syn),'_sh',num2str(NShuffles),'.mat');
load(SortedSynFile);
% tmp
W_4D = W_5D;
H_4D = H_5D;
%
R2 = nan(N_syn,NShuffles,27);
Rmse = nan(N_syn,NShuffles,27);

H_SynSh = nan(NTime,N_syn,NShuffles,24);

for iSyn =1:N_syn
    for iShuffle = 1:NShuffles
        H_IntS1 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'INTACT.0.4-0.4')));%size = [NTime,1]
        H_IntS2 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'INTACT.0.7-0.7')));%size = [NTime,1]
        H_IntS3 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'INTACT.1.0-1.0')));%size = [NTime,1]

        H_SpiS1 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'SPINAL.0.4-0.4')));%size = [NTime,1]
        H_SpiS2 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'SPINAL.0.7-0.7')));%size = [NTime,1]
        H_SpiS3 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'SPINAL.1.0-1.0')));%size = [NTime,1]

        H_IntSpeedComS1 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'AllSpeedsCondINTACT')));
        H_IntSpeedComS2 = squeeze(H_4D(1+NTime:2*NTime,iSyn,iShuffle,ismember(Cases,'AllSpeedsCondINTACT')));
        H_IntSpeedComS3 = squeeze(H_4D(1+2*NTime:3*NTime,iSyn,iShuffle,ismember(Cases,'AllSpeedsCondINTACT')));

        H_SpiSpeedComS1 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'AllSpeedsCondSPINAL')));
        H_SpiSpeedComS2 = squeeze(H_4D(1+NTime:2*NTime,iSyn,iShuffle,ismember(Cases,'AllSpeedsCondSPINAL')));
        H_SpiSpeedComS3 = squeeze(H_4D(1+2*NTime:3*NTime,iSyn,iShuffle,ismember(Cases,'AllSpeedsCondSPINAL')));

        H_CondComIntSpeedS1 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'AllCondSpeed0.4')));
        H_CondComIntSpeedS2 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'AllCondSpeed0.7')));
        H_CondComIntSpeedS3 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'AllCondSpeed1.0')));

        H_CondComSpiSpeedS1 = squeeze(H_4D(1+NTime:2*NTime,iSyn,iShuffle,ismember(Cases,'AllCondSpeed0.4')));
        H_CondComSpiSpeedS2 = squeeze(H_4D(1+NTime:2*NTime,iSyn,iShuffle,ismember(Cases,'AllCondSpeed0.7')));
        H_CondComSpiSpeedS3 = squeeze(H_4D(1+NTime:2*NTime,iSyn,iShuffle,ismember(Cases,'AllCondSpeed1.0')));

        H_CondComIntSpeedComS1 = squeeze(H_4D(1:NTime,iSyn,iShuffle,ismember(Cases,'AllCond')));
        H_CondComSpiSpeedComS1 = squeeze(H_4D(1+NTime:2*NTime,iSyn,iShuffle,ismember(Cases,'AllCond')));
        H_CondComIntSpeedComS2 = squeeze(H_4D(1+2*NTime:3*NTime,iSyn,iShuffle,ismember(Cases,'AllCond')));
        H_CondComSpiSpeedComS2 = squeeze(H_4D(1+3*NTime:4*NTime,iSyn,iShuffle,ismember(Cases,'AllCond')));
        H_CondComIntSpeedComS3 = squeeze(H_4D(1+4*NTime:5*NTime,iSyn,iShuffle,ismember(Cases,'AllCond')));
        H_CondComSpiSpeedComS3 = squeeze(H_4D(1+5*NTime:6*NTime,iSyn,iShuffle,ismember(Cases,'AllCond')));
        

        H_SynSh(1:NTime,iSyn,iShuffle,1) = H_IntS1;
        H_SynSh(1:NTime,iSyn,iShuffle,2) = H_IntS2;
        H_SynSh(1:NTime,iSyn,iShuffle,3) = H_IntS3;
        
        H_SynSh(1:NTime,iSyn,iShuffle,4) = H_SpiS1;
        H_SynSh(1:NTime,iSyn,iShuffle,5) = H_SpiS2;
        H_SynSh(1:NTime,iSyn,iShuffle,6) = H_SpiS3;
        
        H_SynSh(1:NTime,iSyn,iShuffle,7) = H_IntSpeedComS1;
        H_SynSh(1:NTime,iSyn,iShuffle,8) = H_IntSpeedComS2;
        H_SynSh(1:NTime,iSyn,iShuffle,9) = H_IntSpeedComS3;
        
        H_SynSh(1:NTime,iSyn,iShuffle,10) = H_SpiSpeedComS1;
        H_SynSh(1:NTime,iSyn,iShuffle,11) = H_SpiSpeedComS2;
        H_SynSh(1:NTime,iSyn,iShuffle,12) = H_SpiSpeedComS3;
        
        H_SynSh(1:NTime,iSyn,iShuffle,13) = H_CondComIntSpeedS1;
        H_SynSh(1:NTime,iSyn,iShuffle,14) = H_CondComIntSpeedS2;
        H_SynSh(1:NTime,iSyn,iShuffle,15) = H_CondComIntSpeedS3;
        
        H_SynSh(1:NTime,iSyn,iShuffle,16) = H_CondComSpiSpeedS1;
        H_SynSh(1:NTime,iSyn,iShuffle,17) = H_CondComSpiSpeedS2;
        H_SynSh(1:NTime,iSyn,iShuffle,18) = H_CondComSpiSpeedS3;

        H_SynSh(1:NTime,iSyn,iShuffle,19) = H_CondComIntSpeedComS1;
        H_SynSh(1:NTime,iSyn,iShuffle,20) = H_CondComIntSpeedComS2;
        H_SynSh(1:NTime,iSyn,iShuffle,21) = H_CondComIntSpeedComS3;
        H_SynSh(1:NTime,iSyn,iShuffle,22) = H_CondComSpiSpeedComS1;
        H_SynSh(1:NTime,iSyn,iShuffle,23) = H_CondComSpiSpeedComS2;
        H_SynSh(1:NTime,iSyn,iShuffle,24) = H_CondComSpiSpeedComS3;

            %  Intact different speeds
        R2(iSyn,iShuffle,1) = corr(H_IntS1, H_IntS2)^2;
        Rmse(iSyn,iShuffle,1) =sqrt(mean((H_IntS1 - H_IntS2).^2));

        R2(iSyn,iShuffle,2) = corr(H_IntS3, H_IntS2)^2;
        Rmse(iSyn,iShuffle,2) =sqrt(mean((H_IntS3 - H_IntS2).^2));

        R2(iSyn,iShuffle,3) = corr(H_IntS1, H_IntS3)^2;
        Rmse(iSyn,iShuffle,3) =sqrt(mean((H_IntS1 - H_IntS3).^2));
            %  Spinal different speeds
        R2(iSyn,iShuffle,4) = corr(H_SpiS1, H_SpiS2)^2;
        Rmse(iSyn,iShuffle,4) =sqrt(mean((H_SpiS1 - H_SpiS2).^2));

        R2(iSyn,iShuffle,5) = corr(H_SpiS3, H_SpiS2)^2;
        Rmse(iSyn,iShuffle,5) =sqrt(mean((H_SpiS3 - H_SpiS2).^2));

        R2(iSyn,iShuffle,6) = corr(H_SpiS1, H_SpiS3)^2;
        Rmse(iSyn,iShuffle,6) =sqrt(mean((H_SpiS1 - H_SpiS3).^2));

        %       one speed, different conditions
        R2(iSyn,iShuffle,7) = corr(H_IntS1, H_SpiS1)^2;
        Rmse(iSyn,iShuffle,7) =sqrt(mean((H_IntS1 - H_SpiS1).^2));

        R2(iSyn,iShuffle,8) = corr(H_IntS2, H_SpiS2)^2;
        Rmse(iSyn,iShuffle,8) =sqrt(mean((H_IntS2 - H_SpiS2).^2));

        R2(iSyn,iShuffle,9) = corr(H_IntS3, H_SpiS3)^2;
        Rmse(iSyn,iShuffle,9) =sqrt(mean((H_IntS3 - H_SpiS3).^2));

        % Intact one speed with part of all speed combined for this speed
        R2(iSyn,iShuffle,10) = corr(H_IntS1, H_IntSpeedComS1)^2;
        Rmse(iSyn,iShuffle,10) =sqrt(mean((H_IntS1 - H_IntSpeedComS1).^2));
 
        R2(iSyn,iShuffle,11) = corr(H_IntS2, H_IntSpeedComS2)^2;
        Rmse(iSyn,iShuffle,11) =sqrt(mean((H_IntS2 - H_IntSpeedComS2).^2));

        R2(iSyn,iShuffle,12) = corr(H_IntS3, H_IntSpeedComS3)^2;
        Rmse(iSyn,iShuffle,12) =sqrt(mean((H_IntS3 - H_IntSpeedComS3).^2));

        % Spinal one speed with part of all speed combined for this speed
        R2(iSyn,iShuffle,13) = corr(H_SpiS1, H_SpiSpeedComS1)^2;
        Rmse(iSyn,iShuffle,13) =sqrt(mean((H_SpiS1 - H_SpiSpeedComS1).^2));

        R2(iSyn,iShuffle,14) = corr(H_SpiS2, H_SpiSpeedComS2)^2;
        Rmse(iSyn,iShuffle,14) =sqrt(mean((H_SpiS2 - H_SpiSpeedComS2).^2));

        R2(iSyn,iShuffle,15) = corr(H_SpiS3, H_SpiSpeedComS3)^2;
        Rmse(iSyn,iShuffle,15) =sqrt(mean((H_SpiS3 - H_SpiSpeedComS3).^2));

        %first speed, combined conditions
        R2(iSyn,iShuffle,16) = corr(H_IntS1, H_CondComIntSpeedS1)^2;
        Rmse(iSyn,iShuffle,16) =sqrt(mean((H_IntS1 - H_CondComIntSpeedS1).^2));
       
        R2(iSyn,iShuffle,17) = corr(H_SpiS1, H_CondComSpiSpeedS1)^2;
        Rmse(iSyn,iShuffle,17) =sqrt(mean((H_SpiS1 - H_CondComSpiSpeedS1).^2));

                %second speed, combined conditions
        R2(iSyn,iShuffle,18) = corr(H_IntS2, H_CondComIntSpeedS2)^2;
        Rmse(iSyn,iShuffle,18) =sqrt(mean((H_IntS2 - H_CondComIntSpeedS2).^2));
       
        R2(iSyn,iShuffle,19) = corr(H_SpiS2, H_CondComSpiSpeedS2)^2;
        Rmse(iSyn,iShuffle,19) =sqrt(mean((H_SpiS2 - H_CondComSpiSpeedS2).^2));

                %third speed, combined conditions
        R2(iSyn,iShuffle,20) = corr(H_IntS3, H_CondComIntSpeedS3)^2;
        Rmse(iSyn,iShuffle,20) =sqrt(mean((H_IntS3 - H_CondComIntSpeedS3).^2));

        R2(iSyn,iShuffle,21) = corr(H_SpiS3, H_CondComSpiSpeedS3)^2;
        Rmse(iSyn,iShuffle,21) =sqrt(mean((H_SpiS3 - H_CondComSpiSpeedS3).^2));

        % all cond
        R2(iSyn,iShuffle,22) = corr(H_IntS1, H_CondComIntSpeedComS1)^2;
        Rmse(iSyn,iShuffle,22) =sqrt(mean((H_IntS1 - H_CondComIntSpeedComS1).^2));

        R2(iSyn,iShuffle,23) = corr(H_IntS2, H_CondComIntSpeedComS2)^2;
        Rmse(iSyn,iShuffle,23) =sqrt(mean((H_IntS2 - H_CondComIntSpeedComS2).^2));

        R2(iSyn,iShuffle,24) = corr(H_IntS3, H_CondComIntSpeedComS3)^2;
        Rmse(iSyn,iShuffle,24) =sqrt(mean((H_IntS3 - H_CondComIntSpeedComS3).^2));

        R2(iSyn,iShuffle,25) = corr(H_SpiS1, H_CondComSpiSpeedComS1)^2;
        Rmse(iSyn,iShuffle,25) =sqrt(mean((H_SpiS1 - H_CondComSpiSpeedComS1).^2));

        R2(iSyn,iShuffle,26) = corr(H_SpiS2, H_CondComSpiSpeedComS2)^2;
        Rmse(iSyn,iShuffle,26) =sqrt(mean((H_SpiS2 - H_CondComSpiSpeedComS2).^2));

        R2(iSyn,iShuffle,27) = corr(H_SpiS3, H_CondComSpiSpeedComS3)^2;
        Rmse(iSyn,iShuffle,27) =sqrt(mean((H_SpiS3 - H_CondComSpiSpeedComS3).^2));

    end
end
%
R2Cases = {'IntactS1&IntactS2','IntactS3&IntactS2','IntactS1&IntactS3',...
    'SpinalS1&SpinalS2','SpinalS3&SpinalS2','SpinalS1&SpinalS3',...
    'IntactS1&SpinalS1','IntactS2&SpinalS2','IntactS3&SpinalS3',...
    'IntactS1&IntactCombinedSpeedS1','IntactS2&IntactCombinedSpeedS2',...
    'IntactS3&IntactCombinedSpeedS3','SpinalS1&SpinalCombinedSpeedS1',...
    'SpinalS2&SpinalCombinedSpeedS2','SpinalS3&SpinalCombinedSpeedS3',...
    'IntactS1&CombinedCondIntactS1','SpinalS1&CombinedCondSpinalS1',...
    'IntactS2&CombinedCondIntactS2','SpinalS2&CombinedCondSpinalS2',...
    'IntactS3&CombinedCondIntactS3','SpinalS3&CombinedCondSpinalS3',...
    'IntactS1&AllCombinedIntactS1','IntactS2&AllCombinedIntactS2',...
    'IntactS3&AllCombinedIntactS3','SpinalS1&AllCombinedSpinalS1',...
    'SpinalS2&AllCombinedSpinalS2','SpinalS3&AllCombinedSpinalS3'};

TR2Vars = {'Syn','Case','R2','RMSE'};
TR2 = cell2table(cell(0,length(TR2Vars)),'VariableNames',TR2Vars);

TmR2Vars = {'Syn','Case','meanR2','sdR2','meanRMSE','sdRMSE'};
TmR2 = cell2table(cell(0,length(TmR2Vars)),'VariableNames',TmR2Vars);

for iCase =1:length(R2Cases)
    for iSyn = 1:N_syn
        CR2_1 = [repmat({iSyn},NShuffles,1),repmat(R2Cases(iCase),NShuffles,1),...
            num2cell(R2(iSyn,:,iCase)'),num2cell(Rmse(iSyn,:,iCase)')];
        TR2_1 = cell2table(CR2_1,'VariableNames',TR2Vars);
        TR2 = [TR2;TR2_1];

        CmR2_1 = [{iSyn},R2Cases(iCase),{mean(R2(iSyn,:,iCase),2)},...
            {std(R2(iSyn,:,iCase),[],2)},{mean(Rmse(iSyn,:,iCase),2)},...
            {std(Rmse(iSyn,:,iCase),[],2)}];
        TmR2_1 = cell2table(CmR2_1,'VariableNames',TmR2Vars);
        TmR2 = [TmR2;TmR2_1];
    end
end
%
writetable(TR2,'../Data/R2BetweenC');
writetable(TmR2,'../Data/meanR2BetweenC');


%%
ColorD = {'b','g','r'};
t = linspace(0,100,101);
for iSyn = 1:N_syn
    subplot(N_syn,1,iSyn)
    hold on
    for iCase = 1:3
        mH = mean(H_SynSh(:,iSyn,:,iCase),3);
        plot(t,mH,'Color',ColorD{iCase},'LineStyle','-')
        mH = mean(H_SynSh(:,iSyn,:,iCase+3),3);
        plot(t,mH,'Color',ColorD{iCase},'LineStyle','--')
    end

    hold off
    ylim([0,1])
    xlim([0,100])
    R2_Int_str1 = strcat('R^2_I_N_T_0_._4_-_I_N_T_0_._7 = ',num2str(mean(R2(iSyn,:,1)),2));
    R2_Int_str2 = strcat('R^2_I_N_T_0_._7_-_I_N_T_1_._0 = ',num2str(mean(R2(iSyn,:,2)),2));
    R2_Int_str3 = strcat('R^2_I_N_T_0_._4_-_I_N_T_1_._0 = ',num2str(mean(R2(iSyn,:,3)),2));

    R2_Spi_str1 = strcat('R^2_S_P_I_0_._4_-_S_P_I_0_._7 = ',num2str(mean(R2(iSyn,:,4)),2));
    R2_Spi_str2 = strcat('R^2_S_P_I_0_._7_-_S_P_I_1_._0 = ',num2str(mean(R2(iSyn,:,5)),2));
    R2_Spi_str3 = strcat('R^2_S_P_I_0_._4_-_S_P_I_1_._0 = ',num2str(mean(R2(iSyn,:,6)),2));


    text(0,.9,{R2_Int_str1;R2_Int_str2;R2_Int_str3})
    text(80,.9,{R2_Spi_str1;R2_Spi_str2;R2_Spi_str3})
end
%%
FigNameC_PDF = strcat('../Figures/CorrR2IntactSpinal.pdf');
f1 = gcf;
exportgraphics(f1,FigNameC_PDF,'ContentType','vector')

FigNameC_EPS = strcat('../Figures/CorrR2IntactSpinal.eps');
exportgraphics(f1,FigNameC_EPS,'ContentType','vector')

%%
figure('Name','INTACTCom')
for iSyn = 1:N_syn
    subplot(N_syn,1,iSyn)
    hold on
    for iCase = 1:3
        mH = mean(H_SynSh(:,iSyn,:,iCase),3);
        plot(t,mH,'Color',ColorD{iCase},'LineStyle','-')
        mH = mean(H_SynSh(:,iSyn,:,iCase+6),3);
        plot(t,mH,'Color',ColorD{iCase},'LineStyle','--')
    end

    hold off
    ylim([0,1])
    xlim([0,100])
    R2_Int_str1 = strcat('R^2_I_N_T_0_._4_-_I_N_T_0_._7 = ',num2str(mean(R2(iSyn,:,1)),2));
    R2_Int_str2 = strcat('R^2_I_N_T_0_._7_-_I_N_T_1_._0 = ',num2str(mean(R2(iSyn,:,2)),2));
    R2_Int_str3 = strcat('R^2_I_N_T_0_._4_-_I_N_T_1_._0 = ',num2str(mean(R2(iSyn,:,3)),2));

    R2_IntCom_str1 = strcat('R^2_I_N_T_0_._4_-_C_O_M_-_I_N_T_0_._4 = ',num2str(mean(R2(iSyn,:,10)),2));
    R2_IntCom_str2 = strcat('R^2_I_N_T_0_._7_-_C_O_M_-_I_N_T_0_._7 = ',num2str(mean(R2(iSyn,:,11)),2));
    R2_IntCom_str3 = strcat('R^2_I_N_T_1_._0_-_C_O_M_-_I_N_T_1_._0 = ',num2str(mean(R2(iSyn,:,12)),2));


    text(0,.9,{R2_Int_str1;R2_Int_str2;R2_Int_str3})
    text(80,.9,{R2_IntCom_str1;R2_IntCom_str2;R2_IntCom_str3})
end
%%
FigNameC_PDF = strcat('../Figures/CorrR2IntactCom.pdf');
f1 = gcf;
exportgraphics(f1,FigNameC_PDF,'ContentType','vector')

FigNameC_EPS = strcat('../Figures/CorrR2IntactCom.eps');
exportgraphics(f1,FigNameC_EPS,'ContentType','vector')

%%

figure('Name','SPINALCom')
for iSyn = 1:N_syn
    subplot(N_syn,1,iSyn)
    hold on
    for iCase = 1:3
        mH = mean(H_SynSh(:,iSyn,:,iCase+3),3);
        plot(t,mH,'Color',ColorD{iCase},'LineStyle','-')
        mH = mean(H_SynSh(:,iSyn,:,iCase+9),3);
        plot(t,mH,'Color',ColorD{iCase},'LineStyle','--')
    end

    hold off
    ylim([0,1])
    xlim([0,100])
    R2_Spi_str1 = strcat('R^2_S_P_I_0_._4_-_S_P_I_0_._7 = ',num2str(mean(R2(iSyn,:,4)),2));
    R2_Spi_str2 = strcat('R^2_S_P_I_0_._7_-_S_P_I_1_._0 = ',num2str(mean(R2(iSyn,:,5)),2));
    R2_Spi_str3 = strcat('R^2_S_P_I_0_._4_-_S_P_I_1_._0 = ',num2str(mean(R2(iSyn,:,6)),2));

    R2_SpiCom_str1 = strcat('R^2_S_P_I_0_._4_-_C_O_M_-_S_P_I_0_._4 = ',num2str(mean(R2(iSyn,:,13)),2));
    R2_SpiCom_str2 = strcat('R^2_S_P_I_0_._7_-_C_O_M_-_S_P_I_0_._7 = ',num2str(mean(R2(iSyn,:,14)),2));
    R2_SpiCom_str3 = strcat('R^2_S_P_I_1_._0_-_C_O_M_-_S_P_I_1_._0 = ',num2str(mean(R2(iSyn,:,15)),2));


    text(0,.9,{R2_Spi_str1;R2_Spi_str2;R2_Spi_str3})
    text(80,.9,{R2_SpiCom_str1;R2_SpiCom_str2;R2_SpiCom_str3})
end
%%
FigNameC_PDF = strcat('../Figures/CorrR2SpinalCom.pdf');
f1 = gcf;
exportgraphics(f1,FigNameC_PDF,'ContentType','vector')

FigNameC_EPS = strcat('../Figures/CorrR2SpinalCom.eps');
exportgraphics(f1,FigNameC_EPS,'ContentType','vector')