clear variables
close all
DataFolder = '..\DataHindTied\';

Conditions = {'INTACT','SPINAL'};
TreadmillCondition = 'TIED';

[TClean,CleanMuscles,NCleanMuscles,N_Allmuscles,SideSpeeds,NSideSpeeds,...
    MusMap,MusclesNew] = loadAllSpeedsAllCondTable2(DataFolder,TreadmillCondition);

MusNames_str = char(join(CleanMuscles,','));

N_syn = 5;
NShuffles = 50;
NTime = 101;
TempSortByPeaks = false;% true;
NHbase_repeat = 10;%  How many time repeat nnmf to find base of synergies Hbase/Wbase
NBins =10; % for stat tables
doSave = true;%false;%
doSavePlots = true;%false;%
Nskip = NShuffles; %NShuffles for mean only

TW_vars = {'Condition','SideSpeed','ContraSideSpeed','Syn','Muscle','W'};
TH_vars = {'Condition','SideSpeed','ContraSideSpeed','Syn','tBin','H'};

%
[doPlot,doPlotEach,doDispR2,pos,FontSz,FontSzLabel,CapSz,...
    FontWeightLabel,NConditions,CindCond,...
    CindCondInBetween,CtCond,Cx2] = SetPlotParametersAllCondAllSPeed2(TreadmillCondition,NTime);
%
%adjustment of order of synergies
Hmaps = repmat(1:N_syn,NConditions+1+2+3,1);
if N_syn ==8 &&N_Allmuscles==15 &&NShuffles ==50
    Hmaps(1,:) = circshift(1:N_syn,1);%0.4, intact
    Hmaps(2,:) = circshift(1:N_syn,1);%0.7, intact
    Hmaps(3,:) = circshift(1:N_syn,1);%1.0, intact
    Hmaps(4,:) = circshift(1:N_syn,1);%0.4, spinal
    %Hmaps(5,:) = circshift(1:N_syn,1);%0.7, spinal

    Hmaps(8,:) = [1,2,8,3,4,5,6,7];%all speeds,intact
    %Hmaps(9,:) = circshift(1:N_syn,1);%all speeds, Spinal
    % 8,1,2,3,4,5,6,7
    % 6,1,7,8,2,3,4,5
    Hmaps(9,:) = [6,1,7,8,2,3,4,5];%all speeds, Spinal

    Hmaps(10,:) = circshift(1:N_syn,1);%0.4, all conditions
    %    Hmaps(11,:) = [2,3,5,4,6,7,8,1];%0.7, all conditions
    Hmaps(11,:) = [8,3,1,2,5,4,6,7];%0.7, all conditions
end
if N_syn ==5 &&N_Allmuscles==15 &&NShuffles ==50
    Hmaps(4,:) = circshift(1:N_syn,1);%0.4, Spinal
    Hmaps(9,:) = [4,1,2,3,5];%all speeds, Spinal
    Hmaps(10,:) = circshift(1:N_syn,1);%0.4, all conditions
end
if N_syn ==5 &&N_Allmuscles==16 &&NShuffles ==1000
    Hmaps(4,:) = circshift(1:N_syn,1);%0.4, Spinal
    Hmaps(5,:) = circshift(1:N_syn,1); %0.7,Spinal
    %Hmaps(6,:) = circshift(1:N_syn,1); %1,0,Spinal
    %     Hmaps(NConditionsPlus,:) = circshift(1:N_syn,1);
    %     %Hmaps(NConditionsPlus,:) = [1,2,3,5,4];%[5,4,3,1,2];%[1,3,4,2,5]; %all Speeds, Spinal
    %Hmaps(NConditionsPlus,:) = circshift(1:N_syn,1);
end
if N_syn ==7
    if strcmp(Condition,'INTACT')
        Hmaps(1,:) = circshift(1:N_syn,1);%0.4, INTACT
        Hmaps(2,:) = circshift(1:N_syn,1); %0.7,INTACT
        Hmaps(3,:) = circshift(1:N_syn,1); %1,0,INTACT

        % Hmaps(NConditionsPlus,:) = [1,2,4,5,3]; %all speeds intact
    else
        Hmaps(1,:) = circshift(1:N_syn,1);%0.4, Spinal
        Hmaps(2,:) = circshift(1:N_syn,1); %0.7,Spinal
        %Hmaps(3,:) = circshift(1:N_syn,1); %1,0,Spinal
        %     Hmaps(NConditionsPlus,:) = circshift(1:N_syn,1);
        %     %Hmaps(NConditionsPlus,:) = [1,2,3,5,4];%[5,4,3,1,2];%[1,3,4,2,5]; %all Speeds, Spinal
    end
    Hmaps(NConditionsPlus,:) = circshift(1:N_syn,1);
end


%Spacial
HmapsS = repmat(1:N_syn,NConditions+1+2+3,1);


emg_set_tm_mN = nan(NTime,N_Allmuscles,NConditions);
Sw2Plots = nan(1,NConditions);
Nmus_expected = N_Allmuscles;%length(unique(CTO{1,iSel}.Muscle));
NCondMuscles = N_Allmuscles;
NTimeAll = NTime*NConditions;

%temporal sorted
W_mD = zeros(N_syn,Nmus_expected,NConditions);
W_sdD = zeros(N_syn,Nmus_expected,NConditions);
H_mD = zeros(NTimeAll,N_syn,NConditions);
H_sdD = zeros(NTimeAll,N_syn,NConditions);

W_5D = zeros(N_syn,Nmus_expected,NShuffles,NConditions);
H_5D = zeros(NTimeAll,N_syn,NShuffles,NConditions);

W_4D = zeros(N_syn,Nmus_expected,NShuffles,NConditions);
H_4D = zeros(NTimeAll,N_syn,NShuffles,NConditions);

W_6D = zeros(N_syn,Nmus_expected,NShuffles,NConditions);
H_6D = zeros(NTimeAll,N_syn,NShuffles,NConditions);


%Spacially sorted
W_mDs = zeros(N_syn,Nmus_expected,NConditions);
W_sdDs = zeros(N_syn,Nmus_expected,NConditions);
H_mDs = zeros(NTimeAll,N_syn,NConditions);
H_sdDs = zeros(NTimeAll,N_syn,NConditions);

Cases = cell(1,NConditions);
TW = cell2table(cell(0,length(TW_vars)),'VariableNames',TW_vars);
TH = cell2table(cell(0,length(TH_vars)),'VariableNames',TH_vars);

TWs = cell2table(cell(0,length(TW_vars)),'VariableNames',TW_vars);
THs = cell2table(cell(0,length(TH_vars)),'VariableNames',TH_vars);

iCase = 1;% will be increamenting across two nested cycles to get every
% combination of SideSpeed and ContraSideSpeed

Hbases  = nan(NTime,N_syn,NConditions); %holds Hbase for each condition
Hbases2  = nan(NTime,N_syn,NConditions); %holds Hbase for each condition
Hbases3  = nan(NTime,N_syn,NConditions); %holds Hbase for each condition

hConds = [];
for iCond = 1:length(Conditions)
    Condition = Conditions{iCond};
    TCond = TClean(ismember(TClean.Condition,Condition),:);
    for iSideSpeed = 1:length(SideSpeeds)
        SideSpeed = SideSpeeds{iSideSpeed};
        % possible speeds of opposite belt for choosen this side speed
        ContraSideSpeeds = unique(TCond(ismember(TCond.SideSpeed,SideSpeed),:).ContraSideSpeed);

        for iContraSideSpeed = 1:length(ContraSideSpeeds)
            ContraSideSpeed = ContraSideSpeeds{iContraSideSpeed};

            CondName = strcat(Condition,'.',SideSpeed,'-',ContraSideSpeed);
            Cases{iCase} = CondName;
            NameStr = strcat(CondName,'-syn',num2str(N_syn));
            % NameStrS = strcat(CondName,'-syn',num2str(N_syn),'Spacial sort');

            [TCase,CondMuscles,MusMapCond,NCondAllMuscles,H_all,W_all,V_all,R2,Rmse]...
                = loadOneSpeedsTable5(DataFolder,Condition,TreadmillCondition,...
                N_syn,SideSpeed,ContraSideSpeed,NShuffles,N_Allmuscles,TCond,...
                CleanMuscles,MusclesNew,NCondMuscles,CondName);

            [Sw2Plot,SwSD2Plot] = GetSwTime4Plots(TCase);

            % find Hbase for temp sorting
            [Hbase, Hbase_set,MV,MI,Hbase_setSorted,R2Hbase,RmseHbase] = ...
                findBasicTempPatterns3(H_all,NTime,N_syn,NHbase_repeat,NameStr);
            %[Hbase, Hbase_set,MV,MI] = findBasicTempPatterns2(H_all,NTime,N_syn,NHbase_repeat,NameStr);
            %Hbases(:,:,iCase) = Hbase;

            %find Wbase for spacial sorting
            [Wbase,Wbase_set] = findBasicSpacialPatterns(W_all,NCondMuscles,N_syn,NHbase_repeat);

            H_all3D_unsorted = reshape(permute(H_all,[3,2,1]),NTime,N_syn,[]); %make rectangular matrix with all patterns
            W_all3D_unsorted = reshape(permute(W_all,[3,2,1]),N_syn,NCondMuscles,[]); %make rectangular matrix with all patterns

            NTime_method = 0;

            [HbaseS,H_all3D_s,W_all3D_s] = sortByTempBasePatterns(Hbase,N_syn,...
                Hmaps(iCase,:),TempSortByPeaks,NShuffles,H_all3D_unsorted,W_all3D_unsorted,NTime,NTime_method);
            %
            % W_5D = zeros(N_syn,Nmus_expected,NShuffles,NConditions);
            % H_5D = zeros(NTimeAll,N_syn,NShuffles,NConditions);
            %Hbases(:,:,iCase) = HbaseS;

            W_5D(:,:,:,iCase) = W_all3D_s;
            H_5D(1:NTime,:,:,iCase) = H_all3D_s;


            [WbaseS,H_all3D_ss,W_all3D_ss] = sortBySpacialBasePatterns(Wbase,N_syn,...
                NShuffles,H_all3D_unsorted,W_all3D_unsorted,HmapsS(iCase,:));


            [TW,TH,TWs,THs,TWS1D,THS1D,TWS1Ds,THS1Ds] = fillStatTables5(TW_vars,...
                TH_vars,TW,TH,TWs,THs,NShuffles,N_syn,NCondMuscles,CondMuscles,...
                W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,Condition,SideSpeed,...
                ContraSideSpeed,NBins,NTime,SideSpeeds,Conditions);

            % prep to plot, shift H numbers based on Condition/CondName
            [emg_set_tm_m,emg_set_tm_sd,F_all,F_alls,F_tm_m,F_tm_sd,F_tms_m,...
                F_tms_sd,W_sm_m,W_sm_sd,W_sms_m,W_sms_sd,Hs_m,Hs_sd,Hss_m,Hss_sd,...
                H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs]...
                = prepPlots(V_all,NTime,NCondMuscles, NShuffles,N_Allmuscles,MusMapCond,...
                W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,SideSpeed,ContraSideSpeed,iCase,...
                H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs);%,Hmaps

            % plot
            [hConds(iCase),Sw2Plots] = plotMusPattSyn2(CleanMuscles,NameStr,pos,MusMapCond,...
                emg_set_tm_m,emg_set_tm_sd,Sw2Plot,iCase,NTime,MusclesNew,R2,Rmse,...
                F_tm_m,F_tm_sd,SwSD2Plot,N_syn,H_mD,H_sdD,H_all3D_s,doPlotEach,...
                W_mD,W_sdD,N_Allmuscles,Sw2Plots,1,doSavePlots);

            iCase = iCase+1;
        end
    end
end

% combined all speed synergies
CondName = 'AllCond';
SideSpeed = 'AllSpeeds';
ContraSideSpeed = SideSpeed;
Cases{iCase} = CondName;
NameStr = strcat(CondName,'-syn',num2str(N_syn),'_sh',num2str(NShuffles));
NameStrS = strcat(CondName,'-syn',num2str(N_syn),'_sh',num2str(NShuffles),'Spacial sort');

[TCase,CondMuscles,MusMapCond,NCondAllMuscles,H_all,W_all,V_all,R2,Rmse]...
    = loadOneSpeedsTable5(DataFolder,Condition,TreadmillCondition,...
    N_syn,SideSpeed,ContraSideSpeed,NShuffles,N_Allmuscles,TClean,...
    CleanMuscles,MusclesNew,NCondMuscles,CondName);

[Sw2Plot,SwSD2Plot] = GetSwTime4Plots(TCase);

[Hbase, Hbase_set,MV,MI,Hbase_setSorted,R2Hbase,RmseHbase] = findBasicTempPatterns3(H_all,NTimeAll,N_syn,NHbase_repeat,NameStr);

[Wbase,Wbase_set] = findBasicSpacialPatterns(W_all,NCondMuscles,N_syn,NHbase_repeat);

H_all3D_unsorted = reshape(permute(H_all,[3,2,1]),NTimeAll,N_syn,[]); %make rectangular matrix with all patterns
W_all3D_unsorted = reshape(permute(W_all,[3,2,1]),N_syn,NCondMuscles,[]); %make rectangular matrix with all patterns


NTime_method = 6; % 0- whole time, 1- 1st part, 2- mean of 3 parts,6 - mean of 6 parts
[HbaseS,H_all3D_s,W_all3D_s] = sortByTempBasePatterns(Hbase,N_syn,...
    Hmaps(iCase,:),TempSortByPeaks,NShuffles,H_all3D_unsorted,W_all3D_unsorted,NTime,NTime_method);

W_5D(:,:,:,iCase) = W_all3D_s;
H_5D(1:NTimeAll,:,:,iCase) = H_all3D_s;

[WbaseS,H_all3D_ss,W_all3D_ss] = sortBySpacialBasePatterns(Wbase,N_syn,...
    NShuffles,H_all3D_unsorted,W_all3D_unsorted,HmapsS(iCase,:));
% tables before adding combined all speed/all conditions synergies
TW_WoAll = TW;
TH_WoAll = TH;
[TW,TH,TWs,THs,TWS1D,THS1D,TWS1Ds,THS1Ds] = ...
    fillStatTables5(TW_vars,TH_vars,TW,TH,TWs,THs,NShuffles,N_syn,...
    NCondMuscles,CondMuscles,W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,...
    CondName,SideSpeed,ContraSideSpeed,NBins,NTime,SideSpeeds,Conditions);

% prep to plot
[emg_set_tm_m,emg_set_tm_sd,F_all,F_alls,F_tm_m,F_tm_sd,F_tms_m,...
    F_tms_sd,W_sm_m,W_sm_sd,W_sms_m,W_sms_sd,Hs_m,Hs_sd,Hss_m,Hss_sd,...
    H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs]...
    = prepPlots(V_all,NTimeAll,NCondMuscles, NShuffles,N_Allmuscles,MusMapCond,...
    W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,SideSpeed,ContraSideSpeed,iCase,...
    H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs);%,Hmaps

% plot
[hConds(iCase),Sw2Plots] = plotMusPattSyn2(CleanMuscles,NameStr,pos,MusMapCond,...
    emg_set_tm_m,emg_set_tm_sd,Sw2Plot,iCase,NTimeAll,MusclesNew,R2,Rmse,...
    F_tm_m,F_tm_sd,SwSD2Plot,N_syn,H_mD,H_sdD,H_all3D_s,doPlotEach,...
    W_mD,W_sdD,N_Allmuscles,Sw2Plots,NTimeAll/NTime,doSavePlots);
iCase = iCase+1;

% v3 adding new combinations, 1 - all speeds one condition
SideSpeed = 'AllSpeeds';
ContraSideSpeed = SideSpeed;
for iCond = 1:length(Conditions)
    Condition = Conditions{iCond};
    CaseName =  strcat('AllSpeedsCond',Condition);
    %SideSpeed =strcat('AllSpeedsCond',Condition);
    %ContraSideSpeed = SideSpeed;
    Cases{iCase} = CaseName;
    NameStr = strcat(CaseName,'-syn',num2str(N_syn),'_sh',num2str(NShuffles));
    NameStrS = strcat(CaseName,'-syn',num2str(N_syn),'_sh',num2str(NShuffles),'Spacial sort');

    [TCase,CondMuscles,MusMapCond,NCondAllMuscles,H_all,W_all,V_all,R2,Rmse]...
        = loadOneSpeedsTable5(DataFolder,Condition,TreadmillCondition,...
        N_syn,SideSpeed,ContraSideSpeed,NShuffles,N_Allmuscles,TClean,...
        CleanMuscles,MusclesNew,NCondMuscles,Condition);

    [Sw2Plot,SwSD2Plot] = GetSwTime4Plots(TCase);

    [Hbase, Hbase_set,MV,MI,Hbase_setSorted,R2Hbase,RmseHbase] = ...
        findBasicTempPatterns3(H_all,NTime*length(SideSpeeds),N_syn,NHbase_repeat,NameStr);

    [Wbase,Wbase_set] = findBasicSpacialPatterns(W_all,NCondMuscles,N_syn,NHbase_repeat);

    H_all3D_unsorted = reshape(permute(H_all,[3,2,1]),NTime*length(SideSpeeds),N_syn,[]); %make rectangular matrix with all patterns
    W_all3D_unsorted = reshape(permute(W_all,[3,2,1]),N_syn,NCondMuscles,[]); %make rectangular matrix with all patterns

    NTime_method = length(SideSpeeds); % 0- whole time, 1- 1st part, 2/3- mean of 2/3 parts,6 - mean of 6 parts

    [HbaseS,H_all3D_s,W_all3D_s] = sortByTempBasePatterns(Hbase,N_syn,...
        Hmaps(iCase,:),TempSortByPeaks,NShuffles,H_all3D_unsorted,W_all3D_unsorted,NTime,NTime_method);

    W_5D(:,:,:,iCase) = W_all3D_s;
    H_5D(1:NTime*length(SideSpeeds),:,:,iCase) = H_all3D_s;

    [WbaseS,H_all3D_ss,W_all3D_ss] = sortBySpacialBasePatterns(Wbase,N_syn,...
        NShuffles,H_all3D_unsorted,W_all3D_unsorted,HmapsS(iCase,:));

    [TW,TH,TWs,THs,TWS1D,THS1D,TWS1Ds,THS1Ds] = ...
        fillStatTables5(TW_vars,TH_vars,TW,TH,TWs,THs,NShuffles,N_syn,...
        NCondMuscles,CondMuscles,W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,...
        Condition,SideSpeed,ContraSideSpeed,NBins,NTime,SideSpeeds,Conditions);
    % prep to plot
    [emg_set_tm_m,emg_set_tm_sd,F_all,F_alls,F_tm_m,F_tm_sd,F_tms_m,...
        F_tms_sd,W_sm_m,W_sm_sd,W_sms_m,W_sms_sd,Hs_m,Hs_sd,Hss_m,Hss_sd,...
        H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs]...
        = prepPlots(V_all,NTime*length(SideSpeeds),NCondMuscles, NShuffles,N_Allmuscles,MusMapCond,...
        W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,SideSpeed,ContraSideSpeed,iCase,...
        H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs);%,Hmaps

    % plot
    [hConds(iCase),Sw2Plots] = plotMusPattSyn3(CleanMuscles,NameStr,pos,MusMapCond,...
        emg_set_tm_m,emg_set_tm_sd,Sw2Plot,iCase,NTime*length(SideSpeeds),MusclesNew,R2,Rmse,...
        F_tm_m,F_tm_sd,SwSD2Plot,N_syn,H_mD,H_sdD,H_all3D_s,doPlotEach,...
        W_mD,W_sdD,N_Allmuscles,Sw2Plots,length(SideSpeeds),doSavePlots);

    iCase = iCase+1;
end
% all Conditions one speed
Condition = 'AllCond';
for iSideSpeed = 1:length(SideSpeeds)
    SideSpeed = SideSpeeds{iSideSpeed};
    ContraSideSpeed = SideSpeed;

    CaseName = strcat('AllCondSpeed',SideSpeed);
    Cases{iCase} = CaseName;
    NameStr = strcat(CaseName,'-syn',num2str(N_syn),'_sh',num2str(NShuffles));
    %    NameStrS = strcat(CaseName,'-syn',num2str(N_syn),'_sh',num2str(NShuffles),'Spacial sort');

    [TCase,CondMuscles,MusMapCond,NCondAllMuscles,H_all,W_all,V_all,R2,Rmse]...
        = loadOneSpeedsTable5(DataFolder,Condition,TreadmillCondition,...
        N_syn,SideSpeed,ContraSideSpeed,NShuffles,N_Allmuscles,TClean,...
        CleanMuscles,MusclesNew,NCondMuscles,Condition);

    [Sw2Plot,SwSD2Plot] = GetSwTime4Plots(TCase);

    [Hbase, Hbase_set,MV,MI,Hbase_setSorted,R2Hbase,RmseHbase] = ...
        findBasicTempPatterns3(H_all,NTime*length(Conditions),N_syn,NHbase_repeat,NameStr);

    [Wbase,Wbase_set] = findBasicSpacialPatterns(W_all,NCondMuscles,N_syn,NHbase_repeat);

    H_all3D_unsorted = reshape(permute(H_all,[3,2,1]),NTime*length(Conditions),N_syn,[]); %make rectangular matrix with all patterns
    W_all3D_unsorted = reshape(permute(W_all,[3,2,1]),N_syn,NCondMuscles,[]); %make rectangular matrix with all patterns

    NTime_method = length(Conditions); % 0- whole time, 1- 1st part, 2/3- mean of 2/3 parts,6 - mean of 6 parts

    [HbaseS,H_all3D_s,W_all3D_s] = sortByTempBasePatterns(Hbase,N_syn,...
        Hmaps(iCase,:),TempSortByPeaks,NShuffles,H_all3D_unsorted,W_all3D_unsorted,NTime,NTime_method);

    W_5D(:,:,:,iCase) = W_all3D_s;
    H_5D(1:NTime*length(Conditions),:,:,iCase) = H_all3D_s;

    [WbaseS,H_all3D_ss,W_all3D_ss] = sortBySpacialBasePatterns(Wbase,N_syn,...
        NShuffles,H_all3D_unsorted,W_all3D_unsorted,HmapsS(iCase,:));

    [TW,TH,TWs,THs,TWS1D,THS1D,TWS1Ds,THS1Ds] = ...
        fillStatTables5(TW_vars,TH_vars,TW,TH,TWs,THs,NShuffles,N_syn,...
        NCondMuscles,CondMuscles,W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,...
        Condition,SideSpeed,ContraSideSpeed,NBins,NTime,SideSpeeds,Conditions);
    % prep to plot
    [emg_set_tm_m,emg_set_tm_sd,F_all,F_alls,F_tm_m,F_tm_sd,F_tms_m,...
        F_tms_sd,W_sm_m,W_sm_sd,W_sms_m,W_sms_sd,Hs_m,Hs_sd,Hss_m,Hss_sd,...
        H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs]...
        = prepPlots(V_all,NTime*length(Conditions),NCondMuscles, NShuffles,N_Allmuscles,MusMapCond,...
        W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,SideSpeed,ContraSideSpeed,iCase,...
        H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs);%,Hmaps

    % plot
    [hConds(iCase),Sw2Plots] = plotMusPattSyn3(CleanMuscles,NameStr,pos,MusMapCond,...
        emg_set_tm_m,emg_set_tm_sd,Sw2Plot,iCase,NTime*length(Conditions),MusclesNew,R2,Rmse,...
        F_tm_m,F_tm_sd,SwSD2Plot,N_syn,H_mD,H_sdD,H_all3D_s,doPlotEach,...
        W_mD,W_sdD,N_Allmuscles,Sw2Plots,length(Conditions),doSavePlots);

    iCase = iCase+1;

end

%% combined figure
FigName = strcat('_CombAllCondAllSpeed_syn_',num2str(N_syn));
UsedConditions = 1:6;
CombCondition = 7;
[hComb] = plotCombineFigure3(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);
%
FigName = strcat('_CombIntactAllSpeed_syn_',num2str(N_syn));
UsedConditions = 1:3;
CombCondition = 8;
[hComb] = plotCombineFigure3(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);

%
FigName = strcat('_CombSpinalAllSpeed_syn_',num2str(N_syn));
UsedConditions = 4:6;
CombCondition = 9;
[hComb] = plotCombineFigure3(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);

%
FigName = strcat('_CombAllCondSpeed0.4_syn_',num2str(N_syn));
UsedConditions = [1,4];
CombCondition = 10;
[hComb] = plotCombineFigure3(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);

%
FigName = strcat('_CombAllCondSpeed0.7_syn_',num2str(N_syn));
UsedConditions = [2,5];
CombCondition = 11;
[hComb] = plotCombineFigure3(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);
%
FigName = strcat('_CombAllCondSpeed1.0_syn_',num2str(N_syn));
UsedConditions = [3,6];
CombCondition = 12;
[hComb] = plotCombineFigure3(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);
%
FigName = strcat('_CombAllCond-0.4-0.7-1.0_syn_',num2str(N_syn));
UsedConditions = [10,11,12];
CombCondition = [];
[hComb] = plotCombineFigure4(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);
%
FigName = strcat('_CombINTACT-SPINAL-AllSpeeds_syn_',num2str(N_syn));
UsedConditions = [8,9];
CombCondition = [];
[hComb] = plotCombineFigure4(FigName,N_syn,H_5D,W_5D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition);


%%
if doSave
    NameStr = strcat('AllCondAllSpeeds_syn-',num2str(N_syn),'_sh',num2str(NShuffles));
    StatTablesFile = strcat('../Data/Tables_AllCondAllSpeeds_syn-',num2str(N_syn),'_sh',num2str(NShuffles),'.mat');
    save(StatTablesFile,'TH','TW','TH_WoAll','TW_WoAll')
    writetable(TH,strcat('../Data/TH_',NameStr,'.txt'))
    writetable(TW,strcat('../Data/TW_',NameStr,'.txt'))
    writetable(TH_WoAll,strcat('../Data/TH_WoAll_',NameStr,'.txt'))
    writetable(TW_WoAll,strcat('../Data/TW_WoAll_',NameStr,'.txt'))

    SortedSynFile = strcat('../Data/Sorted_syn-',num2str(N_syn),'_sh',num2str(NShuffles),'.mat');
    save(SortedSynFile,'TH','TW','W_5D','H_5D','Cases','Conditions',...
        'MusMap','MusclesNew','CleanMuscles','NTimeAll','NTime')

end

