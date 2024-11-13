function [hComb] = plotCombineFigure4(FigName,N_syn,H_4D,W_4D,...
    Nskip,Sw2Plots,SwSD2Plot,NTime,NShuffles,CtCond,CindCond,...
    CindCondInBetween,Cx2,FontWeightLabel,FontSzLabel,FontSz,CapSz,MusclesNew,...
    MusMap,doSavePlots,UsedConditions,CombCondition)

Nmus_expected = length(MusclesNew);
NConditions = length(UsedConditions);
switch NConditions
    case 2 % for each condition, three speeds combined
        Ncycles = 3;
    case 3 % for each speed, two conditions combined
        Ncycles = 2;
end
NTimeAll = NTime*Ncycles;

hComb =figure('Name',FigName);
%hComb =figure('Name',strcat('syn_',num2str(N_syn)));
pos = 1000*[0.01,    0,    1,    0.5];
set(hComb,'position',pos);
NCol = 3;
pos_sH = zeros(N_syn,4);
pos_sW = zeros(N_syn,4);

darkBrown = [0.8500, 0.3250, 0.0980];
lightBrown = [0.9290, 0.6940, 0.1250];
lightBlue = [91, 207, 244] / 255;
lightGreen = [0.7608    0.9294    0.7333];%[91, 244, 207] / 255;
lightRed = [1.0000    0.5686    0.5686];% [244,91, 207] / 255;

ColorNames = {'Blue','Green','Red','Brown'};
ColorD = [{'b','g','r'},repmat({lightBrown},1,12-NConditions)];
ColorInd = [{[0,0,0.5],[0,0.5,0],[0.5,0,0]},repmat({darkBrown},1,12-NConditions)];

for i = 1:N_syn
    % W_5D = zeros(N_syn,Nmus_expected,NShuffles,NConditions);
    % H_5D = zeros(NTimeAll,N_syn,NShuffles,NConditions);

    %HSyn = squeeze(H_4D(:,i,:,:));
    HSyn = squeeze(H_4D(:,i,:,UsedConditions));

    Hm = squeeze(mean(HSyn,2)); %across all Shuffles
    Hsd = squeeze(std(HSyn,[],2));

    Hsh = HSyn(:,1:Nskip:end,:);%every Nskip of shuffled data

    HmPlusSD = Hm+Hsd;
    HmMinusSD = Hm-Hsd;

    inBetweenPlus = [Hm',fliplr(HmPlusSD')]';
    inBetweenMinus = [HmMinusSD',fliplr(Hm')]';

    hAxisH(i) = subplot(N_syn,NCol,NCol*i);% for H
    pos_sH(i,:) = get( hAxisH(i), 'Position' );

    v_offset = 0.01;%for plot black around mean
    Hmp = Hm+v_offset;
    Hmm = Hm-v_offset;

    hold on
    for iCond = 1:NConditions
        UsedCondition = UsedConditions(iCond);
        switch UsedCondition
            case 8
                UsedSingleConditions = 1:3;
            case 9
                UsedSingleConditions = 4:6;
            case 10
                UsedSingleConditions = [1,4];
            case 11
                UsedSingleConditions = [2,5];
            case 12
                UsedSingleConditions = [3,6];
        end
        for iCyc =1:Ncycles
            %swing onset
            Sw2Plot = Sw2Plots(UsedSingleConditions(iCyc)) + NTime*(iCyc-1);
            line([Sw2Plot Sw2Plot],[0 1],'Color',ColorD{iCond},'LineStyle','-','LineWidth',1)
            % line([Sw2Plot+SwSD2Plot Sw2Plot+SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',.5)
            % line([Sw2Plot-SwSD2Plot Sw2Plot-SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',.5)
    
            %cycles separation
            line([NTime*(iCyc-1) NTime*(iCyc-1)],[0 1],'Color','k','LineStyle','-','LineWidth',1)
        end
    end

    ConditionsID = UsedConditions;
    for iCond = 1:NConditions
        iCondID = ConditionsID(iCond);
        itCondID = iCond;
        if iCond >NConditions
            itCondID = iCondID;
        end
        if Nskip~=NShuffles
            plot(CtCond{iCondID},squeeze(Hsh(CindCond{iCondID},:,iCond)),'color',ColorD{iCond},'LineWidth',0.01);
        end

        plot(CtCond{iCondID},Hm(CindCond{iCondID},iCond),'Color',ColorInd{iCond},'LineWidth',0.6);

        plot1 = fill(Cx2{iCondID}, inBetweenPlus(CindCondInBetween{iCondID},iCond),ColorD{iCond});
        plot1.FaceAlpha = .2;
        plot1.EdgeColor ='none';

        plot2 = fill(Cx2{iCondID}, inBetweenMinus(CindCondInBetween{iCondID},iCond),ColorD{iCond});
        plot2.FaceAlpha = .2;
        plot2.EdgeColor ='none';
    end

    ylabel(strcat(['C_' num2str(i)]),'fontweight',FontWeightLabel, 'FontSize', FontSzLabel);
    set(gca,'box','off')
    if(i~=N_syn)
        %             set(gca,'XColor','none');
        set(gca,'XTickLabel',{[]});
    end
    set(gca, 'FontSize', FontSz,'TickDir','out');
    axis([0 NTimeAll 0 1.2]);

    hold off

    %WSyn = squeeze(W_4D(i,MusMap,:,:));
    WSyn = squeeze(W_4D(i,MusMap,:,UsedConditions));

    Wm = squeeze(mean(WSyn,2)); %across all Shuffles
    Wsd = squeeze(std(WSyn,[],2));

    Wsh = WSyn(:,1:Nskip:end,:);%every Nskip of shuffled data

    hAxisW(i) = subplot(N_syn,NCol,[NCol*(i-1)+1,NCol*(i-1)+2]);% for W
    pos_sW(i,:) = get( hAxisW(i), 'Position' );

    b= bar(Wm,'FaceColor','none');

    set(b,'BarWidth',1);    % The bars will now touch each other
    hold on;

    [numgroups,numbars] = size(Wm);%N_mus, N_dir
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for j = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        %errorbar(x, Wm(MusMap,j), Wsd(MusMap,j), 'k', 'linestyle', 'none','CapSize',CapSz);
        errorbar(x, Wm(:,j), Wsd(:,j), 'k', 'linestyle', 'none','CapSize',CapSz);
    end
    set(gca,'box','off','TickDir','out')
    %         if(i~=N_syn)
    %             set(gca,'XColor','none');
    %         end
    axis([0.5 Nmus_expected+0.5 0 1.2]);
    ylabel(strcat(['W_' num2str(i)]),'fontweight',FontWeightLabel, 'FontSize', FontSzLabel)
    set(gca,'XTick',linspace(1,Nmus_expected,Nmus_expected),'XTickLabel',MusclesNew, 'FontSize', FontSz);

    %shuffles
    if Nskip~=NShuffles
        Wsh2 = reshape(Wsh,Nmus_expected,[]);
        bsh= bar(squeeze(Wsh2));
        for iCond =1:NConditionsPlus
            set(bsh(1+(iCond-1)*NShuffles/Nskip:iCond*NShuffles/Nskip),'FaceColor',ColorD{iCond});
        end
        set(bsh,'BarWidth',1);
    end
    b= bar(Wm);
    for iCond = 1:NConditions % [UsedConditions,CombCondition] %
        %        set(b(iW),'FaceColor','k','EdgeColor',[1 1 1]);
        if Nskip~=NShuffles
            set(b(iCond),'FaceColor','none','EdgeColor','k','LineWidth',0.5);
        else
            set(b(iCond),'FaceColor',ColorD{iCond},'EdgeColor','k','LineWidth',0.5);
        end
    end
    set(b,'BarWidth',1);    % The bars will now touch each other

    hold off;
end

% adjust space between synergy plot panels
BetweenS = 0.01;
Total_height = pos_sH(1,2)+pos_sH(1,end)- pos_sH(N_syn,2);%NConditionsPlus
NewHeight = (Total_height-N_syn*BetweenS)/N_syn;
for iSyn = N_syn-1:-1:1
    pos_sH(iSyn,2) = pos_sH(iSyn+1,2) + NewHeight+BetweenS; %shift down
    pos_sW(iSyn,2) = pos_sW(iSyn+1,2) + NewHeight+BetweenS; %shift down
end
pos_sH(:,4) = repmat(NewHeight,N_syn,1); %increase height ,5
pos_sW(:,4) = repmat(NewHeight,N_syn,1); %increase height ,5
for iSyn = 1:N_syn
    set( hAxisH(iSyn), 'Position', pos_sH(iSyn,:));
    set( hAxisW(iSyn), 'Position', pos_sW(iSyn,:));
end

if N_syn<10
    sgtitle(FigName(2:end-6))
else
    sgtitle(FigName(2:end-7))
end

if doSavePlots
    combFigName = strcat('../Figures/Fig_',FigName,'_Selected_',num2str(NShuffles/Nskip),'.pdf');
    print(combFigName,'-dpdf','-bestfit','-r0','-painters')%);
end