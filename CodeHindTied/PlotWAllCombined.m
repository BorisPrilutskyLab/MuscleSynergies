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
pos = 1000*[0.01,    0,    1,    0.9];
Muscles = CleanMuscles;
N_Allmuscles = length(Muscles);
NameStr = 'Combined Synegies';
h =figure('Name',NameStr);
set(h,'position',pos)
NCol = 4;
FontSz = 6;% FontSz = 10;
FontSzLabel =12;
CapSz = 1;
FontWeightLabel ='bold';

UsedConditions = 7:12;
doSavePlots = true;
% plot synergies
for iSyn = 1:N_syn
    WSyn = squeeze(W_4D(iSyn,MusMap,:,UsedConditions));
    Wm = squeeze(mean(WSyn,2)); %across all Shuffles
    Wsd = squeeze(std(WSyn,[],2));

    hAxisW(iSyn) = subplot(N_syn,NCol,[NCol*(iSyn-1)+1,NCol*(iSyn-1)+4]);% for W
    pos_sW(iSyn,:) = get( hAxisW(iSyn), 'Position' );

    b= bar(Wm);%,'FaceColor','none'

    set(b,'BarWidth',1);    % The bars will now touch each other
    hold on;

    [numgroups,numbars] = size(Wm);%N_mus, N_dir
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for j = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, Wm(:,j), Wsd(:,j), 'k', 'linestyle', 'none','CapSize',CapSz);
    end
    set(gca,'box','off','TickDir','out')
    axis([0.5 N_Allmuscles+0.5 0 1.2]);
    ylabel(strcat(['W_' num2str(iSyn)]),'fontweight',FontWeightLabel, 'FontSize', FontSzLabel)
    set(gca,'XTick',linspace(1,N_Allmuscles,N_Allmuscles),'XTickLabel',MusclesNew, 'FontSize', FontSz);
    legend(Cases(UsedConditions))
end

if doSavePlots
    f3 = gcf;
    combFigName_PDF = strcat('../Figures/Fig_',NameStr,'.pdf');
    exportgraphics(f3,combFigName_PDF,'ContentType','vector')

    combFigName_EPS = strcat('../Figures/Fig_',NameStr,'.eps');
    exportgraphics(f3,combFigName_EPS,'ContentType','vector')
end
