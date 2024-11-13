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
% Calculate angles between synergies (W) Fig13(12 in rev1)
A_IntS1S2 = zeros(NShuffles,N_syn);
A_IntS3S2 = zeros(NShuffles,N_syn);
A_IntS1S3 = zeros(NShuffles,N_syn);

A_SpiS1S2 = zeros(NShuffles,N_syn);
A_SpiS3S2 = zeros(NShuffles,N_syn);
A_SpiS1S3 = zeros(NShuffles,N_syn);

A_IntSpiS1 = zeros(NShuffles,N_syn);
A_IntSpiS2= zeros(NShuffles,N_syn);
A_IntSpiS3 = zeros(NShuffles,N_syn);

A_IntS1IntCom = zeros(NShuffles,N_syn);
A_IntS2IntCom = zeros(NShuffles,N_syn);
A_IntS3IntCom = zeros(NShuffles,N_syn);

A_SpiS1SpiCom = zeros(NShuffles,N_syn);
A_SpiS2SpiCom = zeros(NShuffles,N_syn);
A_SpiS3SpiCom = zeros(NShuffles,N_syn);

for iSyn =1:N_syn
    for iShuffle = 1:NShuffles 
    W_IntS1 = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'INTACT.0.4-0.4')));%size = [1,Nmus]
    W_IntS2 = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'INTACT.0.7-0.7')));%size = [1,Nmus]
    W_IntS3 = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'INTACT.1.0-1.0')));%size = [1,Nmus]

    W_SpiS1 = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'SPINAL.0.4-0.4')));%size = [1,Nmus]
    W_SpiS2 = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'SPINAL.0.7-0.7')));%size = [1,Nmus]
    W_SpiS3 = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'SPINAL.1.0-1.0')));%size = [1,Nmus]

    W_IntCom = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'AllSpeedsCondINTACT')));%size = [1,Nmus]
    W_SpiCom = squeeze(W_4D(iSyn,:,iShuffle,ismember(Cases,'AllSpeedsCondSPINAL')));%size = [1,Nmus]

    nW_IntS1 = sqrt(sumsqr(W_IntS1));% norm of W_IntS1
    nW_IntS2 = sqrt(sumsqr(W_IntS2));
    nW_IntS3 = sqrt(sumsqr(W_IntS3));

    nW_SpiS1 = sqrt(sumsqr(W_SpiS1));% norm of W_SpiS1
    nW_SpiS2 = sqrt(sumsqr(W_SpiS2));
    nW_SpiS3 = sqrt(sumsqr(W_SpiS3));

    nW_IntCom = sqrt(sumsqr(W_IntCom));
    nW_SpiCom = sqrt(sumsqr(W_SpiCom));
    
    dot_IntS1S2 = W_IntS1*W_IntS2';
    cos_IntS1S2 = dot_IntS1S2/(nW_IntS1*nW_IntS2);
    A_IntS1S2(iShuffle,iSyn) = acos(cos_IntS1S2)*180/pi;

    dot_IntS3S2 = W_IntS3*W_IntS2';
    cos_IntS3S2 = dot_IntS3S2/(nW_IntS3*nW_IntS2);
    A_IntS3S2(iShuffle,iSyn) = acos(cos_IntS3S2)*180/pi;

    dot_IntS1S3 = W_IntS1*W_IntS3';
    cos_IntS1S3 = dot_IntS1S3/(nW_IntS1*nW_IntS3);
    A_IntS1S3(iShuffle,iSyn) = acos(cos_IntS1S3)*180/pi;

    dot_SpiS1S2 = W_SpiS1*W_SpiS2';
    cos_SpiS1S2 = dot_SpiS1S2/(nW_SpiS1*nW_SpiS2);
    A_SpiS1S2(iShuffle,iSyn) = acos(cos_SpiS1S2)*180/pi;

    dot_SpiS3S2 = W_SpiS3*W_SpiS2';
    cos_SpiS3S2 = dot_SpiS3S2/(nW_SpiS3*nW_SpiS2);
    A_SpiS3S2(iShuffle,iSyn) = acos(cos_SpiS3S2)*180/pi;

    dot_SpiS1S3 = W_SpiS1*W_SpiS3';
    cos_SpiS1S3 = dot_SpiS1S3/(nW_SpiS1*nW_SpiS3);
    A_SpiS1S3(iShuffle,iSyn) = acos(cos_SpiS1S3)*180/pi;

   
    dot_IntSpiS1 = W_IntS1*W_SpiS1';
    cos_IntSpiS1 = dot_IntSpiS1/(nW_IntS1*nW_SpiS1);
    A_IntSpiS1(iShuffle,iSyn) = acos(cos_IntSpiS1)*180/pi;

    dot_IntSpiS2 = W_IntS2*W_SpiS2';
    cos_IntSpiS2 = dot_IntSpiS2/(nW_IntS2*nW_SpiS2);
    A_IntSpiS2(iShuffle,iSyn) = acos(cos_IntSpiS2)*180/pi;

    dot_IntSpiS3 = W_IntS3*W_SpiS3';
    cos_IntSpiS3 = dot_IntSpiS3/(nW_IntS3*nW_SpiS3);
    A_IntSpiS3(iShuffle,iSyn) = acos(cos_IntSpiS3)*180/pi;
    
    dot_IntS1IntCom = W_IntS1*W_IntCom';
    cos_IntS1IntCom = dot_IntS1IntCom/(nW_IntS1*nW_IntCom);
    A_IntS1IntCom(iShuffle,iSyn) = acos(cos_IntS1IntCom)*180/pi;

    dot_IntS2IntCom = W_IntS2*W_IntCom';
    cos_IntS2IntCom = dot_IntS2IntCom/(nW_IntS2*nW_IntCom);
    A_IntS2IntCom(iShuffle,iSyn) = acos(cos_IntS2IntCom)*180/pi;

    dot_IntS3IntCom = W_IntS3*W_IntCom';
    cos_IntS3IntCom = dot_IntS3IntCom/(nW_IntS3*nW_IntCom);
    A_IntS3IntCom(iShuffle,iSyn) = acos(cos_IntS3IntCom)*180/pi;

    dot_SpiS1SpiCom = W_SpiS1*W_SpiCom';
    cos_SpiS1SpiCom = dot_SpiS1SpiCom/(nW_SpiS1*nW_SpiCom);
    A_SpiS1SpiCom(iShuffle,iSyn) = acos(cos_SpiS1SpiCom)*180/pi;

    dot_SpiS2SpiCom = W_SpiS2*W_SpiCom';
    cos_SpiS2SpiCom = dot_SpiS2SpiCom/(nW_SpiS2*nW_SpiCom);
    A_SpiS2SpiCom(iShuffle,iSyn) = acos(cos_SpiS2SpiCom)*180/pi;

    dot_SpiS3SpiCom = W_SpiS3*W_SpiCom';
    cos_SpiS3SpiCom = dot_SpiS3SpiCom/(nW_SpiS3*nW_SpiCom);
    A_SpiS3SpiCom(iShuffle,iSyn) = acos(cos_SpiS3SpiCom)*180/pi;

    end
end
%
CBetweens = {'IntS1S2','IntS3S2','IntS1S3','SpiS1S2','SpiS3S2','SpiS1S3',...
    'IntSpiS1','IntSpiS2','IntSpiS3','IntS1IntCom','IntS2IntCom',...
    'IntS3IntCom','SpiS1SpiCom','SpiS2SpiCom','SpiS3SpiCom'};
TVars = {'Syn','Angle','Between'};
TAngles = cell2table(cell(0,length(TVars)),'VariableNames',TVars);
for iBet = 1:length(CBetweens)
    Bet_m = genvarname(strcat('A_',CBetweens{iBet}));
    eval(['A_m =', Bet_m, ';']);
    
    for iSyn = 1:N_syn
        Syns = repmat({iSyn},NShuffles,1);
        As = num2cell(A_m(:,iSyn));
        Bs = repmat(CBetweens(iBet),NShuffles,1);
        T = table(Syns,As,Bs,'VariableNames',TVars);
        TAngles = [TAngles;T];
    end
end
%
%writetable(TAngles,'../Data/AnglesBetweenW');
%%
doPlotShuffles = false;%true;%
NplotCases = 9;
ColorG = ones(NplotCases,3);
for iCol = 1:NplotCases
    ColorG(iCol,:) = (iCol*(1./(NplotCases+1)))*ones(1,3);
end

figure('Name',strcat('Angles between Conditions W_syn_',num2str(N_syn)));
CapSz = 4;
if NplotCases==6
    Am = [mean(A_IntS1S2)',mean(A_IntS3S2)',mean(A_IntS1S3)',...
        mean(A_SpiS1S2)',mean(A_SpiS3S2)',mean(A_SpiS1S3)'];
    Asd = [std(A_IntS1S2)',std(A_IntS3S2)',std(A_IntS1S3)',...
        std(A_SpiS1S2)',std(A_SpiS3S2)',std(A_SpiS1S3)'];
    A_sh = [A_IntS1S2',A_IntS3S2',A_IntS1S3',A_SpiS1S2',A_SpiS3S2',A_SpiS1S3'];
end
if NplotCases==9
    Am = [mean(A_IntS1S2)',mean(A_IntS3S2)',mean(A_IntS1S3)',...
        mean(A_SpiS1S2)',mean(A_SpiS3S2)',mean(A_SpiS1S3)',...
        mean(A_IntSpiS1)', mean(A_IntSpiS2)', mean(A_IntSpiS3)'];
    Asd = [std(A_IntS1S2)',std(A_IntS3S2)',std(A_IntS1S3)',...
        std(A_SpiS1S2)',std(A_SpiS3S2)',std(A_SpiS1S3)',...
        std(A_IntSpiS1)', std(A_IntSpiS2)', std(A_IntSpiS3)'];
    A_sh = [A_IntS1S2',A_IntS3S2',A_IntS1S3',A_SpiS1S2',A_SpiS3S2',A_SpiS1S3',...
        A_IntSpiS1', A_IntSpiS2', A_IntSpiS3'];
end
if doPlotShuffles
    bsh = bar(A_sh);
    hold on;
    for iCase = 1:NplotCases
        set(bsh(1+(iCase-1)*NShuffles:iCase*NShuffles),'FaceColor',ColorG(iCase,:));
    end
    b = bar(Am,'FaceColor','none');
else
    b = bar(Am,'FaceColor','none');
    for iCase = 1:NplotCases
        set(b(iCase),'FaceColor',ColorG(iCase,:));
    end
    hold on
end

set(b,'BarWidth',1);    % The bars will now touch each other

[numgroups,numbars] = size(Am);%N_mus, N_dir
groupwidth = min(0.8, numbars/(numbars+1.5));
for j = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, Am(:,j), Asd(:,j), 'k', 'linestyle', 'none','CapSize',CapSz);
end
set(gca,'box','off','TickDir','out')
% text(2-3/7,60,'IntS1S2')%,'Color',ColorG(1,:)
% text(2-2/7,55,'IntS3S2')%,'Color',ColorG(2,:)
% text(2-1/7,50,'IntS1S3')%,'Color',ColorG(3,:)
% text(2,60,'SpiS1S2')%,'Color',ColorG(4,:)
% text(2+1/7,55,'SpiS3S2')%,'Color',ColorG(5,:)
% text(2+2/7,50,'SpiS1S3')%,'Color',ColorG(6,:)
set(gca,'YTick',linspace(0,60,4),'XTickLabel',{'Syn#1','Syn#2','Syn#3','Syn#4','Syn#5'}, 'FontSize', 20);
hold off
hLg = legend(CBetweens{1:NplotCases});%,'Location','northwest'
%%
if doPlotShuffles
FigNameAngSh_pdf = strcat('../Figures/AnglesWsh_v',num2str(NplotCases),'.pdf');
f_A = gcf;
exportgraphics(f_A,FigNameAngSh_pdf,'ContentType','vector')

FigNameAngSh_eps = strcat('../Figures/AnglesWsh_v',num2str(NplotCases),'.eps');
exportgraphics(f_A,FigNameAngSh_eps,'ContentType','vector')
else

FigNameAng_pdf = strcat('../Figures/AnglesW_v',num2str(NplotCases),'.pdf');
f_A = gcf;
exportgraphics(f_A,FigNameAng_pdf,'ContentType','vector')

FigNameAng_eps = strcat('../Figures/AnglesW_v',num2str(NplotCases),'.eps');
exportgraphics(f_A,FigNameAng_eps,'ContentType','vector')
end

%%
ColorD = {'b','g','r'};
CaseNames = cell(1,3);
doPlotShuffles = false; %true;
NplotCases =3;
iCond = 2;

Condition = Conditions{iCond};

figure('Name',strcat('Angles between individual Cases and combine W_syn_',Condition,'-Syn-',num2str(N_syn)));
switch Condition
    case 'INTACT'
       A1 = A_IntS1IntCom;
       A2 = A_IntS2IntCom;
       A3 = A_IntS3IntCom;
       CaseNames{1} = 'Intact0.4&IntactCombined';
       CaseNames{2} = 'Intact0.7&IntactCombined';
       CaseNames{3} = 'Intact1.0&IntactCombined';
    case 'SPINAL'
       A1 = A_SpiS1SpiCom;
       A2 = A_SpiS2SpiCom;
       A3 = A_SpiS3SpiCom;
       CaseNames{1} = 'Spinal0.4&SpinalCombined';
       CaseNames{2} = 'Spinal0.7&SpinalCombined';
       CaseNames{3} = 'Spinal1.0&SpinalCombined';
    otherwise
       A1 = A_IntSpiS1;
       A2 = A_IntSpiS2;
       A3 = A_IntSpiS3;
       CaseNames{1} = 'Intact0.4&Spinal0.4';
       CaseNames{2} = 'Intact0.7&Spinal0.7';
       CaseNames{3} = 'Intact1.0&Spinal1.0';
end
Am = [mean(A1)',mean(A2)',mean(A3)'];
Asd = [std(A1)',std(A2)',std(A3)'];
% bsh = bar([A1',A2',A3']);
% hold on;
% for iCase =1:NplotCases
%     set(bsh(1+(iCase-1)*NShuffles:iCase*NShuffles),'FaceColor',ColorD{iCase});
% end
%b= bar(Am,'FaceColor','none');
if doPlotShuffles
    bsh = bar([A1',A2',A3']);
    hold on;
    for iCase = 1:NplotCases
        set(bsh(1+(iCase-1)*NShuffles:iCase*NShuffles),'FaceColor',ColorD{iCase});
    end
    b = bar(Am,'FaceColor','none');
else
    b = bar(Am,'FaceColor','none');
    for iCase = 1:NplotCases
        set(b(iCase),'FaceColor',ColorD{iCase});
    end
    hold on
end
set(b,'BarWidth',1);    % The bars will now touch each other

[numgroups,numbars] = size(Am);%N_mus, N_dir
groupwidth = min(0.8, numbars/(numbars+1.5));
for j = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, Am(:,j), Asd(:,j), 'k', 'linestyle', 'none','CapSize',CapSz);
end
set(gca,'box','off','TickDir','out')
if doPlotShuffles
    text(1-0.25,100,CaseNames(1),'Color',ColorD{1})
    text(1,97,CaseNames(2),'Color',ColorD{2})
    text(1+0.25,94,CaseNames(3),'Color',ColorD{3})
else
    hLg2 = legend(CaseNames);
end

%bar([A_LU;A_UD;A_LD])
set(gca,'YTick',linspace(0,60,4),'XTickLabel',{'Syn#1','Syn#2','Syn#3','Syn#4','Syn#5'}, 'FontSize', 20);
%set(gca,'YTick',linspace(0,60,4),'XTickLabel',{'L&U','U&D','L&D'}, 'FontSize', 20);
hold off
%%
if doPlotShuffles
    FigNameAngSh_pdf = strcat('../Figures/',Condition,'&Combined_AnglesWsh_v1.pdf');
    f_A = gcf;
    exportgraphics(f_A,FigNameAngSh_pdf,'ContentType','vector')

    FigNameAngSh_eps = strcat('../Figures/',Condition,'&Combined_AnglesWsh_v1.eps');
    exportgraphics(f_A,FigNameAngSh_eps,'ContentType','vector')
else

    FigNameAng_pdf = strcat('../Figures/',Condition,'&Combined_AnglesW_v1.pdf');
    f_A = gcf;
    exportgraphics(f_A,FigNameAng_pdf,'ContentType','vector')

    FigNameAng_eps = strcat('../Figures/',Condition,'&Combined_AnglesW_v1.eps');
    exportgraphics(f_A,FigNameAng_eps,'ContentType','vector')
end


