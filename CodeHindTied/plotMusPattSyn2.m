function [h,Sw2Plots] = plotMusPattSyn2(CleanMuscles,NameStr,pos,MusMapCond,...
    emg_set_tm_m,emg_set_tm_sd,Sw2Plot,iCase,NTime,MusclesNew,R2,Rmse,...
    F_tm_m,F_tm_sd,SwSD2Plot,N_syn,H_mD,H_sdD,H_all3D_s,doPlotEach,W_mD,W_sdD,...
    N_Allmuscles,Sw2Plots,Ncycles,doSavePlots)
%  if doPlot
mus_axis = [0 NTime 0 1];
syn_axis = [0 NTime 0 1.2];
if ~contains(NameStr,'All')
    Sw2Plots(iCase) = Sw2Plot;
end

NTime1Cycle = NTime/Ncycles;


Muscles = CleanMuscles;
h =figure('Name',NameStr);
set(h,'position',pos)

mus_plot_rows = length(CleanMuscles);
mus_plot_cols = 1;
NCol = 4;
%        max_m = max(emg_set_tm_m);
max_m_p_sd = zeros(1,length(CleanMuscles));
for i = 1:length(CleanMuscles)
    % m_id = ismember(CondMuscles,CleanMuscles(MusMap(i)));% MusMap(i); %use map if need to sort
    m_id = MusMapCond(i);
    if ~isnan(MusMapCond(i))

        subplot(mus_plot_rows,NCol,NCol*i)
        m_p_sd = emg_set_tm_m(:,m_id)+emg_set_tm_sd(:,m_id);
        m_m_sd = emg_set_tm_m(:,m_id)-emg_set_tm_sd(:,m_id);
        max_m_p_sd(m_id) = max(m_p_sd);
        %collect emgs for different conditions
        %emg_set_tm_mN(:,i,iCond) = emg_set_tm_m(:,m_id)/max_m_p_sd(m_id);
        plot(emg_set_tm_m(:,m_id)/max_m_p_sd(m_id),'k','LineWidth',3); % muscles
        hold on;
        plot(m_p_sd/max_m_p_sd(m_id),'Color', [0.5 0.5 0.5],'LineWidth',2);
        plot(m_m_sd/max_m_p_sd(m_id),'Color', [0.5 0.5 0.5],'LineWidth',2);
        plot(F_tm_m(:,m_id)/max_m_p_sd(m_id),'k.','LineWidth',2), hold off; % restored from Syn and weights

        %     line([Sw2Plot Sw2Plot],[0 max(mus_p_std(:,iMus))],'Color','k','LineStyle','-','LineWidth',2)
        %     line([Sw2Plot+SwSD2Plot Sw2Plot+SwSD2Plot],[0 max(mus_p_std(:,iMus))],'Color','k','LineStyle','--','LineWidth',1)
        %     line([Sw2Plot-SwSD2Plot Sw2Plot-SwSD2Plot],[0 max(mus_p_std(:,iMus))],'Color','k','LineStyle','--','LineWidth',1)
        if contains(NameStr,'All')
            for iC = 1:iCase-1
                Sw2Plot = Sw2Plots(iC) + (NTime/(iCase-1))*(iC-1);
            line([Sw2Plot Sw2Plot],[0 1],'Color','k','LineStyle','-','LineWidth',2)
            line([Sw2Plot+SwSD2Plot Sw2Plot+SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
            line([Sw2Plot-SwSD2Plot Sw2Plot-SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
            

    for iCyc =1:Ncycles
            line([iCyc*NTime1Cycle iCyc*NTime1Cycle],[0 1],'Color','k','LineStyle','-','LineWidth',1)
    end
            % line([NTime1Cycle NTime1Cycle],[0 1],'Color','k','LineStyle','-','LineWidth',1)
            % line([2*NTime1Cycle 2*NTime1Cycle],[0 1],'Color','k','LineStyle','-','LineWidth',1)

                
            end
        else
            line([Sw2Plot Sw2Plot],[0 1],'Color','k','LineStyle','-','LineWidth',2)
            line([Sw2Plot+SwSD2Plot Sw2Plot+SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
            line([Sw2Plot-SwSD2Plot Sw2Plot-SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
        end
        %line([SwSD2Plot, SwSD2Plot],[0 1],'Color','k','LineStyle','--');
        axis(mus_axis);
        %            yticks([0, max_m_p_sd(m_id)/2, max_m_p_sd(m_id)])
        %yticks([0, 1])
        yticks([0,0.5, 1])
        yticklabels({'0',num2str(max_m_p_sd(m_id)/2,1)})
        set(gca,'XColor','none');
        %'XAxisLocation','top',
        set(gca,'YAxisLocation','right','Color','none');

        %            max_m_p_sd(m_id)
        % set(gca,'XColor','none','YColor','none','Visible','off');

        %        text(0,1,strcat(Muscles(m_id),{'    '},'R^2 = ',num2str(R2(m_id),2),'; RMSE = ', num2str(Rmse(m_id),2)))
        text(0,1,strcat(MusclesNew(i),{'    '},'R^2 = ',num2str(R2(m_id),2),'; RMSE = ', num2str(Rmse(m_id),2)))
        hold off;
    end

end

% plot synergies
for i = 1:N_syn

    subplot(N_syn,NCol,NCol*i-1)% for H
    plot(H_mD(1:NTime,i,iCase),'k','LineWidth',2); %Hs_m(:,i)
    hold on;
    %     plot(Hs_m(:,i)+Hs_sd(:,i),'LineWidth',1,'Color','k');% [0.5 0.5 0.5] -grey
    %     plot(Hs_m(:,i)-Hs_sd(:,i),'LineWidth',1,'Color','k');
    plot(H_mD(1:NTime,i,iCase) + H_sdD(1:NTime,i,iCase),'LineWidth',1,'Color','k');% [0.5 0.5 0.5] -grey
    plot(H_mD(1:NTime,i,iCase) - H_sdD(1:NTime,i,iCase),'LineWidth',1,'Color','k');
    if doPlotEach ~= 0
        plot(squeeze(H_all3D_s(:,i,1:doPlotEach:end)),Color=[0.5 0.5 0.5])
    end
        if contains(NameStr,'All')
            for iC = 1:iCase-1
                Sw2Plot = Sw2Plots(iC) + (NTime/(iCase-1))*(iC-1);
            line([Sw2Plot Sw2Plot],[0 1],'Color','k','LineStyle','-','LineWidth',2)
            line([Sw2Plot+SwSD2Plot Sw2Plot+SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
            line([Sw2Plot-SwSD2Plot Sw2Plot-SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)

    for iCyc =1:Ncycles
            line([iCyc*NTime1Cycle iCyc*NTime1Cycle],[0 1],'Color','k','LineStyle','-','LineWidth',1)
    end
            % line([NTime/3 NTime/3],[0 1],'Color','k','LineStyle','-','LineWidth',1)
            % line([2*NTime/3 2*NTime/3],[0 1],'Color','k','LineStyle','-','LineWidth',1)
                
            end
        else
            line([Sw2Plot Sw2Plot],[0 1],'Color','k','LineStyle','-','LineWidth',2)
            line([Sw2Plot+SwSD2Plot Sw2Plot+SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
            line([Sw2Plot-SwSD2Plot Sw2Plot-SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
        end
        % if NTime>200
        %     line([NTime/3 NTime/3],[0 1],'Color','k','LineStyle','-','LineWidth',1)
        %     line([2*NTime/3 2*NTime/3],[0 1],'Color','k','LineStyle','-','LineWidth',1)
        % end
    % line([Sw2Plot Sw2Plot],[0 1],'Color','k','LineStyle','-','LineWidth',2)
    % line([Sw2Plot+SwSD2Plot Sw2Plot+SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
    % line([Sw2Plot-SwSD2Plot Sw2Plot-SwSD2Plot],[0 1],'Color','k','LineStyle','--','LineWidth',1)
    %text(0,1,strcat(['H syn# = ' num2str(i)]))
    ylabel(strcat(['C_' num2str(i)]))

    hChildren = get(gca,'Children');%reverse order of plotting (first plot will be in front)
    set(gca,'Children',hChildren(end:-1:1))

    hold off;
    axis(syn_axis);


    %title(strcat(['H syn# = ' num2str(i)]))
    subplot(N_syn,NCol,[NCol*(i-1)+1,NCol*(i-1)+2])% for W
    %            b= bar(W_sm_m(i,MusMap));set(b,'FaceColor','k','EdgeColor',[1 1 1]);
    b= bar(W_mD(i,:,iCase));set(b,'FaceColor','k','EdgeColor',[1 1 1]);
    hold on;
    %            errorbar(W_sm_m(i,MusMap), W_sm_sd(i,MusMap), 'k', 'linestyle', 'none');
    errorbar(W_mD(i,:,iCase), W_sdD(i,:,iCase), 'k', 'linestyle', 'none');
    hold off;
    ylabel(strcat(['W_' num2str(i)]))
    %title(strcat(['W syn# = ' num2str(i)]))
    axis([0 N_Allmuscles+1 0 1.2]);
    %            set(gca,'XTick',linspace(1,N_Allmuscles,N_Allmuscles),'XTickLabel',Muscles(MusMap));
    set(gca,'XTick',linspace(1,N_Allmuscles,N_Allmuscles),'XTickLabel',MusclesNew);
end
hold off;
if doSavePlots
    NShuffles = size(H_all3D_s,3);
    combFigName = strcat('../Figures/Fig_',NameStr,'_Selected_',num2str(NShuffles/doPlotEach),'.pdf');
    print(combFigName,'-dpdf','-bestfit','-r0','-painters')%);
end
