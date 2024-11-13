function    [TCase,CondMuscles,MusMapCond,NCondAllMuscles,H_all,W_all,V_all,R2,Rmse]...
    = loadOneSpeedsTable5(DataFolder,Condition,TreadmillCondition,N_syn,...
    SideSpeed,ContraSideSpeed,NShuffles,N_Allmuscles,TClean,CleanMuscles,MusclesNew,NCondMuscles,CondName)

if contains(SideSpeed,'AllSpeeds')
        if contains(CondName,'AllCond')
        % f_name = strcat(DataFolder,'CondAllSpeedsAll_',TreadmillCondition,'_Syn-',num2str(N_syn),...
        %     '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
        f_name = strcat(DataFolder,TreadmillCondition,'CondAllSpeedAll','_Syn-',num2str(N_syn),...
            '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
        else
        f_name = strcat(DataFolder,TreadmillCondition,'Cond',Condition,'SpeedAll','_Syn-',num2str(N_syn),...
            '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
        end

    % f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',num2str(N_syn),...
    %     '_SideSpeedAll_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
    % if contains(CondName,'AllSpeedsCond')
    %     f_name = strcat(DataFolder,TreadmillCondition,'Cond',Condition,'SpeedAll','_Syn-',num2str(N_syn),...
    %         '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
    % elseif contains(CondName,'AllCond')
    %     f_name = strcat(DataFolder,'CondAllSpeedsAll_',TreadmillCondition,'_Syn-',num2str(N_syn),...
    %         '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
    % end
elseif contains(CondName,'AllCond')
    f_name = strcat(DataFolder,TreadmillCondition,'CondAllSpeed',SideSpeed,'_Syn-',num2str(N_syn),...
        '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');

    % f_name = strcat(DataFolder,TreadmillCondition,'Cond',Condition,'Speed',SideSpeed,'_Syn-',num2str(N_syn),...
    %     '_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');

else %separate condtions/speeds
    f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',num2str(N_syn),'_SideSpeed-',SideSpeed,...
        '_ContraSideSpeed-',ContraSideSpeed,'_NShuffles-',...
        num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
end

load(f_name,'W_all','H_all','V_all','R2','Rmse');
%[~,NCondMuscles,~] = size(W_all);%NShuffles,N_muscles,N_syn
%[~,~,NTime] = size(H_all);%NShuffles,N_syn,NTime

if contains(SideSpeed,'AllSpeeds')
    if contains(CondName,'AllCond')
        TCase = TClean;
    else
        TCase = TClean(ismember(TClean.Condition,Condition),:);
    end
else
    if contains(CondName,'AllCond')
        TCase = TClean(ismember(TClean.SideSpeed,SideSpeed)&...
            ismember(TClean.ContraSideSpeed,ContraSideSpeed),:);
    else
        TCase = TClean(ismember(TClean.SideSpeed,SideSpeed)&...
            ismember(TClean.ContraSideSpeed,ContraSideSpeed)&...
            ismember(TClean.Condition,Condition),:);
    end
end
%number of samples for each muscle from 'CleanMuscles'
NCondAllMuscles = cellfun(@(x) sum(ismember(TCase.Muscle,x)), CleanMuscles);

%muscles we have for this condition
CondMuscles = CleanMuscles(NCondAllMuscles~=0);

MusMapCond = nan(1,N_Allmuscles);
for i = 1:N_Allmuscles
    indx = find(ismember(CondMuscles,MusclesNew(i)));
    if ~isempty(indx)
        MusMapCond(i) = indx;
    end
end

% MusNames_str = char(join(CleanMuscles,','));