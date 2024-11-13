clear variables
close all
DataFolder = '..\DataHindTied\';

Conditions = {'INTACT','SPINAL'};
TreadmillCondition = 'TIED';
Speeds = {'0.4','0.7','1.0'};

[TClean,CleanMuscles,NCleanMuscles,N_Allmuscles,SideSpeeds,NSideSpeeds,...
    MusMap,MusclesNew] = loadAllSpeedsAllCondTable2(DataFolder,TreadmillCondition);
%MusNames_str = char(join(CleanMuscles,','));
%
Cats = unique(TClean.Cat);
ANCyclesCatCondSpeed = zeros(length(Cats),length(Conditions),length(Speeds)); 
for iCat = 1:length(Cats)
    Cat = Cats(iCat);
    for iCond = 1:length(Conditions)
        Condition = Conditions(iCond);
        for iSpeed = 1:length(Speeds)
            Speed = Speeds(iSpeed);
            T1 = TClean(ismember(TClean.Cat,Cat)&ismember(TClean.Condition,Condition)&...
                ismember(TClean.SideSpeed,Speed),:);
            Connectors = unique(T1.Conn);
            NCyclesCatCondSpeed = 0; 
            for iConn = 1:length(Connectors)
                Connector = Connectors(iConn);
                T2 = T1(ismember(T1.Conn,Connector),:);
                Sides = unique(T2.MSide);
                for iSide = 1:length(Sides)
                    Side = Sides(iSide);
                    T3 = T2(ismember(T2.MSide,Side),:);
                    Cycles = unique(T3.CycleN);
                    NCyclesCatCondSpeed = NCyclesCatCondSpeed +length(Cycles);
                end
            end
            ANCyclesCatCondSpeed(iCat,iCond,iSpeed) = NCyclesCatCondSpeed;
        end
    end
end

ANCyclesCatCondSpeed2 = reshape(ANCyclesCatCondSpeed,[length(Cats),length(Conditions)*length(Speeds)]);
ANCyclesCatCondSpeed1 = permute(ANCyclesCatCondSpeed,[1,3,2]);
ANCyclesCatCondSpeed3 = reshape(ANCyclesCatCondSpeed1,[length(Cats),length(Conditions)*length(Speeds)]);
%
TVars = {'Animal','Sex','Mass(kg)'};
for iCond = 1:length(Conditions)
    Condition = Conditions(iCond);
    for iSpeed = 1:length(Speeds)
        Speed = Speeds(iSpeed);
        TVars = [TVars,strcat(Condition,'-',Speed,'(m/s)')];

    end
end
%
T4 = cell2table([Cats,cell(9,2),num2cell(ANCyclesCatCondSpeed3)],'VariableNames',TVars);
writetable(T4,'..\Data\TCycles.txt')
%
TVars2Keep = {'Cat','Conn','MSide','CycleN','CycleTime','StanceTime','DF','SideSpeed','Condition'};
TKin_red = TClean(:,TVars2Keep);
TKin = unique(TKin_red);
SwingTime = TKin.CycleTime - TKin.StanceTime;

Tkin2 = addvars(TKin,SwingTime,'After','StanceTime');

writetable(Tkin2,'..\Data\TKin.txt')
