function  [TW,TH,TWs,THs,TWS1D,THS1D,TWS1Ds,THS1Ds] = ...
    fillStatTables5(TW_vars,TH_vars,TW,TH,TWs,THs,NShuffles,N_syn,...
    NCondMuscles,CondMuscles,W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,...
    Condition,SideSpeed,ContraSideSpeed,NBins,NTime,SideSpeeds,Conditions)
% fill out tables for statistics
%TW_vars = {'Syn','Muscle','Direction','W'};
%TH_vars = {'Syn','tBin','Direction','H'};
Nsamples=  NShuffles;
% ConditionV = repmat({Condition},Nsamples,1);
%SideSpeeds = {'0.4','0.7','1.0'}; %tmp
for iSyn = 1:N_syn
    TWS1D = cell2table(cell(0,length(TW_vars)),'VariableNames',TW_vars);
    THS1D = cell2table(cell(0,length(TH_vars)),'VariableNames',TH_vars);

    TWS1Ds = cell2table(cell(0,length(TW_vars)),'VariableNames',TW_vars);
    THS1Ds = cell2table(cell(0,length(TH_vars)),'VariableNames',TH_vars);

    % W
    ConditionV = repmat({Condition},Nsamples,1);
    SideSpeedV =repmat({SideSpeed},Nsamples,1);
    ContraSideSpeedV =  repmat({ContraSideSpeed},Nsamples,1);

    for iMus = 1:NCondMuscles
        Muscle = CondMuscles(iMus);
        W =num2cell(squeeze(W_all3D_s(iSyn,iMus,:)));
        Syns = repmat({iSyn},Nsamples,1);

        TWM = table(ConditionV,SideSpeedV,ContraSideSpeedV,Syns,repmat(Muscle,Nsamples,1),W,'VariableNames',TW_vars);
        TW = [TW;TWM];
        TWS1D = [TWS1D;TWM];
        %------------spatial sorting------------------
        Ws =num2cell(squeeze(W_all3D_ss(iSyn,iMus,:)));
        TWMs = table(ConditionV,SideSpeedV,ContraSideSpeedV,Syns,repmat(Muscle,Nsamples,1),Ws,'VariableNames',TW_vars);
        TWs = [TWs;TWMs];
        TWS1Ds = [TWS1Ds;TWMs];
        %-
    end
    for iBin = 1:NBins
        if contains(SideSpeed,'AllSpeeds')
            if contains(Condition,'AllCond')
                %--'AllSpeedsAllCond'
                for iCond = 1:length(Conditions)
                    Cond = strcat(Condition,Conditions{iCond});
                    ConditionV = repmat({Cond},Nsamples,1);
                    for iSpeed = 1:length(SideSpeeds)

                        Speed = strcat(SideSpeed,SideSpeeds{iSpeed});
                        ContraSpeed = strcat(SideSpeed,SideSpeeds{iSpeed});

                        SideSpeedV =repmat({Speed},Nsamples,1);
                        ContraSideSpeedV =  repmat({ContraSpeed},Nsamples,1);

                        rangeBin = (1+(iBin-1)*NBins:NTime/NBins+(iBin-1)*NBins) + NTime*(iSpeed-1) + length(SideSpeeds)*NTime*(iCond-1);
                        % rangeBin = (1+(iBin-1)*NBins:NTime/NBins+(iBin-1)*NBins) + NTime*(iSpeed-1);
                        % rangeBin = (1+(iBin-1)*NBins:NTime/NBins+(iBin-1)*NBins) + NTime*(iCond-1);

                        HBin = num2cell(squeeze(mean(H_all3D_s(rangeBin,iSyn,:),1)));
                        THM = table(ConditionV,SideSpeedV,ContraSideSpeedV,Syns,repmat({iBin},Nsamples,1),HBin,'VariableNames',TH_vars);
                        TH = [TH;THM];
                        THS1D = [THS1D;THM];
                    end
                end
                %--'AllSpeedsAllCond'
            else
                %--'AllSpeeds, different cond'
                for iSpeed = 1:length(SideSpeeds)
                    Speed = strcat(SideSpeed,SideSpeeds{iSpeed});
                    ContraSpeed = strcat(SideSpeed,SideSpeeds{iSpeed});

                    ConditionV = repmat({Condition},Nsamples,1);
                    SideSpeedV =repmat({Speed},Nsamples,1);
                    ContraSideSpeedV =  repmat({ContraSpeed},Nsamples,1);
                    rangeBin = (1+(iBin-1)*NBins:NTime/NBins+(iBin-1)*NBins) + NTime*(iSpeed-1);
                    HBin = num2cell(squeeze(mean(H_all3D_s(rangeBin,iSyn,:),1)));
                    THM = table(ConditionV,SideSpeedV,ContraSideSpeedV,Syns,repmat({iBin},Nsamples,1),HBin,'VariableNames',TH_vars);
                    TH = [TH;THM];
                    THS1D = [THS1D;THM];
                end
                %--'AllSpeeds, different cond'
            end 
        else
            if contains(Condition,'AllCond')
                % --'AllCond, different speeds '
                for iCond = 1:length(Conditions)
                    Cond = strcat(Condition,Conditions{iCond});
                    ConditionV = repmat({Cond},Nsamples,1);
                    SideSpeedV =repmat({SideSpeed},Nsamples,1);
                    ContraSideSpeedV =  repmat({ContraSideSpeed},Nsamples,1);
                    rangeBin = (1+(iBin-1)*NBins:NTime/NBins+(iBin-1)*NBins) + NTime*(iCond-1);
                    HBin = num2cell(squeeze(mean(H_all3D_s(rangeBin,iSyn,:),1)));
                    THM = table(ConditionV,SideSpeedV,ContraSideSpeedV,Syns,repmat({iBin},Nsamples,1),HBin,'VariableNames',TH_vars);
                    TH = [TH;THM];
                    THS1D = [THS1D;THM];
                end
                % --'AllCond, different speeds'
            else
                %separate Cond/speeds
                ConditionV = repmat({Condition},Nsamples,1);
                SideSpeedV =repmat({SideSpeed},Nsamples,1);
                ContraSideSpeedV =  repmat({ContraSideSpeed},Nsamples,1);
                rangeBin = 1+(iBin-1)*NBins:NTime/NBins+(iBin-1)*NBins;
                HBin = num2cell(squeeze(mean(H_all3D_s(rangeBin,iSyn,:),1)));
                THM = table(ConditionV,SideSpeedV,ContraSideSpeedV,Syns,repmat({iBin},Nsamples,1),HBin,'VariableNames',TH_vars);
                TH = [TH;THM];
                THS1D = [THS1D;THM];
                %separate Cond/speeds
            end
        end

        %------------spatial sorting------------------
        HBins = num2cell(squeeze(mean(H_all3D_ss(rangeBin,iSyn,:),1)));
        THMs = table(ConditionV,SideSpeedV,ContraSideSpeedV,Syns,repmat({iBin},Nsamples,1),HBins,'VariableNames',TH_vars);
        THs = [THs;THMs];
        THS1Ds = [THS1Ds;THMs];
        %-
    end
end