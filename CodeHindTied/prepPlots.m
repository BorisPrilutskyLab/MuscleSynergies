function   [emg_set_tm_m,emg_set_tm_sd,F_all,F_alls,F_tm_m,F_tm_sd,F_tms_m,...
    F_tms_sd,W_sm_m,W_sm_sd,W_sms_m,W_sms_sd,Hs_m,Hs_sd,Hss_m,Hss_sd,...
    H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs]...
    = prepPlots(V_all,NTime,NCondMuscles, NShuffles,N_Allmuscles,MusMapCond,...
    W_all3D_s,H_all3D_s,W_all3D_ss,H_all3D_ss,SideSpeed,ContraSideSpeed,iCond,...
    H_mD,H_sdD,H_mDs,H_sdDs,W_mD,W_sdD,W_mDs,W_sdDs)%,Hmaps


% prep to plot
emg_set_time_mus_rest = permute(V_all,[3,2,1]);%reshape(emg_set_all,NTime,N_muscles,[]);
emg_set_tm_m = mean(emg_set_time_mus_rest,3);
emg_set_tm_sd = std(emg_set_time_mus_rest,[],3);

F_all = zeros(NTime,NCondMuscles, NShuffles);
F_alls = zeros(NTime,NCondMuscles, NShuffles);


for i_r =1: NShuffles
    W = squeeze(W_all3D_s(:,:,i_r))';%11 *Nsyn
    H = squeeze(H_all3D_s(:,:,i_r))';%Nsyn * 100
    F = W*H;                      %11 *100
    F_all(:,:,i_r) = F';%100x11x2000

    Ws = squeeze(W_all3D_ss(:,:,i_r))';%11 *Nsyn
    Hs = squeeze(H_all3D_ss(:,:,i_r))';%Nsyn * 100
    Fs = Ws*Hs;                      %11 *100
    F_alls(:,:,i_r) = Fs';%100x11x2000

end

F_tm = reshape(F_all,NTime,NCondMuscles,[]);
F_tm_m = mean(F_tm,3);
F_tm_sd = std(F_tm,[],3);

F_tms = reshape(F_alls,NTime,NCondMuscles,[]);
F_tms_m = mean(F_tms,3);
F_tms_sd = std(F_tms,[],3);

W_sm_m = mean(W_all3D_s,3);
W_sm_sd = std(W_all3D_s,[],3);

W_sms_m = mean(W_all3D_ss,3);
W_sms_sd = std(W_all3D_ss,[],3);

Hs_m = squeeze(mean(H_all3D_s,3));
Hs_sd = squeeze(std(H_all3D_s,[],3));

Hss_m = squeeze(mean(H_all3D_ss,3));
Hss_sd = squeeze(std(H_all3D_ss,[],3));

for iMus = 1:N_Allmuscles
    if ~isnan(MusMapCond(iMus))

        W_mD(:,iMus,iCond) = W_sm_m(:,MusMapCond(iMus));
        W_sdD(:,iMus,iCond) = W_sm_sd(:,MusMapCond(iMus));

        % W_mD(:,iMus,iCond) = W_sm_m(Hmaps(iCond,:),MusMapCond(iMus));
        % W_sdD(:,iMus,iCond) = W_sm_sd(Hmaps(iCond,:),MusMapCond(iMus));
        % 
        % 
        % % if strcmp(SideSpeed,'1.0')
        % %     W_mD(:,iMus,iCond) = circshift(W_sm_m(:,MusMapCond(iMus)),1,1);
        % %     W_sdD(:,iMus,iCond) = circshift(W_sm_sd(:,MusMapCond(iMus)),1,1);
        % % elseif strcmp(ContraSideSpeed,'1.0')
        % %     % switch 1 and 5 synergies
        % %     W_mD(:,iMus,iCond) = W_sm_m([5,2,3,4,1],MusMapCond(iMus));
        % %     W_sdD(:,iMus,iCond) = W_sm_sd([5,2,3,4,1],MusMapCond(iMus));
        % %
        % % else
        % %     W_mD(:,iMus,iCond) = W_sm_m(:,MusMapCond(iMus));
        % %     W_sdD(:,iMus,iCond) = W_sm_sd(:,MusMapCond(iMus));
        % % end


        W_mDs(:,iMus,iCond) = W_sms_m(:,MusMapCond(iMus));
        W_sdDs(:,iMus,iCond) = W_sms_sd(:,MusMapCond(iMus));

        %         W_4D(:,iMus,:) = W_all3D_s(:,MusMapCond(iMus),:);
        %         W_4Ds(:,iMus,:) = W_all3D_ss(:,MusMapCond(iMus),:);

    end
end

H_mD(1:NTime,:,iCond) = Hs_m;
H_sdD(1:NTime,:,iCond) = Hs_sd;

% H_mD(1:NTime,:,iCond) = Hs_m(:,Hmaps(iCond,:));
% H_sdD(1:NTime,:,iCond) = Hs_sd(:,Hmaps(iCond,:));
% 
% % if strcmp(SideSpeed,'1.0')
% %     H_mD(1:NTime,:,iCond) = circshift(Hs_m,1,2);
% %     H_sdD(1:NTime,:,iCond) = circshift(Hs_sd,1,2);
% % elseif strcmp(ContraSideSpeed,'1.0')
% %     % switch 1 and 5 synergies
% %     H_mD(1:NTime,:,iCond) = Hs_m(:,[5,2,3,4,1]);
% %     H_sdD(1:NTime,:,iCond) = Hs_sd(:,[5,2,3,4,1]);
% % else
% %     H_mD(1:NTime,:,iCond) = Hs_m;
% %     H_sdD(1:NTime,:,iCond) = Hs_sd;
% % end

H_mDs(1:NTime,:,iCond) = Hss_m;
H_sdDs(1:NTime,:,iCond) = Hss_sd;

% H_4D = H_all3D_s;
% H_4Ds = H_all3D_ss;

% plot muscles and restored from synergies and weight patterns
%figure(1)%,set(figure(1), 'WindowStyle', 'docked');
