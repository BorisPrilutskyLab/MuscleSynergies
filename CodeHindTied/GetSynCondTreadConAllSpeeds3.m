clear variables
%DataFolder = 'D:\Dropbox (GaTech)\ForeLimbEMGsFrigon\Data\';
DataFolder = 'E:\DropBox(Gatech)\Dropbox (GaTech)\ForeLimbEMGsFrigon\Data\';

N_syn = 5;
Nreplicates = 30;
NTime= 101;
NShuffles = 1000;

TreadmillConditions = {'SPLIT','TIED'};
TreadmillCondition = 'TIED';
Conditions = {'INTACT','SPINAL'};
SideSpeeds = {'0.4','0.7','1.0'};
NCondMuscles = 16;
for iCond = 1:length(Conditions)
    Condition = Conditions{iCond};%'SPINAL';%'INTACT';%


    V_allSpeeds = [];
    TCondSpeeds = [];
    for iSideSpeed = 1:length(SideSpeeds)
        % iSideSpeed = 1;
        SideSpeed = SideSpeeds{iSideSpeed};
        ContraSideSpeed = SideSpeed; %for tied

        f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',num2str(N_syn),'_SideSpeed-',SideSpeed,...
            '_ContraSideSpeed-',ContraSideSpeed,'_NShuffles-',...
            num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
        load(f_name)
        V_allSpeeds = cat(3,V_allSpeeds,V_all);%V_all(Nshuffles,NMus,NTime)
        TCondSpeeds = [TCondSpeeds;TCond];
    end

    %
    emg_sets = permute(V_allSpeeds,[2,3,1]);%to be consistent with previous code

    W_all = zeros(NShuffles,NCondMuscles,N_syn);
    H_all = zeros(NShuffles,N_syn,NTime*length(SideSpeeds));
    V_all = zeros(NShuffles,NCondMuscles,NTime*length(SideSpeeds));
    F_all = zeros(NShuffles,NCondMuscles,NTime*length(SideSpeeds));
    R2 = zeros(NShuffles,1);
    Rmse = zeros(NShuffles,1);
    R2EachMus = zeros(NShuffles,NCondMuscles);
    RmseEachMus = zeros(NShuffles,NCondMuscles);

    tic
    for j_s =1:NShuffles %cycles
        V = squeeze(emg_sets(:,:,j_s));%11x100
        stdev = std(V,0,2);
        V0 = diag(1./stdev)*V;%scale the data to have unit variance of this data set

        [W,H] = nnmf(V0,N_syn,'replicates',Nreplicates,'algorithm','mult');%
        W = diag(stdev)*W;
        max_i=max(W);% vector with max activation values
        for i=1:N_syn
            H(i,:)=H(i,:)*max_i(i);
            W(:,i)=W(:,i)/max_i(i);
        end
        W_all(j_s,:,:) = W; %N_muscles,N_syn
        H_all(j_s,:,:) = H; %N_syn,NTime
        F = W*H;
        F_all(j_s,:,:) = F;%N_muscles,NTime
        V1 = reshape(V',1,[]);
        F1 = reshape(F',1,[]);
        [r2, rmse] = rsquare(V1,F1);
        R2(j_s,1) = r2;
        Rmse(j_s,1)= rmse;
        for i_mus = 1:NCondMuscles
            [r2, rmse] = rsquare(V(i_mus,:),F(i_mus,:));
            R2EachMus(j_s,i_mus) = r2;
            RmseEachMus(j_s,i_mus)= rmse;
        end
        V_all(j_s,:,:) = V;
    end
    toc

    R2_in_1_line = reshape(R2,1,[]);
    disp(strcat(' mean R2 =',num2str(mean(R2(:)))))
    disp(strcat(' min R2 =',num2str(min(R2(:)))))
    disp(strcat(' SD R2 =',num2str(std(R2(:)))))
    %[mean(R2_in_1_line),min(R2_in_1_line),std(R2_in_1_line)]



    f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',...
        num2str(N_syn),'_SideSpeedAll','_NShuffles-',num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');
    %
    CondMuscles = unique(TCondSpeeds.Muscle);
    save(f_name, 'R2','Rmse','V_all','W_all','H_all','CondMuscles','TCondSpeeds');
    save(f_name, 'F_all','R2EachMus','RmseEachMus','-append');

end
