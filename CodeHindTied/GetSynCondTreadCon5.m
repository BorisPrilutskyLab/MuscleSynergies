clear variables
DataFolder = '..\DataHindTied\';

Nreplicates = 20;
NTime= 101;
NShuffles = 50;
Conditions = {'INTACT','SPINAL'};
TreadmillConditions = {'SPLIT','TIED'};

TreadmillCondition = 'TIED';
TableSpeedsFilename = strcat(DataFolder,TreadmillCondition,'_','TablesVicB.mat');
load(TableSpeedsFilename,'TCleanTCond15R')
for N_syn = 1:14
    for iCond = 1:length(Conditions)
        Condition = Conditions{iCond};%'SPINAL';%'INTACT';%

        TClean = TCleanTCond15R(ismember(TCleanTCond15R.Condition,Condition),:);

        SideSpeeds = unique(TClean.SideSpeed);%Speed of treadmill of side of limb

        for iSideSpeed = 1:length(SideSpeeds)
            SideSpeed = SideSpeeds{iSideSpeed};

            % possible speeds of opposite belt for choosen this side speed
            ContraSideSpeeds = unique(TClean(ismember(TClean.SideSpeed,SideSpeed),:).ContraSideSpeed);

            for iContraSideSpeed = 1:length(ContraSideSpeeds)
                ContraSideSpeed = ContraSideSpeeds{iContraSideSpeed};

                % one condition is unique combination of SideSpeed and ContraSideSpeed
                TCond = TClean(ismember(TClean.SideSpeed,SideSpeed)&...
                    ismember(TClean.ContraSideSpeed,ContraSideSpeed),:);

                CondMuscles = unique(TCond.Muscle);
                NCondMuscles = length(CondMuscles);
                emg_sets = nan(NTime,NCondMuscles,NShuffles);
                for iMus = 1:length(CondMuscles)
                    Muscle = CondMuscles{iMus};
                    TMus = TCond(ismember(TCond.Muscle,Muscle),:);

                    EMGMusSet = nan(NTime,NShuffles);
                    for iRep = 1:floor(NShuffles/height(TMus))+1
                        EMGMus = GetEMGOneMusCatsMaximized2(TMus);%v2- normalized across cond&speed

                        i1 = 1 + (iRep-1)*height(TMus);
                        i2 = iRep*height(TMus);
                        EMGMusSet(:,i1:i2) = EMGMus;
                    end
                    emg_sets(:,iMus,:) =  EMGMusSet(:,1:NShuffles);

                end

                W_all = zeros(NShuffles,NCondMuscles,N_syn);
                H_all = zeros(NShuffles,N_syn,NTime);
                V_all = zeros(NShuffles,NCondMuscles,NTime);
                F_all = zeros(NShuffles,NCondMuscles,NTime);
                R2 = zeros(NShuffles,1);
                Rmse = zeros(NShuffles,1);
                R2EachMus = zeros(NShuffles,NCondMuscles);
                RmseEachMus = zeros(NShuffles,NCondMuscles);

                tic
                for j_s =1:NShuffles %cycles
                    V = squeeze(emg_sets(:,:,j_s))';%11x100
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
                disp(strcat(num2str(N_syn),'.',TreadmillCondition,'.',Condition,'.',SideSpeed))
                disp(strcat(' mean R2 =',num2str(mean(R2(:)))))
                disp(strcat(' min R2 =',num2str(min(R2(:)))))
                disp(strcat(' SD R2 =',num2str(std(R2(:)))))

                f_name = strcat(DataFolder,Condition,'_',TreadmillCondition,'_Syn-',num2str(N_syn),'_SideSpeed-',SideSpeed,...
                    '_ContraSideSpeed-',ContraSideSpeed,'_NShuffles-',...
                    num2str(NShuffles),'_NMus-',num2str(NCondMuscles),'.mat');

                save(f_name, 'R2','Rmse','V_all','W_all','H_all','TCond');
                save(f_name, 'F_all','R2EachMus','RmseEachMus','CondMuscles','-append');
            end
        end
    end
end