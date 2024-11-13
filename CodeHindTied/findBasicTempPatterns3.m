function [Hbase, Hbase_set,MV,MI,Hbase_setSorted,R2,Rmse] = ...
    findBasicTempPatterns3(H_all,NTime,N_syn,NHbase_repeat,NameStr)
        %  finds basic Temporal patterns by using nnmf on already found patterns
        Ndims = ndims(H_all);
    switch Ndims
        case 4
            H_allPerm = permute(H_all,[4,3,2,1]);
        case 3
            H_allPerm = permute(H_all,[3,2,1]);
        case 2
            H_allPerm = permute(H_all,[2,1]);
        otherwise
            H_allPerm = H_all;
    end
        H_all_unsorted = reshape(H_allPerm,NTime,[]); %make rectangular matrix with all patterns
        [~,NSamples] = size(H_all_unsorted);
W_all = nan(NHbase_repeat,NTime,N_syn);
H_all = nan(NHbase_repeat,N_syn,NSamples);
F_all = nan(NHbase_repeat,NTime,NSamples);
R2 = nan(NHbase_repeat,1);
Rmse = nan(NHbase_repeat,1);
Nreplicates = 30;
        tic
        for ir = 1:NHbase_repeat
            V = H_all_unsorted;
    stdev = std(V,0,2);
    V0 = diag(1./stdev)*V;%scale the data to have unit variance of this data set

    [W,H] = nnmf(V0,N_syn,'replicates',Nreplicates,'algorithm','mult');%
    W = diag(stdev)*W;
    max_i=max(W);% vector with max activation values
    for i=1:N_syn
        H(i,:)=H(i,:)*max_i(i);
        W(:,i)=W(:,i)/max_i(i);
    end
    W_all(ir,:,:) = W; %N_muscles,N_syn
    H_all(ir,:,:) = H; %N_syn,NTime
    F = W*H;
    F_all(ir,:,:) = F;%N_muscles,NTime
    V1 = reshape(V',1,[]);
    F1 = reshape(F',1,[]);
    [r2, rmse] = rsquare(V1,F1);
    R2(ir,1) = r2;
    Rmse(ir,1)= rmse;
        end
        toc
        Hbase_set = permute(W_all,[2,3,1]);
        disp(strcat('temp basic pattern found'))

        % find if temp basic patterns split in several sets
        NNonSplitSets = nan(1,NHbase_repeat);
        Hbase_setSorted = nan(size(Hbase_set));
        MV_all = nan(NHbase_repeat,N_syn,NHbase_repeat);
        MI_all= nan(NHbase_repeat,N_syn,NHbase_repeat);
        for ir2 = 1:NHbase_repeat
            [Hbase_setSorted(:,:,ir2),MV_all(:,:,ir2),MI_all(:,:,ir2),NUnique] = CorrHbaseSets(ir2,Hbase_set,NHbase_repeat,N_syn);
            NNonSplitSets(ir2) = sum(NUnique==N_syn);
        end
        
        [maxCorrV,maxCorrID] = max(R2);
        [minRmseV,minRmseID] = min(Rmse);

        [mSet,iSet] =max(NNonSplitSets);

        disp(strcat('maxCorrID= ',num2str(maxCorrID),', minRmseID = ',num2str(minRmseID),...
            ', nonSplitSetsAtmaxCorr = ',...
            num2str(NNonSplitSets(maxCorrID)),', maxnonSplitSets = ',num2str(mSet)));
        
        

        Hbase = Hbase_setSorted(:,:,maxCorrID);
        MV = MV_all(:,:,maxCorrID);
        MI = MI_all(:,:,maxCorrID);

        % Hbase_set = zeros(NTime,N_syn,NHbase_repeat);
        % tic
        % for ir = 1:NHbase_repeat
        %     [Hbase,~] = nnmf(H_all_unsorted,N_syn,'replicates',5); %Finds basic patterns
        %     Hbase_set(:,:,ir)= Hbase./max(Hbase);
        % end
        % toc
        % disp(strcat('temp basic pattern found'))

        % find if temp basic patterns split in several sets
       %  NNonSplitSets = nan(1,NHbase_repeat);
       %  Hbase_setSorted = nan(size(Hbase_set));
       %  MV_all = nan(NHbase_repeat,N_syn,NHbase_repeat);
       %  MI_all= nan(NHbase_repeat,N_syn,NHbase_repeat);
       %  for ir2 = 1:NHbase_repeat
       %      [Hbase_setSorted(:,:,ir2),MV_all(:,:,ir2),MI_all(:,:,ir2),NUnique] = CorrHbaseSets(ir2,Hbase_set,NHbase_repeat,N_syn);
       %      NNonSplitSets(ir2) = sum(NUnique==N_syn);
       % end
       %  [mSet,iSet] =max(NNonSplitSets);
       %  Hbase = Hbase_setSorted(:,:,iSet);
       %  MV = MV_all(:,:,iSet);
       %  MI = MI_all(:,:,iSet);

        % Hbase1 = squeeze(Hbase_set(:,:,1));
        % MI = zeros(NHbase_repeat,N_syn);
        % MV = MI;
        % idx1 = []; 
        % idx2 = [];
        % Hbase_setSorted = Hbase_set;
        % for ir = 1:NHbase_repeat
        %     [MV(ir,:), MI(ir,:)] = max(corr(Hbase1,squeeze(Hbase_set(:,:,ir))).^2);
        %     [~,idxu,idxMI] = unique(MI(ir,:));% Unique values
        %     %[C,ia,ic] = unique; C = A(ia),A = C(ic)
        % 
        %     if numel(idxu)<N_syn
        %         idx2 = [idx2,ir];
        %     else
        %         idx1 = [idx1,ir];
        %         Hbase_setSorted(:,MI(ir,:),ir)= Hbase_set(:,:,ir);
        %     end
        % end
        disp(strcat(NameStr,'_MaxR2(mean/SD) for finding temp base '))
        disp(strcat(num2str(mean(MV(2:end,:)))))
        disp(strcat(num2str(std(MV(2:end,:)))))
        disp(num2str(NNonSplitSets))

        disp(strcat('restored H from HbaseR2.-?',num2str(R2')))
% 
%         if isempty(idx2) % set didn't split
%             Hbase = mean(Hbase_setSorted,3);
%             disp('Basic temporal set didn"t splits');
%         else
%             disp('Basic temporal set splits');
%         end
 disp(' ') 