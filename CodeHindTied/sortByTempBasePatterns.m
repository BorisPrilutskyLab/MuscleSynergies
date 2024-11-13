function  [HbaseS,H_all3D_s,W_all3D_s] = sortByTempBasePatterns(Hbase,...
    N_syn,HMap,TempSortByPeaks,NShuffles,H_all3D_unsorted,W_all3D_unsorted,NTime,NTime_method)

IHbaseSorted = SortByPeaks(Hbase,N_syn,NTime,NTime_method);
HbaseS = Hbase(:,IHbaseSorted);
HbaseS = HbaseS(:,HMap);%Change order to be consisted with our previous paper

NTimeCond = size(Hbase,1);

% sort temp patterns
% H_all3D_unsorted = reshape(permute(H_all,[4,3,2,1]),NTime,N_syn,[]); %make rectangular matrix with all patterns
% W_all3D_unsorted = reshape(permute(W_all,[4,3,2,1]),N_syn,NCondMuscles,[]); %make rectangular matrix with all patterns

H_all3D_s = H_all3D_unsorted; % sorted by temp base patterns
W_all3D_s = W_all3D_unsorted;
if TempSortByPeaks
    for ir = 1: NShuffles
        H  = squeeze(H_all3D_unsorted(1:NTimeCond,:,ir));%100x5
        IHSorted = SortByPeaks(H,N_syn);
        H_all3D_s(:,:,ir)= H_all3D_unsorted(:,IHSorted,ir);
        W_all3D_s(:,:,ir)= W_all3D_unsorted(IHSorted,:,ir);
    end
else
    MI = zeros( NShuffles,N_syn);
    MV = zeros( NShuffles,N_syn);
    idx1 = [];
    idx2 =[];
    for ir = 1: NShuffles
        H  = squeeze(H_all3D_unsorted(1:NTimeCond,:,ir));%100x5
        R2ir =corr(HbaseS,H).^2;
        [Max_Corr,Ind_Corr] = GetMax(R2ir);
        MI(ir,:) = Ind_Corr;
        MV(ir,:)= Max_Corr;
        H_all3D_s(:,:,ir)= H_all3D_unsorted(:,MI(ir,:),ir);
        W_all3D_s(:,:,ir)= W_all3D_unsorted(MI(ir,:),:,ir);

    end
    disp(strcat('MaxR2(mean/SD) for temp sort'))
    disp(strcat(num2str(mean(MV))))
    disp(strcat(num2str(std(MV))))
    disp(' ')
end

