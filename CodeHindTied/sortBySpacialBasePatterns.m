function       [WbaseS,H_all3D_ss,W_all3D_ss] = sortBySpacialBasePatterns(Wbase,N_syn,...
    NShuffles,H_all3D_unsorted,W_all3D_unsorted,HMap)
WbaseS = Wbase;% todo: sort spacial if needed
%
H_all3D_ss = H_all3D_unsorted; % sorted by spatial base patterns
W_all3D_ss = W_all3D_unsorted;


MI = zeros( NShuffles,N_syn);
MV = zeros( NShuffles,N_syn);
idx1 = []; idx2 =[];
for ir = 1: NShuffles % all sets and repetitions together
    W  = squeeze(W_all3D_unsorted(:,:,ir))';%16x5 -- one
    R2ir =corr(WbaseS,W).^2;
    [Max_Corr,Ind_Corr] = GetMax(R2ir);
    MI(ir,:) = Ind_Corr;
    MV(ir,:)= Max_Corr;
    H_all3D_ss(:,:,ir)= H_all3D_unsorted(:,MI(ir,:),ir);
    W_all3D_ss(:,:,ir)= W_all3D_unsorted(MI(ir,:),:,ir);

end
disp(strcat('MaxR2(mean/SD) for spatial sort '))
disp(strcat(num2str(mean(MV))))
disp(strcat(num2str(std(MV))))
disp(' ')

% find mean H  to sort by peaks of H:: for spacial sorting
Hmean = mean(H_all3D_ss,3);
[~, imaxHmean] = max(Hmean); %first max first after sorting
[~,IHmeanSorted] = sort(imaxHmean,'ascend');

% first temp syn has max at end
if(imaxHmean(IHmeanSorted(1))>10&&imaxHmean(IHmeanSorted(N_syn))>90)
    % if Hbase(1,IHbaseSorted(N_syn)')>0.5
    IHmeanSorted = circshift(IHmeanSorted,1);
end
HmeanS = Hmean(:,IHmeanSorted);

H_all3D_ss = H_all3D_ss(:,IHmeanSorted,:);
W_all3D_ss = W_all3D_ss(IHmeanSorted,:,:);

HmeanS = HmeanS(:,HMap);%Change order to be consisted with our previous paper
H_all3D_ss = H_all3D_ss(:,HMap,:);
W_all3D_ss = W_all3D_ss(HMap,:,:);
