function [Wbase,Wbase_set] = findBasicSpacialPatterns(W_all,NCondMuscles,N_syn,NHbase_repeat)
tic
%  finds basic spacial patterns by using nnmf on already found patterns
W_all_unsorted = reshape(permute(W_all,[3,4,2,1]),NCondMuscles,[]); %make rectangular matrix with all patterns
Wbase_set = zeros(NCondMuscles,N_syn,NHbase_repeat);
tic
for ir = 1:NHbase_repeat
    [Wbase,~] = nnmf(W_all_unsorted,N_syn,'replicates',5); %Finds basic patterns
    Wbase_set(:,:,ir)= Wbase./max(Wbase);
end
disp(strcat('spatial basic pattern found '))
toc
        % find if spatial basic patterns split in several sets
        Wbase1 = squeeze(Wbase_set(:,:,1));
        MI = zeros(NHbase_repeat,N_syn);
        MV= MI;
        idx1 =[]; 
        idx2 = [];
        Wbase_setSorted = Wbase_set;
        for ir = 1:NHbase_repeat
            [MV(ir,:), MI(ir,:)] = max(corr(Wbase1,squeeze(Wbase_set(:,:,ir))).^2);
            [~,idxu,idxMI] = unique(MI(ir,:));% Unique values
            if numel(idxu)<N_syn
                idx2 = [idx2,ir];
            else
                idx1 = [idx1,ir];
                Wbase_setSorted(:,MI(ir,:),ir)= Wbase_set(:,:,ir);
            end
        end
        if isempty(idx2) % set didn't split
            Wbase = mean(Wbase_setSorted,3);
            disp('Basic spatial set did"t splits');
        else
            Wbase = mean(Wbase_setSorted(:,:,ismember(1:NHbase_repeat,idx2)),3);
            disp('Basic spatial set splits');
        end
        disp(strcat('MaxR2(mean/SD) for finding spatial base '))
        disp(strcat(num2str(mean(MV(2:end,:)))))
        disp(strcat(num2str(std(MV(2:end,:)))))

        disp(' ')