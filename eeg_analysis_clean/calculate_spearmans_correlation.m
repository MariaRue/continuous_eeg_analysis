function [corr_across_subjects] = calculate_spearmans_correlation(data,allTau_reduced_dm,nS)

%first, rank data across subjects
[~,ai] = sort(allTau_reduced_dm);sl_sorted = [1:nS]'; sl_sorted(ai) = sl_sorted; %to use ranks
sl_sorted_dm = sl_sorted - mean(sl_sorted);

    data_ranked = data(:,:);
    [~,ai] = sort(data_ranked,'ascend');
    tmp = 1:nS;
    for i = 1:size(data_ranked,2); 
        data_ranked(ai(:,i),i) = tmp; 
    end
    rrr = corr(data_ranked,sl_sorted_dm);
    
    corr_across_subjects = reshape(rrr,[size(data,2),size(data,3), size(data,4)]);

end 