data = squeeze(betas{3}(40,:,:,1)); 
mean_data = mean(data,2);


smoothed_toplot = conv2(ones(10,1)/10,1,mean_data,'same');

subplot(2,3,6)
plot(time_idx(3).timebins,smoothed_toplot)
title('subject 11')


%% plot data for 10min runs 

data = squeeze(betas{3}(40,:,:,1)); 
  set1 = 1;
    set2 = 2;
for t = 1:3

    mean_data = mean(data(:,set1:set2),2);
    
    
subplot(3,1,t)
smoothed_toplot = conv2(ones(10,1)/10,1,mean_data,'same');


plot(time_idx(3).timebins,smoothed_toplot)
title('subject 12')
      set1 = set1+2;
    set2 = set2+2;
end 