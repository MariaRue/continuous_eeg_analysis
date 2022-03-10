dataFreqShort = corrResponseSelectedData.subjectLevel.FreqShort; 
dataFreqLong = corrResponseSelectedData.subjectLevel.FreqLong; 
dataRareShort = corrResponseSelectedData.subjectLevel.RareShort; 
dataRareLong = corrResponseSelectedData.subjectLevel.RareLong; 


%% 

figure
for sj = 1:24 
    hold on 
    
    plot(dataRareLong{sj}.time, dataRareLong{sj}.avg); 
    
    hold off 
end 
title('Rare Long')
xlabel('time') 
ylabel('effect size')

figure
for sj = 1:24 
    hold on 
    
    plot(dataRareShort{sj}.time, dataRareShort{sj}.avg); 
    
    hold off 
end
title('Rare Short')
xlabel('time') 
ylabel('effect size')

figure
for sj = 1:24 
    hold on 
    
    plot(dataFreqShort{sj}.time, dataFreqShort{sj}.avg); 
    
    hold off 
end 
title('Freq Short')
xlabel('time') 
ylabel('effect size')

figure
for sj = 1:24 
    hold on 
    
    plot(dataFreqLong{sj}.time, dataFreqLong{sj}.avg); 
    
    hold off 
end 
title('Freq Long')
xlabel('time') 
ylabel('effect size')

%% single subject 

figure

for sj = 1:24 
   subplot(4,6,sj)
   
 

    plot(dataFreqShort{sj}.time, dataFreqShort{sj}.avg); 
    
     if sj == 1 
       title('freq short')
       
   else 
      
        title(num2str(sj))
    
       
   end 
    xlabel('time') 
ylabel('effect size')

    
end 

figure

for sj = 1:24 
   subplot(4,6,sj)
   

   
    plot(dataFreqLong{sj}.time, dataFreqLong{sj}.avg); 
    
  
    
    xlabel('time') 
ylabel('effect size')
   if sj == 1
 title('freq long')

   else 
      
        title(num2str(sj)) 
   end 
    
end 

figure

for sj = 1:24 
   subplot(4,6,sj)
   


    plot(dataRareShort{sj}.time, dataRareShort{sj}.avg); 
    
   if sj == 1
      title('rare short') 
   else 
    title(num2str(sj))   
   end
    
    xlabel('time') 
ylabel('effect size')

    
end 

figure

for sj = 1:24 
   subplot(4,6,sj)
   


    plot(dataRareLong{sj}.time, dataRareLong{sj}.avg); 
    

       if sj == 1
     title('rare long')  
   else 
       
      title(num2str(sj))     
   end 
    xlabel('time') 
ylabel('effect size')

    
end 



