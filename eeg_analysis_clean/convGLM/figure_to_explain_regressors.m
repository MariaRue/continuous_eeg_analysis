cl = cbrewer('seq','Greens',36);
cl = cl([6:6:36],:); 
colormap(cl)
x = linspace(1,10,100);
y = [ones(10,1)*-1;ones(10,1);ones(10,1)*3;ones(10,1)*-2;ones(10,1)*1;ones(10,1)*2;ones(10,1)*3;ones(10,1)*1;ones(10,1)*-3;ones(10,1)*-2];
 y2=y(1:2:end);
 x2 = x(1:2:end);
 
jumpIdx = y(2:end) ~= y(1:end-1);
jump = zeros(length(x),1); 
jump(find(jumpIdx)+1) =1; 
x3 = x(find(jump)); 
y3 = y(find(jump));

figure 
subplot(5,1,1)
hold on 
plot(x,y,'k','LineWidth',3) 
plot(x2,y2,'gx','MarkerSize',8)
plot(x3,y3,'rx','MarkerSize',10)
ylim([-3.5 3.5])
xticklabels(0:9)
tidyfig; 
% jump image sc 
subplot(5,1,2)
jumpImage(1,:) = jump; 
jumpImage(2,:) = [0,jump(1:end-1)'];
jumpImage(3,:) = [0,0,jump(1:end-2)'];

imagesc(jumpImage)
%colorbar; colormap(cl); caxis([0 5]); 
xticklabels(1:9)
caxis([0 5])
tidyfig; 

subplot(5,1,3)
level = abs(y(logical(jump))); 
cohlevel = zeros(length(x),1); 
cohlevel(find(jump)) = level; 

levelImage(1,:) = cohlevel'; 
levelImage(2,:) = [0,cohlevel(1:end-1)']; 
levelImage(3,:) = [0,0,cohlevel(1:end-2)']; 
imagesc(levelImage); %colorbar; colormap(cl); caxis([0 5]); 
xticklabels(1:9)
caxis([0 5])
tidyfig; 



subplot(5,1,4)
pes_diff = abs(diff([0,y(logical(jump))']));
pes = zeros(length(x),1);
pes(find(jump)) = pes_diff; 

PEImage(1,:) = pes'; 
PEImage(2,:) = [0,pes(1:end-1)'];
PEImage(3,:) = [0,0,pes(1:end-2)'];

imagesc(PEImage); %colorbar; colormap(cl); caxis([0 5]); 
xticklabels(1:9)
caxis([0 5])
tidyfig; 


subplot(5,1,5)
absY = abs(y); 

absImage(1,:) = absY'; 
absImage(2,:) = [0,absY(1:end-1)']; 
absImage(3,:) = [0,0,absY(1:end-2)']; 

imagesc(absImage); colorbar; colormap(cl); caxis([0 5]); 
xticklabels(1:9)
tidyfig; 