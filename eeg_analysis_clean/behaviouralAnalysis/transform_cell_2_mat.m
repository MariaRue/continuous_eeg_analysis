function [cellNew] = transform_cell_2_mat(cell)
cellNew = [];
for cellCount = 1:length(cell)
    
    cellNew = [cellNew;cell{cellCount}];
    
end 

end 
