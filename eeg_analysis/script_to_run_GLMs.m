%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; 
design_matrix_type = 'absolute_coherences';
source_density = 1;
step_3_continuous_GLM(design_matrix_type,source_density)

clear all; 
design_matrix_type = 'jumps_and_jump_PEs';
source_density = 1;
step_3_continuous_GLM(design_matrix_type,source_density)

clear all; 
design_matrix_type = 'jumps_plus_absolute';
source_density = 1;
step_3_continuous_GLM(design_matrix_type,source_density)

clear all; 
design_matrix_type = 'split_PE';
source_density = 1;
step_3_continuous_GLM(design_matrix_type,source_density)



