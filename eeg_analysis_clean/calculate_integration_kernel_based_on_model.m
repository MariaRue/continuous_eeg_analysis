function [Model] = calculate_integration_kernel_based_on_model(params, selectedSubjects)

for subject = 1:length(selectedSubjects)
    
    for condition = 1:4
        
        Model(subject,condition,:) = calculate_exponential_model_for_estimated_parameters(params.parameters(condition,:,selectedSubjects(subject)),params.time);
    
    end 
    
end 


end 