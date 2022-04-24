function [Model] = calculate_integration_kernel_based_on_model(params, nS)

for subject = 1:nS
    
    for condition = 1:4
        
        Model{subject,condition} = calculate_exponential_model_for_estimated_parameters(params.parameters(condition,:,subject),params.time{condition,subject});
    
    end 
    
end 


end 