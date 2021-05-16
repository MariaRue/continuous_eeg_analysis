function [se] = calculate_se_convGLM(data)

for t = 1:length(data{1}.time)
    for chan = 1:length(data{1}.label)
        
        for sj = 1:length(data)
            
            points_for_se_calc_at_t(sj) = data{sj}.avg(chan,t);
            
            
        end
        se(chan,t) = nanstd(points_for_se_calc_at_t)./sqrt(length(data));
        
    end
end
end
