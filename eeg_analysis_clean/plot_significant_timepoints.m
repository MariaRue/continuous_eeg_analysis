function plot_significant_timepoints(time,SignificantTimePoints,scale)

    plot(time(SignificantTimePoints), ones(length(SignificantTimePoints),1).* scale, '*k')

end 