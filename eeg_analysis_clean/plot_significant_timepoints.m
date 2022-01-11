function plot_significant_timepoints(time,SignificantTimePoints)

    plot(time(SignificantTimePoints), ones(length(SignificantTimePoints),1), '*k')

end 