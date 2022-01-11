function create_topo_plot(data,startTime, endTime, zlim, electrodesForPermTest,colour)
    
        cfg = [];
        cfg.highlight = 'on';
        cfg.highlightchannel = electrodesForPermTest;
        cfg.highlightsymbol = '^';
        cfg.xlim = [startTime endTime];  % time window for which we create a topo plot
        cfg.zlim = zlim; %colorbar scale
        cfg.layout = 'easycapM1.mat';
        cfg.colormap = flip(colour);
   
        ft_topoplotER(cfg,data); colorbar
     
end 