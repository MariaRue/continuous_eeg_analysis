function dataAverage = select_single_subject_channel_data(channels,data)

    cfg = [];
    cfg.channel = channels;
    cfg.avgoverchan = 'yes';
    cfg.nanmean = 'yes';
    dataAverage = ft_selectdata(cfg,  data);


end