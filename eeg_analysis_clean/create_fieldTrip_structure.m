function data_ft = create_fieldTrip_structure(timeBins, new_labels, data, elecStructure)


data_ft.time = timeBins; 
data_ft.label = new_labels; 
data_ft.avg = squeeze(data); 
data_ft.dimord = 'chan_time'; 
data_ft.elec = elecStructure; 





end 