function data = select_data_for_reliability( allBetas, selSensors, conIdx, chanlabels )

sensIdx = match_strings(chanlabels, selSensors);
data = allBetas(:, sensIdx, :, :, conIdx);

% average across blocks of the condition (if needed), then across sensors
if numel(conIdx) > 1
    data = mean(data, 5);
end
data = squeeze(mean(data, 2));
data = permute(data, [1 3 2]);

end