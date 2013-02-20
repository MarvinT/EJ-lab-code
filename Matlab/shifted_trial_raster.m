function shifted_trial_raster(start, stop, spikes, cell_idxs, trigger)

if length(trigger) == 1
    trigger = repmat(trigger, size(cell_idxs));
end

psth_r = [];

for i = 1:length(cell_idxs)
    spk_times = spikes{cell_idxs(i)} - trigger(i);
    spk_idxs = find(spk_times >= start & spk_times <= stop);
    psth_r = [psth_r; 1000*spk_times(spk_idxs), repmat(length(cell_idxs)-i,[length(spk_idxs),1]);];
end

if ~isempty(psth_r)
    plot(psth_r(:,1),psth_r(:,2),'k.');%, 'MarkerSize',10
    axis([start*1000 stop*1000 0 length(cell_idxs)]);
end

end