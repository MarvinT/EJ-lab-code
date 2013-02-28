function sig_str = pop_motion_signal(velocity, spikes, indices1, indices2, x_pos, trigger, stop, tau)
%                   pop_motion_signal(v(i), datarun{2}.spikes,
%                   cell_indices1, cell_indices2, cell_x_pos, tr(k), stop,
%                   run_opt.tau)

pairs = zeros(2, length(indices2) * (length(indices2) - 1) / 2);
counter = 1;
for i = 2:length(indices2)
    for j = 1:i-1
        pairs(:,counter) = [i; j];
        counter = counter + 1;
    end
end

% trial start @ t=0
spks = cellfun(@(s) s - trigger, spikes(indices2), 'UniformOutput', false);
% only consider spikes that occured in the trial
spks = cellfun(@(s) s(s >= 0 & s <= stop), spks, 'UniformOutput', false);
% circularly shift spikes by dt
spks_shifted_L = cellfun(@(s, x) s - x / velocity, spks, num2cell(x_pos(indices1)), 'UniformOutput', false);
spks_shifted_L = cellfun(@(s) mod(s, stop), spks_shifted_L, 'UniformOutput', false);
spks_shifted_R = cellfun(@(s, x) s + x / velocity, spks, num2cell(x_pos(indices1)), 'UniformOutput', false);
spks_shifted_R = cellfun(@(s) mod(s, stop), spks_shifted_R, 'UniformOutput', false);
% replicate spikes before and after trial to minimize artifacts of spikes
% shifting circularly across the border
spks_shifted_L = cellfun(@(s) [s - stop; s; s + stop], spks_shifted_L, 'UniformOutput', false);
spks_shifted_R = cellfun(@(s) [s - stop; s; s + stop], spks_shifted_R, 'UniformOutput', false);
% filter responses
flt_rsp_shifted_L = cellfun(@(s) filtered_response(s, tau), spks_shifted_L, 'UniformOutput', false);
flt_rsp_shifted_R = cellfun(@(s) filtered_response(s, tau), spks_shifted_R, 'UniformOutput', false);
% and calculate the integral
str = cellfun(@(frL1, frR1, frL2, frR2) ... 
    integral(@(t) frR1(t) .* frR2(t) - frL1(t) .* frL2(t), 0, stop, 'AbsTol', 1e-3, 'RelTol', 1e-4), ...
    flt_rsp_shifted_L(pairs(1,:)), flt_rsp_shifted_R(pairs(1,:)), flt_rsp_shifted_L(pairs(2,:)), flt_rsp_shifted_R(pairs(2,:)));
% and sum them
sig_str = sum(str);

end