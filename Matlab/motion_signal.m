function [str, flt_rsp1, flt_rsp2, flt_rsp1_shifted, flt_rsp2_shifted, spks_1_shifted, spks_2_shifted] = motion_signal(velocity, spks_1, spks_2, dx, trigger, trial_length, tau)
% trial start @ t=0
spks_1 = spks_1 - trigger;
spks_2 = spks_2 - trigger;
% only consider spikes that occured in the trial
spks_1 = spks_1(spks_1 >= 0 & spks_1 <= trial_length);
spks_2 = spks_2(spks_2 >= 0 & spks_2 <= trial_length);

if abs(dx / velocity) > trial_length / 2 ||isempty(spks_1) || isempty(spks_2)
    str = 0;
    flt_rsp1 = @(t) 0;
    flt_rsp2 = @(t) 0;
    flt_rsp1_shifted = @(t) 0;
    flt_rsp2_shifted = @(t) 0;
    spks_1_shifted = [];
    spks_2_shifted = [];
    return
end

% circularly shift spikes by dt
spks_1_shifted = spks_1 - dx / velocity;
spks_1_shifted = sort(mod(spks_1_shifted, trial_length));
spks_2_shifted = spks_2 - dx / velocity;
spks_2_shifted = sort(mod(spks_2_shifted, trial_length));
% replicate spikes before and after trial to minimize artifacts of spikes
% shifting circularly across the border
spks_1_shifted = [spks_1_shifted(ceil(end/2):end) - trial_length; spks_1_shifted; spks_1_shifted(1:floor(end/2)) + trial_length];
spks_2_shifted = [spks_2_shifted(ceil(end/2):end) - trial_length; spks_2_shifted; spks_2_shifted(1:floor(end/2)) + trial_length];
% filter responses
flt_rsp1 = filtered_response(spks_1, tau);
flt_rsp2 = filtered_response(spks_2, tau);
flt_rsp1_shifted = filtered_response(spks_1_shifted, tau);
flt_rsp2_shifted = filtered_response(spks_2_shifted, tau);
% and return the integral
str = integral(@(t) flt_rsp1(t) .* flt_rsp2_shifted(t) - flt_rsp2(t) .* flt_rsp1_shifted(t), 0, trial_length, ...
    'AbsTol', 1e-1, 'RelTol', 1e-2);
end