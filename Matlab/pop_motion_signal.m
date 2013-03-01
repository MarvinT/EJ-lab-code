function sig_str = pop_motion_signal(velocity, spikes, indices1, indices2, x_pos, trigger, stop, tau)
%                   pop_motion_signal(v(i), datarun{2}.spikes,
%                   cell_indices1, cell_indices2, cell_x_pos, tr(k), stop,
%                   run_opt.tau)

flt_rsp_shifted_L = cell(size(indices2));
flt_rsp_shifted_R = cell(size(indices2));

for i=1:length(indices2)
    s = spikes{indices2(i)} - trigger;% trial start @ t=0
    s = s(s >= 0 & s <= stop);% only consider spikes that occured in the trial
    % circularly shift spikes by dt
    sL = mod(s - x_pos(indices1(i)) / velocity, stop);
    sR = mod(s + x_pos(indices1(i)) / velocity, stop);
    % replicate spikes before and after trial to minimize artifacts of spikes
    % shifting circularly across the border
    sL = [sL - stop; sL; sL + stop];
    sR = [sR - stop; sR; sR + stop];
    % filter responses
    flt_rsp_shifted_L{i} = filtered_response(sL, tau);
    flt_rsp_shifted_R{i} = filtered_response(sR, tau);
end

sig_str = 0;
for i = 2:length(indices2)
    for j = 1:i-1
        sig_str = sig_str + integral(@(t) flt_rsp_shifted_R{i}(t) .* flt_rsp_shifted_R{j}(t) - flt_rsp_shifted_L{i}(t) .* flt_rsp_shifted_L{j}(t), 0, stop, 'AbsTol', 1e-3, 'RelTol', 1e-4);
    end
end
end