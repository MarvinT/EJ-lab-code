function sig_str = pop_motion_signal(velocity, spikes, indices1, indices2, x_pos, trigger, stop, tau, tol)
%                   pop_motion_signal(v(i), datarun{2}.spikes,
%                   cell_indices1, cell_indices2, cell_x_pos, tr(k), stop,
%                   run_opt.tau)
if nargin < 9
    tol = 1e-3;
end

sig_str = 0;

pairs = zeros(2, length(indices2) * (length(indices2) - 1) / 2);
counter = 1;
for i = 2:length(indices2)
    for j = 1:i-1
        pairs(:,counter) = [i; j];
        counter = counter + 1;
    end
end

for j = 1:length(pairs)
    spks_1 = spikes{indices2(pairs(1,j))};
    spks_2 = spikes{indices2(pairs(2,j))};
    dx = x_pos(indices1(pairs(2,j))) - x_pos(indices1(pairs(1,j)));
    sig_str = sig_str + motion_signal(velocity, spks_1, spks_2, dx, trigger, stop, tau, tol);
end

end