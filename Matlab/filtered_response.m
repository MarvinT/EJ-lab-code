function fun = filtered_response(spk_times, tau)
% spk_times are a [nx1] array

fun = @(t) sum(exp((-(repmat(spk_times, size(t)) - repmat(t, size(spk_times))) .^ 2) ./ (2 .* tau .^2)), 1);
%fun = @(t) arrayfun(fun2, t);
%fun = @(t) response_function(t, spk_times, tau);
end

function rsp = response_function(t, spk_times, tau)
rsp = sum(exp((-(repmat(spk_times, size(t)) - repmat(t, size(spk_times))) .^ 2) ./ (2 .* tau .^2)), 1);
end