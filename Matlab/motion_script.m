% This script calculates the motion signals for a data run displays several
% nice guis for visualization.
% Marvin Thielk 2013
% mthielk@salk.edu


% vary the options set in run_opt to turn different behavior on and off
run_opt.load = true; % T/F
run_opt.remote = true; % T/F
run_opt.data_run = 19; % 12-19
run_opt.cell_type = 'Off parasol'; % on/off parasol, on/off midget
run_opt.config_num = 1; % 1-4
run_opt.raster = false; % T/F
run_opt.trial_raster = false; % T/F
run_opt.trial_raster_shift = false; % T/F
run_opt.manual_speed_tuning = false; % T/F
run_opt.velocity_lim = 150; % >0
run_opt.auto_speed_tuning = false; % T/F
run_opt.tau = .01; % tuning parameter
run_opt.pop_speed_tuning = false; % T/F
run_opt.tol = 1e-3;
run_opt.savefig = true; % T/F
run_opt.trial_num = 1; % > 0
run_opt.trial_estimate = false; % T/F
run_opt.auto_set = true; % T/F -- note: modifies run_opt
run_opt.trial_estimate_start = 120;
run_opt.data_run_plots = true; % T/F

if run_opt.auto_set
    if run_opt.data_run == 17
        run_opt.velocity_lim = 50;
        run_opt.config_num = 2;
        run_opt.trial_estimate_start = 14.6;
        run_opt.tol = 1e-2;
    elseif run_opt.data_run == 18
        run_opt.velocity_lim = 150;
        run_opt.config_num = 1;
        run_opt.trial_estimate_start = 110;
        run_opt.tol = 1e-3;
    elseif run_opt.data_run == 19
        run_opt.velocity_lim = 300;
        run_opt.config_num = 1;
        run_opt.trial_estimate_start = 203;
        run_opt.tol = 1e-4;
    end
end

if exist('export_fig', 'file') == 7
    addpath export_fig
end

if run_opt.load %load data

    clear datarun tr

    if run_opt.remote 
        datarun{1}.names.rrs_params_path='/snle/analysis/2007-03-27-1/data011-nwpca/data011-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/snle/analysis/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/braid/snle/analysis-archive/Experiments/Array/Analysis/2007-03-27-1/stimuli/s%d', run_opt.data_run);
    else
        datarun{1}.names.rrs_params_path='/Data/2007-03-27-1/data011-nwpca/data011-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Data/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/braid/snle/analysis-archive/Experiments/Array/Analysis/2007-03-27-1/stimuli/s%d', run_opt.data_run);
    end
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);

    datarun=load_data(datarun,opt);

    datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true)); 

    datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0); 
end

tic;

if run_opt.raster || run_opt.trial_raster || ...
        run_opt.trial_raster_shift || run_opt.manual_speed_tuning || ...
        run_opt.auto_speed_tuning || run_opt.pop_speed_tuning || ...
        run_opt.trial_estimate || run_opt.data_run_plots
    clf; set(gcf, 'color', 'white');
    
    cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
    cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
    
    cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits);
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1));
    
    cell_indices1 = cell_indices1(cell_sort_idx);
    cell_indices2 = cell_indices2(cell_sort_idx);
    
    start = 0;
    stop = mean(datarun{2}.triggers(2:2:end) - datarun{2}.triggers(1:2:end));

    tr=datarun{2}.triggers(1:2:end); % triggers mark the beginning and end
    t=find(datarun{2}.stimulus.trial_list==run_opt.config_num);
    tr=tr(t);
end

if run_opt.raster %raster

    k=1; kmin=1; kmax=length(cell_indices2); hk=loop_slider(k,kmin,kmax);

    while k
        if ~ishandle(hk)
            break
        end
        k=round(get(hk,'Value')); 
 
        psth_raster(start,stop,datarun{2}.spikes{cell_indices2(k)}',tr);
        title(sprintf('%d %.2f', datarun{2}.cell_ids(cell_indices2(k)), datarun{1}.vision.sta_fits{cell_indices1(k)}.mean(1) ))

        uiwait;
    end
end

if run_opt.trial_raster
    k=1; kmin=1; kmax=length(tr); hk=loop_slider(k,kmin,kmax);
    while k
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        
        shifted_trial_raster(start, stop, datarun{2}.spikes, cell_indices2, tr(k))
        title(sprintf('trial number %d', k))
        
        uiwait;
    end
end

if run_opt.trial_raster_shift
    k=1; kmin=1; kmax=length(tr); hk=loop_slider_n(k,kmin,kmax,1);
    j=0; jmin=-200; jmax=-jmin; hj=loop_slider_n(j,jmin,jmax,2);
    while true
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        j = round(get(hj, 'Value'));
        
        trigger = tr(k) + (cell_x_pos(cell_indices1) - cell_x_pos(cell_indices1(1))) ./ j;
        shifted_trial_raster(start, stop, datarun{2}.spikes, cell_indices2, trigger)
        title(sprintf('trial number %d with putative speed %d stixels/sec', k, j))
        
        uiwait;
    end
end

if run_opt.manual_speed_tuning
    k=1; kmin=1; kmax=length(tr); hk=loop_slider_n(k,kmin,kmax,1);
    jmin=1; jmax=length(cell_indices2); j=round(.15 * jmax); hj=loop_slider_n(j,jmin,jmax,2);
    hmin=1; hmax=length(cell_indices2); h=round(.9 * hmax); hh=loop_slider_n(h,hmin,hmax,3);
    v=15; vmin=-run_opt.velocity_lim; vmax=run_opt.velocity_lim; hv=loop_slider_n(v,vmin,vmax,4);
    while true
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        j = round(get(hj, 'Value'));
        h = round(get(hh, 'Value'));
        v = get(hv, 'Value');
        
        spks_1 = datarun{2}.spikes{cell_indices2(j)};
        spks_2 = datarun{2}.spikes{cell_indices2(h)};
        dx = cell_x_pos(cell_indices1(h)) - cell_x_pos(cell_indices1(j));
        [sig_str, flt_rsp1, flt_rsp2, flt_rsp1_shifted, flt_rsp2_shifted, spks_1_shifted, spks_2_shifted] = motion_signal(v, spks_1, spks_2, dx, tr(k), stop, run_opt.tau);
        t = linspace(0, stop, 200);
        
        lims = [0 stop 0 4];
        
        subplot(3, 1, 1)
        plot(t, flt_rsp1(t), 'b-', t, flt_rsp2_shifted(t), 'g--', spks_2_shifted, ones(size(spks_2_shifted)), 'g.')
        title(sprintf('motion signal strength between neurons %d and %d in trial number %d', cell_indices2(j), cell_indices2(h), k))
        legend('filtered response of 1st cell', 'filtered response of 2nd cell shifted')
        axis(lims)
        
        subplot(3, 1, 2)
        plot(t, flt_rsp2(t), 'g-', t, flt_rsp1_shifted(t), 'b--', spks_1_shifted, ones(size(spks_1_shifted)), 'b.')
        title(sprintf('velocity = %d   and dx = %d', v, dx))
        legend('filtered response of 2nd cell', 'filtered response of 1st cell shifted')
        axis(lims)
        
        subplot(3, 1, 3)
        plot(t, flt_rsp1(t) .* flt_rsp2_shifted(t) - flt_rsp2(t) .* flt_rsp1_shifted(t))
        title(sprintf('signal strength: %d',sig_str))
        
        uiwait;
    end
end

if run_opt.auto_speed_tuning
    k=1; kmin=1; kmax=length(tr); hk=loop_slider_n(k,kmin,kmax,1);
    jmin=1; jmax=length(cell_indices2); j=round(.15 * jmax); hj=loop_slider_n(j,jmin,jmax,2);
    hmin=1; hmax=length(cell_indices2); h=round(.9 * hmax); hh=loop_slider_n(h,hmin,hmax,3);
    while true
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        j = round(get(hj, 'Value'));
        h = round(get(hh, 'Value'));
        
        v = linspace(1, run_opt.velocity_lim);
        
        sig_str = zeros(size(v));
        for i = 1:length(v)
            spks_1 = datarun{2}.spikes{cell_indices2(j)};
            spks_2 = datarun{2}.spikes{cell_indices2(h)};
            dx = cell_x_pos(cell_indices1(h)) - cell_x_pos(cell_indices1(j));
            sig_str(i) = motion_signal(v(i), spks_1, spks_2, dx, tr(k), stop, run_opt.tau);
        end
        
        subplot(2,1,1)
        plot(v, sig_str)
        title(sprintf('motion signal strength between neurons %d and %d in trial number %d', cell_indices2(j), cell_indices2(h), k))
        xlabel('velocity')
        ylabel('net rightward motion signal')
        
        t = linspace(tr(k), tr(k)+stop, 500);
        trial_spks1 = spks_1(spks_1 >= tr(k) & spks_1 <= (tr(k) + stop));
        trial_spks2 = spks_2(spks_2 >= tr(k) & spks_2 <= (tr(k) + stop));
        flt_rsp1 = filtered_response(trial_spks1, run_opt.tau);
        flt_rsp2 = filtered_response(trial_spks2, run_opt.tau);
        subplot(4,1,3)
        plot(t, flt_rsp1(t), 'b', trial_spks1, ones(size(trial_spks1)), 'k.')
        lims = [tr(k) (tr(k) + stop) 0 max(max(flt_rsp1(t)), max(flt_rsp2(t))) * 1.1];
        axis(lims);
        title(sprintf('Average ISI: %d', (trial_spks1(end) - trial_spks1(1)) / (length(trial_spks1) - 1)));
        xlabel('time')
        
        subplot(4,1,4)
        plot(t, flt_rsp2(t), 'b', trial_spks2, ones(size(trial_spks2)), 'r.')
        axis(lims);
        title(sprintf('Average ISI: %d', (trial_spks2(end) - trial_spks2(1)) / (length(trial_spks2) - 1)));
        xlabel('time')
        
        uiwait;
    end
end

if run_opt.pop_speed_tuning
    if matlabpool('size') <= 0
        matlabpool
    end
    
    v = linspace(1, run_opt.velocity_lim, 50);
    
    sig_str = zeros(size(v));
    parfor i = 1:length(v)
        sig_str(i) = pop_motion_signal(v(i), datarun{2}.spikes, cell_indices1, cell_indices2, cell_x_pos, tr(run_opt.trial_num), stop, run_opt.tau, run_opt.tol*.1);
        
        fprintf('*')
    end
    fprintf('\n')
    plot(v, sig_str)
    xlabel('velocity')
    ylabel('net rightward motion signal')
    title(sprintf('motion signal strength  in trial number %d', run_opt.trial_num))
    if run_opt.savefig
        export_fig(sprintf('figs/%s_data_run_%d_config_%d_trial_%d', run_opt.cell_type, run_opt.data_run, run_opt.config_num, run_opt.trial_num), '-png', '-r300', '-painters')
    end
end

if run_opt.trial_estimate
    if matlabpool('size') <= 0
        matlabpool
    end
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 30, 'LargeScale', 'off');
    estimates = zeros(size(tr));
    parfor i = 1:length(tr)
        estimates(i) = fminunc(@(v) -pop_motion_signal(v, datarun{2}.spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
        fprintf('for trial %d, the estimated speed was %d', i, estimates(i))
    end
    figure()
    hist(estimates)
    xlabel('speed estimate (stixels/sec)')
    ylabel('trials')
    title(sprintf('%s data run %d config %d', run_opt.cell_type, run_opt.data_run, run_opt.config_num))
    if run_opt.savefig
        export_fig(sprintf('figs/%s_data_run_%d_config_%d', run_opt.cell_type, run_opt.data_run, run_opt.config_num), '-png', '-r300', '-painters')
    end
end

if run_opt.data_run_plots
    if matlabpool('size') <= 0
        matlabpool
    end
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol, 'MaxFunEvals', 30, 'LargeScale', 'off');
    cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
    for j=1:length(cell_types)
        cell_type = cell_types{j};
        cell_indices1=get_cell_indices(datarun{1},{cell_type});
        cell_indices2=get_cell_indices(datarun{2},{cell_type});
        
        cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits);
        [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1));
        
        cell_indices1 = cell_indices1(cell_sort_idx);
        cell_indices2 = cell_indices2(cell_sort_idx);
        
        estimates = zeros(size(tr));
        parfor i = 1:length(tr)
            estimates(i) = fminunc(@(v) -pop_motion_signal(v, datarun{2}.spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
            fprintf('for %s cells, in trial %d, the estimated speed was %d', cell_type, i, estimates(i))
        end
        figure()
        hist(estimates)
        xlabel('speed estimate (stixels/sec)')
        ylabel('trials')
        title(sprintf('%s data run %d config %d', cell_type, run_opt.data_run, run_opt.config_num))
        if run_opt.savefig
            export_fig(sprintf('figs/%s_data_run_%d_config_%d', cell_type, run_opt.data_run, run_opt.config_num), '-png', '-r300', '-painters')
        end
    end
end

toc

matlabpool close