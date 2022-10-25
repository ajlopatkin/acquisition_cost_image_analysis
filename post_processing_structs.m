% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post_processing_structs.m
%   AUTHOR: AJL
%   DATE: 2021_01_06
%   UPDATE: 2021_04_28
%   DESCRIPTION: Process timecourses of scanner acquisition cost data
%
%   Use this script to clean up time courses of colony growth for scanner
%   image data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all,
clear, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Define plasmids to analyze and saving directory
pth = pwd + "/pre_processed_data/";
saving_dir = pth + "filtered_structs/";

if ~isdir(saving_dir), mkdir(saving_dir), end
s = dir(pth);
s = s(contains({s.name}',".mat"));


%%%% Flags for analysis
plotting_flag = 1;
save_data_flag = 1;
separate_flag = 0;
overwrite_all_data = 0;
col_max = 40;
thresh = 0.8;
clr = ['bc';'rm'];
count = 1;

if overwrite_all_data
    delete(saving_dir + "*")
end

if plotting_flag
    figure, hold on
    num_row = 10;
    num_col = ceil(length(s)/num_row);
end

%%%% Loop through all plasmids
for p = 1:length(s)
    
    
    read_filename = s(p).folder + "/" + s(p).name;
    load(read_filename);
    write_filename = saving_dir + s(p).name;
    
    if ~overwrite_all_data & isfile(write_filename)
        disp("Skipping " + s(p).name)
        continue
    end
    
    avg_baseline = find_baseline(colony_struct, spot_metadata, col_max);
    
    disp("Currently analyzing file " + s(p).name)
    ind_remove = [];
    spot_remove = [];
    for q = 1:length(colony_struct)
        
        
        curr_spot_id = colony_struct(q).spot_id;
        colnys_in_spot = spot_metadata.real_colony_count(find(spot_metadata.spot_id == curr_spot_id));
        if colnys_in_spot >= col_max & ~contains(s(p).name,'NR1') & ~contains(s(p).name,'RP4')
            ind_remove(end+1) = q;
            spot_remove(end+1) = curr_spot_id;
            continue
        end
        
        curr_timecourse = colony_struct(q).timecourse;
        condition = colony_struct(q).condition;
        time=colony_struct.time;
        processed_timecourse = process_timecourse(curr_timecourse, avg_baseline);
        
        % log transform and smooth out initial points
        log_timecourse = log(processed_timecourse./min(processed_timecourse));
        log_timecourse(1:find((log_timecourse(1:end-20)) == 0, 1, 'last')) = 0;
        
        
        if std(log_timecourse) == 0 | any(isnan(log_timecourse)) | any(~isreal(log_timecourse))
            ind_remove(end+1) = q;
            continue
        end
        % fit curve using estimated values
        lag_guess = time(find(log_timecourse>0,1,'first')); if lag_guess == 1, lag_guess=10; end
        init_guess = [max(log_timecourse), .4, lag_guess];
        opts = optimset('Display','off');
        x = lsqcurvefit(@growth_logistic,init_guess,time',log_timecourse,[],[],opts);
        time2 = linspace(time(1),time(end),10000);
        fitted_timecourse = growth_logistic(x,time2);
        
        
        if contains(s(p).name,'R388') | contains(s(p).name,'pRK100') |...
                contains(s(p).name,'R1drd') | contains(s(p).name,'Rs-a') |...
                contains(s(p).name,'R57b')
            xlag = 48;
        else
            xlag = 24;
        end
        
        if contains(s(p).name,'R1')
            diff_thresh = 3000;
        else
            diff_thresh = 350;
        end
        
        
        if x(3) < 1 | x(3) > xlag | ...                             % reasonable lag time
                all(fitted_timecourse < thresh) | ...               % doesn't reach threshold
                max(diff(processed_timecourse)) > diff_thresh | ... % poor quality curve
                max(abs(diff(log_timecourse))) > 1.5 | ...          % poor quality curve
                sum(diff(processed_timecourse(1:find(time>lag_guess,1,'first')))<0) > 20 | ...  % unreliable growth start
                mean(fitted_timecourse(1:10)) > 0.1 | ...                                       % unreliable growth start
                mean(diff(fitted_timecourse(end-5:end))) > 0.03 | ...                           % unrealible growth fit
                max(abs(diff(fitted_timecourse))) > 1                                           % poor quality curve
            ind_remove(end+1) = q;
            continue
        end
        
        colony_struct(q).processed_timecourse = log_timecourse;
        colony_struct(q).fitted_timecourse = fitted_timecourse';
        colony_struct(q).lag_time = x(3);
        colony_struct(q).growth_rate = x(2);
        colony_struct(q).time_to_thresh = time2(find(fitted_timecourse>=thresh,1,'first'));
        colony_struct(q).max_density = max(log_timecourse);
        colony_struct(q).growth_rate_manual = growth_manual(time, log_timecourse, 3, 0);
        
    end
    
    colony_struct_filtered = colony_struct;
    colony_struct_filtered(ind_remove) = [];
    disp("Loop " + p + "_1: Skipped " + length(ind_remove) + " out of " + q + " colonies.")
    conditions = unique([colony_struct_filtered.condition]');
    
    if plotting_flag
        
        for qq = 1:length(conditions)
            col_struct_temp = colony_struct_filtered(strcmp([colony_struct_filtered.condition]',conditions(qq)));
            subplot(num_row,num_col,count), hold on
            plot(time,[col_struct_temp.processed_timecourse], clr(qq,1), 'linewidth', 1)
            plot(time2,[col_struct_temp.fitted_timecourse],clr(qq,2),'linewidth',1)
            title(strrep(string(extractBefore(s(p).name,"_rep")),"_","-"))
        end
    end
    
    if save_data_flag
        save(write_filename,'colony_struct_filtered','metadata','spot_metadata')
    end
    count=count+1;
    
end

% p.fontsize = 36;
