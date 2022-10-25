% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post_processing_tables.m
%   AUTHOR: AJL
%   DATE: 2021_01_06
%   UPDATE: 2021_09_02
%   DESCRIPTION: Process timecourses of scanner acquisition cost data
%
%   Use this script to clean up time courses of colony growth for scanner
%   image data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Define plasmid_drugs to analyze and saving directory
data_pth = pwd + "/pre_processed_data/filtered_structs/";
files = dir(data_pth);
files = files(contains({files.name}','.mat'),:);
plasmid_drugs = unique(extractBefore( {files.name}','_7'));
plasmid_free_flag = 0;
plasmid_drugs = plasmid_drugs(~contains(plasmid_drugs,'plasmidfree'));
saving_dir = pwd + "/post_processed_data/";
if ~isdir(saving_dir), mkdir(saving_dir), end


%%%% Flags for analysis
plotting_flag = 1;
save_data_flag = 1;
subplot_flag = 1;
separate_flag = 0;
clr = ['bc';'rm'];
overwrite_data_flag = 1;

if overwrite_data_flag
    delete(saving_dir + "*")
end

figure, hold on
st_thresh = 2;
%%%% Loop through all plasmid_drugs
count = 1;
for p = 1:length(plasmid_drugs)
    
    plasmid = plasmid_drugs(p);
    disp("Loop " + p + ": Currently analyzing plasmid " + plasmid)
    s_all = dir(data_pth);
    s_all = s_all(contains({s_all.name}',plasmid));
    
    t = [];
    file_name = saving_dir + plasmid + ".xlsx";
    
    for i = 1:size(s_all,1)
        
        % load current data file
        load(data_pth + s_all(i).name);
        conditions = unique([colony_struct_filtered.condition]');
        
        
        for c = 1:length(conditions)
            
            colony_struct_temp = colony_struct_filtered(strcmp([colony_struct_filtered.condition],conditions(c)));
            time=colony_struct_temp.time;
            
            % inialize variables
            growthrate = [colony_struct_temp.growth_rate];
            manual_growthrate = [colony_struct_temp.growth_rate_manual];
            lagtime = [colony_struct_temp.lag_time];
            time_to_thresh = [colony_struct_temp.time_to_thresh];
            max_size = [colony_struct_temp.max_density];
            
            log_timecourse = [colony_struct_temp.processed_timecourse];
            fitted_timecourse = [colony_struct_temp.fitted_timecourse];
            
            meta_table = array2table(repmat([metadata,conditions(c)],length(growthrate),1),'VariableNames',{'Date','ID','Plasmid','INC','drug', 'transonjugants', 'replicate', 'type'});
            data_table = [meta_table,array2table([manual_growthrate', growthrate',lagtime',time_to_thresh',max_size'],'VariableNames',{'gr_MANUAL', 'gr','lt','ttt','max'})];
            %
            filter_ind = data_table.max > (mean(data_table.max) + st_thresh * std(data_table.max)) | ...
                data_table.max < (mean(data_table.max) - st_thresh * std(data_table.max)) | ...
                data_table.lt > (mean(data_table.lt) + st_thresh * std(data_table.lt)) | ...
                data_table.lt < (mean(data_table.lt) - st_thresh * std(data_table.lt));
            %
             if plotting_flag
                if subplot_flag
                    subplot(8,ceil(length(plasmid_drugs)/8),count), hold on
                end
                
                time2 = linspace(time(1),time(end),size(fitted_timecourse,1));
                if separate_flag
                    subplot(1,2,c), hold on
                    plot(time,log_timecourse(:,~filter_ind), clr(c,1), 'linewidth', 1)
                    plot(time2,fitted_timecourse(:,~filter_ind),clr(c,2),'linewidth',1)
                else
                    plot(time,log_timecourse(:,~filter_ind), clr(c,1), 'linewidth', 1)
                    plot(time2,fitted_timecourse(:,~filter_ind),clr(c,2),'linewidth',1)
                end
                text(.02,0.95,plasmid,'fontSize',12,'Units','Normalized','Interpreter','none')
            end
            t = [t;data_table];
        end
        
        % write new file or append existing file with new data
        if save_data_flag
            writetable(t,file_name)
        end
        
    end
    count = count+1;
end