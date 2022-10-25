% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quantify_images_demo.m
%   AUTHOR: AJL
%   DATE: 2021_01_06
%   UPDATE: 2022_05_04
%   DESCRIPTION: Quantify colony size in scanner imaging data.
%
%   Use this script to create time courses of colony growth for scanner
%   image data. This script uses classifiers created using the train_bag
%   scripts.
%
%   Between experiments, some changes may be needed to the variables below;
%   to set the experiment, change the 'plasmid', 'strain' and 'exp_date'
%   variables to read in the appropriate xlsx file. The processing
%   parameters may need tuning between experiments/dates, specifically
%   mask_time_scale. It is also generally useful to retrain
%   the models if performance seems suspect, especially in the case of
%   image noise classification.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1. Enter 1 into all rows from an experiment in the column single_analyze_flag in the
%%%%% images_master.xlsx file that you want to analyze, and *make sure*
%%%%% that all other rows contain a 0. Only rows with '1' will be analyzed.
%%%%% 2. Click Run at the top
%%%%% 3. Provide input to MATLAB to improve colony identification
%%%%% algorithm:
%%%%%         - Would you like to adjust the colony predictions?:
%%%%%                       Enter y or n
%%%%%         - Please enter the index/indices of the incorrect valid colonies:
%%%%%                       Examine the figure labeled "Predicted VALID
%%%%%                       Colonies" and look for abnormally shaped
%%%%%                       colonies or ones that are bleeding off the
%%%%%                       edge. If any are identified, input numbers that
%%%%%                       should be changed from VALID to INVALID by
%%%%%                       listing them in brackets, separated by commas,
%%%%%                       eg. [4,10,11,20,3]. If NONE are identified, hit
%%%%%                       ENTER
%%%%%         - Please enter the index/indices of the incorrect invalid colonies:
%%%%%                       Examine the figure labeled "Predicted INVALID
%%%%%                       Colonies" and look for colonies that are
%%%%%                       incorrectly assigned as invalid (eg ones that
%%%%%                       look like real colonies with a nice shaped
%%%%%                       circle and minimal edge bleeding). Enter
%%%%%                       numbers as described above
%%%%%         - Would you like to add the colonies to the library?:
%%%%%                       For datasets you feel REALLY good about, enter
%%%%%                       'y', otherwise, enter 'no' (without quotes)
%%%%%         - Enter image library path (leave blank for default):
%%%%%                       Leave blank and enter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% define current experiment details
home_path_flag = 1;                 % change to 0 for running off the hard drive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% define dataset
filename = "images_master.xlsx";
master_data_table = readtable(filename);
master_data_table(master_data_table.single_analyze_flag == 0,:) = [];
master_data_table.Experiment_Plasmid_replicate = string(master_data_table.Date) + ...
    string(master_data_table.Experiment) + ...
    string(master_data_table.plasmid) + ...
    string(master_data_table.replicate);
experiments = unique(master_data_table.Experiment_Plasmid_replicate);
saving_dir = pwd + "/pre_processed_data/";
if ~isfolder(saving_dir), mkdir(saving_dir), end
CONDITION_TO_USE = "T_type";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% define global processing parameters %%%%%%%%%%%%%%%%%%%%%%%
int = 1;                                      % interval for timepoints to process
colony_size_min = 30;                         % cutoff for pixels in a valid colony
colony_size_max = 3000;                       % cutoff for prixels in a valid colony
colony_size_filt = 30;                        % cutoff for pixels in a valid colony
mask_time_scale = [0.3 0.5 0.6 0.7 0.8 0.9];  % which timepoint to take as mask
dilations = 0:2:10;                           % range of dilations to use
use_colony_classifier = 1;                    % turn colony classifier on/off
bb_padding = 15;                              % padding for the bounding box, in pixels
contrast_cutoff = 15;                         % filter images that have small pixel intensity range
st = 1.5;                                     % filter threshold
colony_ceil = 5000;                           % max # of colonies in an image
colony_validation = 1;                        % 0/1 to perform colony classifier validation-feedback
colony_score_cutoff = 5;                      % cutoff to include ambiguous colonies; lower=stricter, -inf to disable
noise_score_cutoff = 5;                       % cutoff to include ambiguous noise images
save_data_flag = 1;                           % 1 for yes, 0 for no
save_checkpoint = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for EXPR = 1:length(experiments)
    disp("Currently analyzing exp #" + EXPR + " out of " + length(experiments) + ": " + experiments(EXPR))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% define current experiment details
    disp("Setting up parameters and loading dataset...")
    data_ind = strcmp(master_data_table.Experiment_Plasmid_replicate,experiments(EXPR));
    data_table = master_data_table(data_ind,:);
    image_type = data_table.image_type{1};
    plasmid = data_table.plasmid{1};
    exp_date = data_table.Date{1};
    exp_id = data_table.Experiment{1};
    INC = data_table.INC{1};
    drug = data_table.drug{1};
    transconjugants = data_table.transconjugants(1);
    replicate = data_table.replicate(1);
    metadata = {exp_date,exp_id,plasmid,INC,drug, transconjugants, replicate};
    im_size = -1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% define path
    pth = data_table.path{1};
    
    s = dir(pth);
    s = s(contains({s.name}',image_type) & ~contains({s.name}','._'));
    
    % sort images by increasing number
    snames = {s.name}';
    snames = split(snames,'-');
    snums = split(snames(:,end),'.');
    [hh s_ind] = sort(double(string(snums(:,1))));
    s = s(s_ind,:);
    
    if length(snums) > length(unique(snums))
        s = dir(pth);
        s = s(contains({s.name}',image_type) & ~contains({s.name}','._'));
    end
    
    % redefine ceiling for extra large images
    if s(end).bytes > 1000000
        colony_ceil = 9000;
    end
    
    % load classifiers
    load('noise_classifier-31-Jan-2021')
    load('colony_classifier-28-Jun-2021')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% define different condition types
    conditions = unique(data_table.(CONDITION_TO_USE));
    num_conditions = length(conditions);
    reps = zeros(num_conditions,1);
    centers = zeros(0,2);
    radii = zeros(0,1);
    condition_idx = zeros(0,1);
    colony_count = zeros(0,1);
    
    % obtain centers and radii from images_master.xlsx
    for k = 1:num_conditions
        
        curr_condition = conditions{k};
        idx = strcmp(data_table.T_type, curr_condition);
        reps(k) = sum(idx);
        centers = [centers;[data_table.center_x(idx),data_table.center_y(idx)]];
        radii = [radii; data_table.radius(idx)];
        colony_count = [colony_count; data_table.number_colonies(idx)];
        condition_idx = [condition_idx;repmat(string(curr_condition),reps(k),1)];
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% determine optimal mask for current dataset
    disp("Finding optimal mask and dilation parameters...")
    colony_limit = 500;
    mask_conditions = cell(1,num_conditions);
    mask_end = cell(1,num_conditions);
    rp_conditions = cell(1,num_conditions);
    col_num = zeros(1,num_conditions);
    col_area = zeros(1,num_conditions);
    mask_timepoint = zeros(1,num_conditions);
    
    for rr = mask_time_scale
        
        % read in final frame and adjust
        curr_frame = [pth,s(round(length(s)*rr)).name];
        im = rgb2gray(imread(curr_frame));
        
        % sometimes images can be 1 pixel off- fix size to avoid this
        if all(im_size == -1)
            im_size = size(im);
        elseif any(size(im) ~= im_size)
            im = imresize(im, im_size);
        end
        
        % loop through all conditions for current timepoint
        for k = 1:num_conditions
            
            % get centers and radii for current condition
            curr_centers = centers(strcmp(condition_idx, conditions{k}),:);
            curr_radii = radii(strcmp(condition_idx, conditions{k}),:);
            
            % generate an empty mask
            mask_curr = zeros(size(im,1),size(im,2));
            
            % loop through all spots for given condition and generate mask
            for r = 1:reps(k)
                [msh_x, msh_y] = meshgrid(1:size(im,1),1:size(im,2));
                mask = ((msh_y - curr_centers(r,1)).^2 + (msh_x - curr_centers(r,2)).^2 < curr_radii(r).^2)';
                roi_mean = mean(im(mask));
                im2 = uint8(mask).*im;
                im2(~mask) = roi_mean;
                im3= ~imbinarize(im2);
                if sum(im3(:))/numel(im3) > 0.8
                    im3= ~imbinarize(im2,'adaptive','ForegroundPolarity','dark');
                end
                mask_curr(logical(im3)) = 1;
            end
            
            % get binary blobs and filter based on size
            rp_curr = regionprops(logical(mask_curr),'Circularity','MajoraxisLength','MinoraxisLength','PixelIdxList','BoundingBox','Area','Centroid');
            rp_filtered = rp_curr([rp_curr.Area] > colony_size_filt);
            area_curr = [rp_filtered.Area]';
            
            % mask conditions
            if col_num(k) > 20
                num = floor(col_num(k)/3);
            else
                num = col_num(k) - 1;
            end
            
            if (isempty(mask_conditions{k}) & any(area_curr > 25000))
                mask_flag = 0;
            else
                mask_flag = 1;
            end
            
            % check to update mask
            if  mask_flag & ((sum(area_curr < colony_ceil) / length(area_curr) * 100) > 70 && ...
                    (size(area_curr,1) >= col_num(k) - 1) && ...
                    (size(rp_curr,1) < colony_limit)) || ...
                    (all(area_curr < colony_ceil) && ...
                    (col_area(k) > 5*sum(area_curr)))||...
                    (rr == mask_time_scale(end) && isempty(mask_conditions{k}))
                col_area(k) = sum(area_curr);
                col_num(k) = size(area_curr,1);
                mask_conditions{k} = mask_curr;
                rp_conditions{k} = rp_curr;
                mask_timepoint(k) = rr;
            end
            
            % store end-state masks for dilation
            mask_end{k} = mask_curr;
        end
    end
    
    % find best dilation parameter
    dilation_diffs = zeros(length(dilations),num_conditions);
    for j = 1:num_conditions
        cnt = 1;
        for k = dilations
            mask_dilated = imdilate(mask_conditions{j}, strel('disk',k));
            mask_dilated = mask_dilated(1:size(mask_end{j},1),1:size(mask_end{j},2));
            dilation_diffs(cnt,j) = abs(sum(reshape(mask_end{j} - mask_dilated, 1, [])));
            cnt = cnt+1;
        end
    end
    
    [~,dilution_idx] = min(dilation_diffs);
    dilation_use = dilations(dilution_idx);
    
    % clear out vars for memory consumption purposes
    clear mask_curr im2 im3 msh_x msh_y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% get valid colonies for each mask
    disp("Processing colonies using colony classifier...")
    final_pixels = cell(1, num_conditions);
    final_bb = cell(1, num_conditions);
    
    % loop through conditions and colonies
    for k = 1:num_conditions
        rp_curr = rp_conditions{k};
        num_cols_curr = zeros(1,length(rp_curr));
        invalid_masks = {};
        valid_masks = {};
        
        for r = 1:length(rp_curr)
            
            % get the current colony mask and dilate it
            mask = zeros(size(mask_conditions{k}));
            mask(rp_curr(r).PixelIdxList) = 1;
            mask = imdilate(mask, strel('disk',dilation_use(k)));
            
            % get the bounding box for the mask, update it
            rp_dilated = regionprops(mask, 'BoundingBox');
            BB = rp_dilated.BoundingBox;
            rp_curr(r).BoundingBox = BB;
            
            % reduce ROI down to only colony to evaluate if valid
            mask_curr = mask_conditions{k};
            
            try
                mask = mask_curr(round(BB(2)):round(BB(2))+BB(4),...
                    round(BB(1)):round(BB(1))+BB(3));
            catch
                mask = mask_curr(floor(BB(2)):floor(BB(2))+BB(4),...
                    floor(BB(1)):floor(BB(1))+BB(3));
            end
            mask = imdilate(mask, strel('disk',dilation_use(k)));
            
            % can use either manual or auto colony classification
            if ~use_colony_classifier
                im_sub = im(round(BB(2)):round(BB(2))+BB(4),...
                    round(BB(1)):round(BB(1))+BB(3));
                figure(1), subplot(1,2,1), imshow(mask)
                subplot(1,2,2), imshow(im_sub)
                set(gcf, 'position', [100 500 500 500])
                n = input("Enter 1 fo r valid colony, 0 for invalid: ");
                num_cols_curr(r) = n;
                close(1)
            else
                [idx, score] = predict(colony_net, 255*uint8(mask));
                
                if idx == 1 && score(2)/score(1) > colony_score_cutoff
                    num_cols_curr(r) = 0;
                    invalid_masks{end+1} = mask;
                else
                    num_cols_curr(r) = 1;
                    valid_masks{end+1} = mask;
                end
            end
        end
        
        % routine for allowing user input on colony predictions
        if use_colony_classifier && colony_validation
            disp("Getting user feedback on colony predictions...")
            
            % get initial valid/invalid colony counts & subplot sizes
            num_valid = sum(num_cols_curr);
            num_invalid = sum(~num_cols_curr);
            valid_subplot_size = floor(sqrt(num_valid)) + 1;
            invalid_subplot_size = floor(sqrt(num_invalid)) + 1;
            
            % plot the valid/invalid colonies
            figure(50), subplot(valid_subplot_size, valid_subplot_size, 1)
            set(gcf, "Name", "Predicted Valid Colonies")
            figure(51), subplot(invalid_subplot_size, invalid_subplot_size, 1)
            set(gcf, "Name", "Predicted Invalid Colonies")
            
            for j = 1:num_valid
                figure(50), subplot(valid_subplot_size, valid_subplot_size,j),...
                    imshow(valid_masks{j});
                title(num2str(j))
            end
            for j = 1:num_invalid
                figure(51), subplot(invalid_subplot_size, invalid_subplot_size, j), ...
                    imshow(invalid_masks{j});
                title(num2str(j))
            end
            
            % get input to determine whether to adjust
            n = input("Would you like to adjust the colony predictions?: ", 's');
            while ~ismember(lower(n), ["n","y","no","yes"])
                n = input("Please enter Y or N: ", 's');
            end
            
            if ismember(lower(n), ["y","yes"])
                
                idx_valid = input("Please enter the index/indices of the incorrect valid colonies: ");
                
                while (~isPositiveIntegerValuedNumeric(idx_valid) && ~isempty(idx_valid)) || any(idx_valid > num_valid)
                    disp("Error: input must be a vector of positive integers, scalar," +...
                        "or empty matrix, and cannot be greater than the number of valid colonies.");
                    idx_valid = input("Please enter the index/indices of the incorrect valid colonies: ");
                end
                
                idx_invalid = input("Please enter the index/indices of the incorrect invalid colonies: ");
                
                while (~isPositiveIntegerValuedNumeric(idx_valid) && ~isempty(idx_valid)) || any(idx_valid > num_valid)
                    disp("Error: input must be a vector of positive integers, scalar," +...
                        "or empty matrix, and cannot be greater than the number of invalid colonies.");
                    idx_invalid = input("Please enter the index/indices of the incorrect invalid colonies: ");
                end
                
                % adjust the booleans appropriately
                original_valid = find(num_cols_curr);
                original_invalid = find(~num_cols_curr);
                num_cols_curr(original_valid(idx_valid)) = 0;
                num_cols_curr(original_invalid(idx_invalid)) = 1;
                
                
                real_valid_masks = valid_masks;
                real_valid_masks(idx_valid) = [];
                real_valid_masks = [real_valid_masks,invalid_masks(idx_invalid)];
                
                real_invalid_masks = invalid_masks;
                real_invalid_masks(idx_invalid) = [];
                real_invalid_masks = [real_invalid_masks,valid_masks(idx_valid)];
            end
            
            close(50); close(51)
            
            % get input on whether to save the images back to the
            % classifier library or not
            n = input("Would you like to add the colonies to the library?: ", 's');
            while ~ismember(lower(n), ["n","y","no","yes"])
                n = input("Please enter Y or N: ", 's');
            end
            
            if ismember(lower(n), ["y","yes"])
                skip_write = 0;
                dest_pth = pwd + "/img_learning/colony/";
                n = input("Enter image library path (leave blank for default): ", 's');
                
                if ~isempty(n)
                    dest_path = n;
                    if ~isfolder(dest_path)
                        disp("Invalid path. Skipping image write-back...")
                        skip_write = 1;
                    end
                else
                    disp("Using default path of "+dest_pth)
                end
                
                if ~skip_write
                    n = dir(dest_pth+"0");
                    n(contains({n.name},'Icon')) = [];
                    n_count = max(cellfun(@(x) str2double(string(regexp(x,"([0-9]+)","match"))),{n(4:end).name}));
                    y = dir(dest_pth+"1");
                    y(contains({y.name},'Icon')) = [];
                    y_count = max(cellfun(@(x) str2double(string(regexp(x,"([0-9]+)","match"))),{y(4:end).name}));
                    
                    for j = 1:length(real_valid_masks)
                        im_sub = real_valid_masks{j};
                        imwrite(im_sub,dest_pth+"1/"+y_count+".jpg","jpg")
                        y_count = y_count+1;
                    end
                    for j = 1:length(real_invalid_masks)
                        im_sub = real_invalid_masks{j};
                        imwrite(im_sub,dest_pth+"0/"+n_count+".jpg","jpg")
                        n_count = n_count+1;
                    end
                end
            end
        end
        
        % filter the colonies based on valid/invalid booleans
        rp_curr = rp_curr(num_cols_curr == 1);
        rp_conditions{k} = rp_curr;
        final_pixels{k} = {rp_curr.PixelIdxList}';
        final_bb{k} = {rp_curr.BoundingBox}';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% initiate collection matxs and loop through time
    if save_checkpoint
        save("run_" + date + ".mat", "num_conditions","final_pixels","s","pth",...
            "final_bb","bb_padding","rp_conditions","centers","radii",...
            "colony_count","data_table","conditions","metadata")
    end
    
    disp("Collecting colony timecourses...")
    filtered_timecourses = cell(1,num_conditions);
    for k = 1:num_conditions
        filtered_timecourses{k} = zeros(size(s,1),length(final_pixels{k}));
    end
    time = ((1:size(s,1))*15 - 15) ./ 60;
    
    for q = 1:int:length(time)
        disp("On loop " + q + " of " + size(s,1))
        
        % load in current time frame
        curr_frame = [pth,s(q).name];
        im = rgb2gray(imread(curr_frame));
        
        % loop through conditions
        for k = 1:num_conditions
            
            % get current params
            curr_pixels = final_pixels{k};
            curr_bbs = final_bb{k};
            curr_timecourse = filtered_timecourses{k};
            
            % loop through colonies
            for r = 1:length(curr_pixels)
                BB = curr_bbs{r};
                im_raw = im(round(BB(2)-bb_padding):round(BB(2)-bb_padding)+BB(4)+2*bb_padding,...
                    round(BB(1)-bb_padding):round(BB(1)-bb_padding)+BB(3)+2*bb_padding);
                
                % create circle mask to apply to colony
                circle_mask = zeros(size(im_raw));
                c_x = round(size(circle_mask,1)/2); c_y = round(size(circle_mask,2)/2);
                radius = min(c_x,c_y);
                [msh_x, msh_y] = meshgrid(1:size(circle_mask,1),1:size(circle_mask,2));
                mask = ((msh_y - c_y).^2 + (msh_x - c_x).^2 < radius.^2)';
                
                % binarize the image
                im_new = ~imbinarize(im_raw);
                [p, score] = predict(noise_net, 255*uint8(im_new));
                score_ratio = score(2)/score(1);
                
                if max(im_raw(:))-min(im_raw(:)) < contrast_cutoff || (p == 1 && score_ratio > noise_score_cutoff)
                    curr_timecourse(q,r) = 0;
                else
                    im_new = mask & bwareaopen(logical(im_new),25);
                    curr_timecourse(q,r) = sum(im_new(:));
                end
            end
            filtered_timecourses{k} = curr_timecourse;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% form the final structures to be saved
    
    % example of pulling all adapted timecourses into a matrix:
    % AD=[colony_struct(strcmp([colony_struct.condition],"adapted")).timecourse];
    %
    % or, pulling all timecourses excluding a given list of spots:
    % col = [colony_struct(~ismember([colony_struct.spot_id],[1,2,3])).timecourse];
    
    disp("Creating and saving final structures...")
    
    % create colony structure
    fields = fieldnames(rp_conditions{1})';
    fields = [fields,{'condition','timecourse','time','spot_id'}];
    fields{2,1} = {};
    colony_struct = struct(fields{:});
    for k = 1:num_conditions
        for j = 1:length(rp_conditions{k})
            curr_struct = rp_conditions{k}(j);
            curr_struct.condition = string(conditions{k});
            curr_struct.timecourse = filtered_timecourses{k}(:,j);
            curr_struct.time = time;
            curr_struct.spot_id = find((curr_struct.Centroid(1)-centers(:,1)).^2 + ...
                (curr_struct.Centroid(2)-centers(:,2)).^2 <= radii.^2, 1, 'first');
            colony_struct(end+1) = curr_struct;
        end
    end
    % create metadata table for spots
    [unique_count, unique_vals] = hist([colony_struct.spot_id], unique([colony_struct.spot_id]));
    final_count = zeros(1,length(centers));
    final_count(ismember(1:length(centers),unique_vals)) = unique_count;
    spot_metadata = [(1:length(centers))', final_count', centers, radii, colony_count];
    spot_metadata = array2table(spot_metadata, 'VariableNames', ...
        {'spot_id','colony_count','center_x','center_y','radius','real_colony_count'});
    spot_metadata.condition = string(data_table.T_type);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% filter and save data %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if save_data_flag
        savename = saving_dir + plasmid + "_" + drug + "_" + exp_id + "_rep" + replicate + ".mat";
        save(savename,'spot_metadata','metadata','colony_struct')
    end
    
    figure; hold on, clr = 'brgkm';
    unique_conditions = unique([colony_struct.condition]');
    for q = 1:length(unique_conditions)
        tc=[colony_struct(strcmp([colony_struct.condition],unique_conditions(q))).timecourse];
        time = colony_struct.time;
        tc = tc(1:int:end,:);
        tcFINAL = tc(find(sum(tc>0,2) == max(sum(tc>0,2)),1,'last'),:);
        ind_tc = tcFINAL > mean(tcFINAL) + st*std(tcFINAL) | tcFINAL < mean(tcFINAL) - st*std(tcFINAL);
        tc = tc(:,~ind_tc);
        plotting_images(time,tc,int,metadata,0,clr(q)), xlabel('Hours'), ylabel('Pixels')
    end
    
    
end
