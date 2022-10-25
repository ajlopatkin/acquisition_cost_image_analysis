function timecourse = process_timecourse(timecourse, varargin)

%perform at most 5 loops of outlier removal
for q = 1:5
    % get outliers based on moving median
    [~, pts_rm] = rmoutliers(timecourse,'movmedian',10);

    if ~any(pts_rm)
        break
    end
    
    % if last point is marked as outlier, set it to the preceding value
    if pts_rm(end)
        timecourse(end) = timecourse(end-1);
        pts_rm(end) = 0;
    end

    % if first point is marked as outlier, set it to 0
    if pts_rm(1)
        timecourse(1) = 0;
        pts_rm(1) = 0;
    end

    % loop through each outlier and correct it
    for qq = find(pts_rm)'
        last_good_pt = find(~pts_rm(1:qq),1,'last');
        next_good_pt = find(~pts_rm(qq:end),1,'first')+qq-1;
        timecourse(pts_rm) = mean([timecourse(last_good_pt),timecourse(next_good_pt)]);
    end
end

% attempt to remove dips
for qq = 1:length(timecourse)

    t1 = timecourse(1:end-1) == 0;
    t2 = timecourse(2:end) == 0;
    ind_zero = find(t1.*t2 ~=0, 1,'last');
    if ~isempty(ind_zero) & ind_zero < 50
        timecourse(1:ind_zero) = 0;
    end
    
    dips = find(timecourse == 0);
    no_dips = find(dips > 30);
    fix = dips(no_dips);
    if fix > 0
        for oo = 1:size(fix, 1)
            timecourse(fix(oo, 1)) = timecourse(fix(oo, 1)-1);
        end
    end
    
end

% perform baseline checking (if we are given an average baseline)
if ~all(timecourse == 0) && nargin == 2
    avg_baseline = varargin{1};
    baseline = min(timecourse(timecourse>0));
    
    % if the timecourse baseline is greater than 1.5x our average baseline
    % across all timecourses, use the average; otherwise, use the baseline
    % from this timecourse
    if baseline >= avg_baseline*1.5
        timecourse(timecourse == 0) = avg_baseline;
    else
        timecourse(timecourse == 0) = baseline;
    end
end

if max(diff(timecourse)) > 500
    timecourse = abs(smooth(timecourse,'rlowess'));
end