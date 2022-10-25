close all, clear all, clc

% name of experiment folder: replace with new experiment
exp_folder = "7.6.014-15_pKpQILD2_carb";

% load in data
pth = "path/to/image/files/" + exp_folder + "/";
s = dir(pth);
s = s(contains({s.name}','.jpg'));

% read in final frame and adjust
curr_frame = pth + s(end-1).name;
im = rgb2gray(imread(curr_frame));

% draw figure
figure; imshow(im)
circle_info = drawcircle('Center',[1000,1000],'Radius',96);

% Manually place circle on each Adapted and Denovo spot and enter 'h' into
% the command line to copy and paste each center / radius into the excel
% file