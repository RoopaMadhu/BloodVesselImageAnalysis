function [] = initializeSeedsAndMask_2R(list,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% 03/02/2018 Modified by Roopa Madhu


% record original directory (and return to it at the end)
od = cd;

[spath, ~, ~] = fileparts(list{t});


% move to directory for saving seeds and mask
cd(spath);
cd('SegmentationData');
cframefoldername = sprintf('frame%04d',t);
cd(cframefoldername);

if exist('seeds.mat') == 2
    %image = temporalfilter(list,t,10);
    image = imread(list{t});
    
    loadseeds = load('seeds.mat');
    seeds = loadseeds.seeds;
    loadmask = load('mask.mat');
    mask = loadmask.mask;
else
    % load image
    %image = temporalfilter(list,t,8);
    image = imread(list{t});

%     % generate seeds to start with
    [seeds] = im2seeds(image,12,0);
    mask = false(size(seeds));
end

% make modifications
[seeds,mask] = modifySeedsAndMaskV2_RM(image,seeds,mask);


% save to current directory
save('seeds','seeds')
save('mask','mask')


% switch back to original directory
cd(od)
end

