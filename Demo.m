clear all;
close all;
clc;

%% Parameter setting

% Default settings for LARK (object detection)
conf.Wsize = 5; % LARK window size
conf.h = 2; % smoothing parameter for LARK
conf.alpha = .13; % LARK sensitivity parameter
conf.colormode = 0; % 0: gray scale, 1: color
conf.interval = 3; % compute Covariance matrix at 3 pixel apart 
alpha = 0.99; % confidence level for significance testing

% Default settings for saliency detection
conf1.Wsize = 3; % LARK window size
conf1.alpha = .42; % LARK sensitivity parameter
conf1.h = 0.2; % smoothing parameter for LARK
conf1.L = 5; % # of LARK in the feature matrix
conf1.N = 3; % size of a center + surrounding region for computing self-resemblance
conf1.sigma = 0.07; % fall-off parameter for self-resemblance
conf1.colormode = 0; % 0: gray scale, 1: color
conf1.interval = 1; % compute Covariance matrix at 1 pixel apart 
conf1.block = [16 16];% block size to reduce search space
conf1.thres = 0.3; % threshold for saliency

%%
query = imread(['Faces/query1.jpg']);

% Compute LARKs from Query
[Q,W_Q] = CompLARK(query,conf);

for k = 1
    if k < 10
        FN = ['Faces/target_0' num2str(k) '.jpg'];
    else
        FN = ['Faces/target_' num2str(k) '.jpg'];
    end
    target = imread(FN);
    
    % Compute saliency map from Target and reduce search space
    tic;
    smap = ComputeSaliencyMap(target,[64 64],conf1); % Resize input images to [64 64]
    smap = imresize(smap,1/conf.interval);
    [block,flag,S,E] = Proto_Object(smap,conf1.block,conf1.thres);
    disp(['Saliency: ' num2str(toc) ' sec']);    
    
    % Compute LARKs from Target
    tic;    
    [T,W_T] = CompLARK(imresize(target,1,'lanczos3'),conf);
    disp(['LARK_target: ' num2str(toc) ' sec']);   
    
    figure(1),
    subplot(1,2,1),sc(cat(3,smap,double(T)),'prob_jet');
    subplot(1,2,2),sc(cat(3,block,double(T)),'prob_jet');

     % Obtain PCA reduced features from Q and T
    tic;
    [F_Q,F_T] = PCAfeature(Q,T,conf.Wsize,W_Q,W_T,0);
    disp(['PCA: ' num2str(toc) ' sec']);
    
    % Multi-scale search with coarse-to-fine search.
    MultiScaleSearch_CoarseToFine(Q,T,F_Q,F_T,S,E,flag,alpha,target,conf1.block);
    %     clear F_Q F_T W_T T
end

