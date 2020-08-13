clear all; 
addpath(genpath(pwd))
%% tester for the three algorithms, which should be run before plotter.m
%% ground truth
x0 = phantom('modified shepp-logan', 512);
[nx, ny] = size(x0);
y0 = fftnc(x0);
rng('default');
rng(777);
% hist1: CG-VAMP hist2: CG-VDAMP hist3: CG-VD-VAMP
hist1_cell = cell(1,3); % since there are 12x, 8x and 4x
hist2_cell = cell(1,3);
hist3_cell = cell(1,3);
cor1_cell= cell(1,3);
cor2_cell = cell(1,3);
cor3_cell = cell(1,3);
%% generate sampling mask
% notice we run three times for different undersampling factors 
target_delta = [1/12;1/8;1/4];
%prob_map = genPDF([nx,ny], 8, target_delta, 2, 0,0);
 %target_delta = 1/2;
 %prob_map = ones(nx,ny).*target_delta;
 timer1 =zeros(3,1); % cg_vamp
 timer2 =zeros(3,1); % cg_vd_vamp
 timer3 =zeros(3,1); % vdamp
 timer1_cg =cell(3,1); % timer for computing total cg running time of cg-vamp
 timer2_cg =cell(3,1);% timer for computing total cg running time of cg-vd-vamp
for i=1:3
    prob_map = genPDF([nx,ny], 8, target_delta(i), 2, 0,0);
    mask = binornd(1, prob_map, nx,ny);
    delta=mean(mask(:));
    %% Generate data
    SNR = 40;
    var0 = mean(abs(y0(:)).^2)/(10^(SNR/10));  
    noise = normrnd(0, sqrt(var0), nx,ny)./sqrt(2) + 1i* normrnd(0, sqrt(var0), nx,ny)./sqrt(2);
    dcoil = mask.*(fftnc(x0) + noise);

    %% CG-VD-VAMP options
    opts.maxIter = 30;
    opts.maxTime = 100;
    opts.verbose = 1;
    opts.scales = 4;
    opts.SURE =1;
    opts.lambda=35;
    opts.saveHist = 1;
    opts.denoiserDiv = 0; % VDAMP-alpha
    scales = 4;
    % cor ij , i = 1,2 j =1,2,3
    %% Run VDAMP 
    t3 = tic;
    [x_hat3,hist3, cor13, cor23] = VDAMP_O(dcoil, mask, prob_map, var0, x0, opts);
    timer3(i) = toc(t3);
    hist3_cell{1,i} =hist3; 
    cor3_cell{1,i}=[cor13,cor23];
    %[x_hat, hist,cor1,cor2]  = VDAMP_O(dcoil,  mask, prob_map, var0, x0, opts);
    %% Run CG-VAMP
    t1 = tic;
    [x_hat1, cor11, cor21,hist1] = CG_VAMP_MRI(dcoil, mask, prob_map, var0, x0, opts);
    timer1(i) = toc(t1);
    hist1_cell{1,i} =hist1; 
    cor1_cell{1,i}=[cor11,cor21];
    timer1_cg{i}= hist1.cg_timer;
    %% Run CG-VD-VAMP
    t2 = tic;
    [x_hat2, cor12, cor22,hist2] = CG_VD_VAMP(dcoil, mask, prob_map, var0, x0, opts);
    timer2(i) = toc(t2);
    hist2_cell{1,i}=hist2;
    cor2_cell{1,i}=[cor12,cor22];
    timer2_cg{i}=hist2.cg_timer;
end






