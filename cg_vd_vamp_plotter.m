close all
%% plotter adapted from vdamp @ Chengpin Luo 2020

% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2020  Charles Millard
% Copyright (C) 2020  Siemens Healthineers


c0 = multiscaleDecomp(x0, opts.scales);


%% BASIC SETTINGS FIGURE 3.1
figure
imshow(abs(x0));
title('original image');
figure
imshow(abs(abs(hist2.r1(:,:,1))));
title('unbiased estimate');
figure
imshow(mask);
title('Polynomial mask')
figure
imshow(prob_map)
title('Probability map')
%% view reconstructed image and its error FIGURE 3.2

figure
imagesc(abs(x_hat1), [0 max(x0(:))]);
colormap(gca, 'gray')
axis image off;
title('CG-VAMP RECONSTRUCTED IMAGE')


figure
imagesc(abs(x_hat1-x0)); 
axis image off;
colormap(gca, 'jet')
colorbar off;
title('CG-VAMP Reconstruction Error')


figure
imagesc(abs(x_hat2), [0 max(x0(:))]);
colormap(gca, 'gray')
axis image off;
title('CG-VD-VAMP RECONSTRUCTED IMAGE')


figure
imagesc(abs(x_hat2-x0)); 
axis image off;
colormap(gca, 'jet')
colorbar off;
title('CG-VD-VAMP Reconstruction Error')

figure
imagesc(abs(x_hat3), [0 max(x0(:))]);
colormap(gca, 'gray')
axis image off;
title('VDAMP RECONSTRUCTED IMAGE')


figure
imagesc(abs(x_hat3-x0)); 
axis image off;
colormap(gca, 'jet')
colorbar off;
title('VDAMP Reconstruction Error')

%% FIGURE 3.3
figure
% CG-VAMP
plot(hist1_cell{1,1}.x_nmse,'LineWidth',1.5,'Marker','x','Color','#0072BD') %12x
hold on 
plot(hist1_cell{1,2}.x_nmse,'LineWidth',1.5,'Marker','x','Color','#D95319') %8x
hold on
plot(hist1_cell{1,3}.x_nmse,'LineWidth',1.5,'Marker','x','Color','#EDB120') %4x
hold on
% CG-VD-VAMP
plot(hist2_cell{1,1}.x_nmse,'LineWidth',1.5,'Marker','o','Color','#0072BD')
hold on 
plot(hist2_cell{1,2}.x_nmse,'LineWidth',1.5,'Marker','o','Color','#D95319')
hold on
plot(hist2_cell{1,3}.x_nmse,'LineWidth',1.5,'Marker','o','Color','#EDB120')
hold on
% VDAMP
plot(hist3_cell{1,1}.x_nmse,'LineWidth',1.5,'Marker','^','Color','#0072BD')
hold on 
plot(hist3_cell{1,2}.x_nmse,'LineWidth',1.5,'Marker','^','Color','#D95319')
hold on
plot(hist3_cell{1,3}.x_nmse,'LineWidth',1.5,'Marker','^','Color','#EDB120')
xlabel('Iterations','FontSize',20);
ylabel('NMSE(dB)','FontSize',20);
legend('CG-VAMP,12x','CG-VAMP,8x','CG-VAMP,4x','CG-VD-VAMP,12x','CG-VD-VAMP,8x',...
'CG-VD-VAMP,4x','VDAMP,12x','VDAMP,8x','VDAMP,4x','FontSize',15,'NumColumns',3);
set(gca,'FontSize',20)


%% PLOT IN THE WAVELET DOMAIN 
%% FIGURE 3.4
figure 
c0_py = pyramid(c0); %pyramid structure of c0
it_choice = [2, 3, 4]; % which iterations to look at
imshow(abs(c0_py));

r_1k_py = pyramid(multiscaleDecomp(hist1.r1(:,:,1),opts.scales));
figure,imshow(abs(r_1k_py));
diff = abs(r_1k_py-c0_py);
figure('Name', 'Difference between r10 and c0');
imagesc(abs(r_1k_py-c0_py),[0,1]); 
axis image off;
colormap(gca, 'jet')
colorbar off;

%% FIGURE 3.5
for iter = 1:numel(it_choice)  
    subplot(3,numel(it_choice),iter);
    ii = it_choice(iter);
    %figure
    r_1k = multiscaleDecomp(hist1.r1(:,:,ii),opts.scales);
    r_1k_py = pyramid(r_1k);
    diff = abs(r_1k_py-c0_py);
    
    if iter == 1
        thr = 0.4*max(diff(:)); % choose clip
        diff(diff>thr)= thr;
        diff = rescale(diff);
    else
        diff(diff>thr) = thr;
        diff = rescale(diff)*max(diff(:))/thr;
    end

    
    imagesc(diff, [0, 1]); colormap jet;
    ax = gca;
    axis(ax,'off')
    if iter ==1
       
       ylabel(ax,'CG-VAMP','FontSize',15);
       set(get(ax,'YLabel'),'Visible','on')
    end
    %title(['k=', num2str(ii-1)]);
    colorbar off;
    %axis image off;% image
    %set(gca,'FontName','times')
end

for iter = 1:numel(it_choice)  
    subplot(3,numel(it_choice),3+iter);
    ii = it_choice(iter);
    %figure
    r_1k = multiscaleDecomp(hist2.r1(:,:,ii),opts.scales);
    r_1k_py = pyramid(r_1k);
    diff = abs(r_1k_py-c0_py);
    
    if iter == 1
        thr = 0.4*max(diff(:)); % choose clip
        diff(diff>thr)= thr;
        diff = rescale(diff);
    else
        diff(diff>thr) = thr;
        diff = rescale(diff)*max(diff(:))/thr;
    end

    
    imagesc(diff, [0, 1]); colormap jet;
    ax = gca;
    axis(ax,'off')
    if iter ==1
       
       ylabel(ax,'CG-VD-VAMP','FontSize',15);
       set(get(ax,'YLabel'),'Visible','on')
    end
    %title(['k=', num2str(ii-1)]);
    colorbar off;
    %axis image off;% image
    
    %set(gca,'FontName','times')
end

for iter = 1:numel(it_choice)  
    subplot(3,numel(it_choice),6+iter);
    ii = it_choice(iter);
    %figure
    r_1k = multiscaleDecomp(hist3.r1(:,:,ii),opts.scales);
    r_1k_py = pyramid(r_1k);
    diff = abs(r_1k_py-c0_py);
    
    if iter == 1
        thr = 0.4*max(diff(:)); % choose clip
        diff(diff>thr)= thr;
        diff = rescale(diff);
    else
        diff(diff>thr) = thr;
        diff = rescale(diff)*max(diff(:))/thr;
    end

    
    imagesc(diff, [0, 1]); colormap jet;
    ax = gca;
    axis(ax,'off')
    if iter ==1
       ylabel(ax,'VDAMP','FontSize',15);
       set(get(ax,'YLabel'),'Visible','on')
    end
    xlabel(ax,['k=', num2str(ii-1)],'FontSize',15);
    set(get(ax,'XLabel'),'Visible','on')
    %title(['k=', num2str(ii-1)]);
    colorbar off;
    %axis image off% image off;
    %set(gca,'FontName','times')
end
%axis off

%% QQ plots FIGURE 3.6, 3.7, 3.8
it_choice = [2, 3, 4]; % which iterations to look at

div = 1; %normalisation
qqylim = [-1, 1];
figure('Name', 'QQ plots');
for iter = 1:numel(it_choice)
    subplot(3,numel(it_choice),iter);
    ii = it_choice(iter);
    r_1k = multiscaleDecomp(hist3.r1(:,:,ii),opts.scales);
    
    diff = (c0{1}.D - r_1k{1}.D)/div; 
    
    p = qqplot(real(diff(:)));
    set(p, 'Marker', '.');
    
    xh = xlabel('Standard normal quantiles');
    set(xh, 'FontName', 'times');
    set(gca,'FontSize',8)
    
    yh = ylabel('Quantiles of input sample');
    set(yh, 'FontName', 'times');
    
    axis square;
    title(['Diagonal, scale 1, k=', num2str(ii-1)], 'FontWeight','bold', 'FontName', 'Times new roman');
    ylim([-1 1]);
    xlim([-5, 5]);
end

for iter = 1:numel(it_choice) 
    subplot(3,numel(it_choice),iter+numel(it_choice));
    ii = it_choice(iter);
    r_1k = multiscaleDecomp(hist3.r1(:,:,ii),opts.scales);
    
    diff = (c0{2}.H - r_1k{2}.H)/div; 
    
    p = qqplot(imag(diff(:)));
    set(p, 'Marker', '.');
    
    xh = xlabel('Standard normal quantiles');
    set(xh, 'FontName', 'times');
    
    yh = ylabel('Quantiles of input sample');
    set(yh, 'FontName', 'times');
    axis square;
    set(gca,'FontSize',8)
    
    title(['Horizonal, scale 2, k=', num2str(ii-1)], 'FontWeight','bold', 'FontName', 'Times new roman');
    ylim([-1 1]);
    xlim([-5, 5]);
end

for iter = 1:numel(it_choice)
    subplot(3,numel(it_choice),iter+2*numel(it_choice));
    ii = it_choice(iter);
    r_1k = multiscaleDecomp(hist3.r1(:,:,ii),opts.scales);
    
    diff = (c0{4}.V - r_1k{4}.V)/div; 
    
    p = qqplot(real(diff(:)));
    set(p, 'Marker', '.');
    xh = xlabel('Standard normal quantiles');
    set(xh, 'FontName', 'times');
    
    yh = ylabel('Quantiles of input sample');
    set(yh, 'FontName', 'times');
    axis square;
    set(gca,'FontSize',8)
    
    title(['Vertical, scale 4, k=', num2str(ii-1)], 'FontWeight','bold', 'FontName', 'Times new roman');
    ylim([-1 1]);
    xlim([-5, 5]);
end


%% INPUT AND OUTPUT Correlations FIGURE 3.9
% 12x undersampling
% cg-vamp
figure
subplot(2,1,1)
for i=1:3
    plot(cor1_cell{i}(:,1),'LineWidth',1.5);
    hold on
end
xlabel('Iterations','FontSize',15);
ylabel('Input Correlation','FontSize',15);
title('CG-VAMP','FontSize',15);
legend('12x','8x','4x','FontSize',15);
subplot(2,1,2)
for i=1:3
    plot(cor1_cell{i}(:,2),'LineWidth',1.5);
    hold on
end
xlabel('Iterations','FontSize',15);
ylabel('Output Correlation','FontSize',15);
legend('12x','8x','4x','FontSize',15);


% cg-vd-vamp
figure
subplot(2,1,1)
for i=1:3
    plot(cor2_cell{i}(:,1),'LineWidth',1.5);
    hold on
end
xlabel('Iterations','FontSize',15);
ylabel('Input Correlation','FontSize',15);
title('CG-VD-VAMP','FontSize',15);
legend('12x','8x','4x','FontSize',15);
subplot(2,1,2)
for i=1:3
    plot(cor2_cell{i}(:,2),'LineWidth',1.5);
    hold on
end
xlabel('Iterations');
ylabel('Output Correlation');
legend('12x','8x','4x');

% vdamp
figure
subplot(2,1,1)
for i=1:3
    plot(cor3_cell{i}(:,1),'LineWidth',1.5);
    hold on
end
xlabel('Iterations','FontSize',15);
ylabel('Input Correlation','FontSize',15);
title('VDAMP','FontSize',15);
legend('12x','8x','4x','FontSize',15);
subplot(2,1,2)
for i=1:3
    plot(cor3_cell{i}(:,2),'LineWidth',1.5);
    hold on
end
xlabel('Iterations','FontSize',15);
ylabel('Output Correlation','FontSize',15);
legend('12x','8x','4x','FontSize',15);


