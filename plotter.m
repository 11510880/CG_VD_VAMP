figure
[Kit,num_band,num_scale]= size(hist.true_err_Cthr);

for j = 1:num_scale
    for k= 1:num_band
        plot(hist.true_err_C(:,k,j));
        hold on
    end
end


% % ||r2k-c0||^2 /N VS ||r2k,i-c0,i||^2/ Ni in each band
% figure
% for j = 1:num_scale
%     for k= 1:num_band
%         plot(hist.true_err_r_2k(:,k,j));
%         hold on
%     end
% end
% ||c1k_hat-c0||^2 /N VS ||c1k,i_hat-c0,i||^2/ Ni in each band
figure
for j = 1:num_scale
    for k= 1:num_band
        plot(hist.true_err_Cthr(:,k,j));
        hold on
    end
end