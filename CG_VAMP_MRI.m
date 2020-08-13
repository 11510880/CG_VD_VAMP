function [x_hat, cor1, cor2,hist] = CG_VAMP_MRI(y, mask, prob_map, sigma_w, x0, opts)
    % CG-VAMP adapted from code shared by Nikolai@ Skuratovs, Nikolajs
    % The original code is from the paper "Upscaling Vector Approximate Message Passing"
    
    
    [nx, ny] = size(y);
    
    if (isfield(opts,'maxIter') && ~isempty(opts.maxIter))
        maxIter = opts.maxIter;
    else
        maxIter = 50;
    end
    
        
    if (isfield(opts,'maxTime') && ~isempty(opts.maxTime))
        maxTime = opts.maxTime;
    else
        maxTime = 60;
    end
    
    if (isfield(opts,'verbose') && ~isempty(opts.verbose))
        verbose = opts.verbose;
    else
        verbose = 0;
    end
    
    if (isfield(opts,'scales') && ~isempty(opts.scales))
        scales = opts.scales;
    else
        scales = 4;
    end
    
    if (isfield(opts,'SURE') && ~isempty(opts.SURE)) 
        sureFlag = opts.SURE;
    else
        sureFlag = 1;
    end
    
    if (isfield(opts,'denoiserDiv') && ~isempty(opts.denoiserDiv)) 
        denoiserDiv = opts.denoiserDiv;
    else
        denoiserDiv = 1;
    end
    
    if (isfield(opts,'stopDelta') && ~isempty(opts.stopDelta)) 
        stopDelta = opts.stopDelta;
    else
        stopDelta = 0;
    end
   
    hist.r1 = zeros(nx,ny,maxIter);
    hist.x = zeros(nx,ny,maxIter);
    hist.true_err_r_1k= zeros(maxIter, 1); % mean(||r1k-c0||.^2)
    hist.belief_std_r_1k = zeros(maxIter, 1);
    hist.true_err_c_1k_hat= zeros(maxIter, 1);
    hist.belief_std_c_1k_hat = zeros(maxIter, 1);% mean(||r1k_hat-c0||.^2)
    hist.true_err_r_2k= zeros(maxIter, 1);
    hist.belief_std_r_2k = zeros(maxIter, 1);% mean(||r2k-c0||.^2)
    hist.true_err_c_2k_hat= zeros(maxIter, 1);% mean(||r2k_hat-c0||.^2)
    hist.belief_std_c_2k_hat = zeros(maxIter, 1);
    hist.Cdiv = zeros(maxIter, 1);
    hist.RMSEyest = zeros(1, maxIter);
    hist.x_mse = zeros(1, maxIter);   
    hist.x_nmse = zeros(1, maxIter);
    hist.cg_timer=zeros(1,maxIter);
    lambda = 0;
    hist.c_1k_hat = zeros(nx,ny,maxIter);      
    c0 = multiscaleDecomp(x0, scales); % wavelet cofficients of original image in cell-form
    c0_vec =cell2vec(c0);
    inv_p = prob_map.^(-1);
    % unbaised initialisation
    m2 = inv_p.*mask.*y;
    r = ifftnc(m2);
    r_1k = multiscaleDecomp(r, scales);
    [r_1k_vec,l] = cell2vec(r_1k);
    h_t = r_1k_vec -c0_vec; %r1k = c0+h_t
    % use the oracle information to initialize sigma_1k
    sigma_1k = mean(h_t'*h_t);
    cor1 = zeros(maxIter,1);
    cor2 = zeros(maxIter,1);
    time_init = tic;
    for iter = 1:maxIter
        hist.r1(:,:,iter) = r;
        if verbose % band-wise mse of unthresholded image
            disp(['                 ITERATION ', num2str(iter)])
            fprintf( ...
                'mean(||r1k_time-x0||^2): true RMSE = %f\n', ...
                sqrt(mean(abs(r(:)-x0(:)).^2)));
            
                
                    fprintf( ...
                        'true RMSE = %f\tmessage std = %f\n', ...
                        sqrt(mean(abs(r_1k_vec(:)-c0_vec(:)).^2)), ...
                        sqrt(sigma_1k));
        end
        % thresholding
        if sureFlag
            [c_1k_hat_vec, df] = SUREsoft(r_1k_vec, sigma_1k);
            %[gz, df] = SUREsoft(z, V)
        else      
            [c_1k_hat_vec, df] = ComplexSoft(r_1k, sigma_1k, lambda);
        end
        % e = df{s}.(band_name) .* var{s}.(band_name);
        % err{s}.(band_name) = mean(e(:)) / 2;
        err = mean(df.*sigma_1k)/2;
        hist.c_1k_hat(:,:, iter) = pyramid(vec2cell(c_1k_hat_vec,l,scales));    
        % check stopping criterion
        if stopDelta >0 && iter> 1
            c_1k_hat_new_vec = c_1k_hat_vec;
            dk = norm(c_1k_hat_old_vec - c_1k_hat_new_vec)/norm(c_1k_hat_new_vec);
            if dk < stopDelta
                disp('Stopping delta reached')
                break
            end
        end
            
        c_1k_hat_old_vec = c_1k_hat_vec;
            
        hist.timer(iter) = toc(time_init);
        if hist.timer(iter) > maxTime
            break
        end
        
        %Onsager correction in wavelet domain the VAMP way
        alpha_1k = mean(df(:))/2;   % divergence          
        r_2k_vec =  c_1k_hat_old_vec - alpha_1k*r_1k_vec;
        Cdiv = 1./(1 - alpha_1k);
        r_2k_vec = r_2k_vec.*Cdiv;
        sigma_2k=mean(abs(r_2k_vec(:)-c0_vec(:)).^2); 
        
        % save hist
        hist.true_err_r_1k(iter) = sqrt(mean(abs(r_1k_vec(:)-c0_vec(:)).^2));
        hist.belief_std_r_1k(iter) = sqrt(sigma_1k);
        hist.true_err_c_1k_hat(iter) = sqrt(mean(abs(c_1k_hat_vec(:)-c0_vec(:)).^2));
        hist.belief_std_c_1k_hat(iter) = sqrt(sigma_1k*alpha_1k);
        hist.true_err_r_2k(iter)=sqrt(mean(abs(r_2k_vec(:)-c0_vec(:)).^2));
        hist.Cdiv(iter) = Cdiv;
            
       
        
        
        
        if verbose % band-wise progress of thresholded image 
             x_est = multiscaleRecon(vec2cell(c_1k_hat_old_vec, l,scales));
            fprintf('mean(||x_hat-x0||^2): true RMSE = %f\n', sqrt(mean(abs(x_est(:)-x0(:)).^2)));
             
           
            fprintf( ...
            'true RMSE = %f\tbelief std = %f', ...
            sqrt(mean(abs(c_1k_hat_old_vec(:)-c0_vec(:)).^2)), sqrt(err));        
        end
      %%  
        %CG-VAMP
        [fxnA,fxnAt] = operatorsGenerate(mask,scales,prob_map,l); % generate A and A'
        z = y(:)-fxnA(r_2k_vec); % z = z-A*r2k
        fxnAAt = @(z) fxnA(fxnAt(z)); % generate AA'       
        q_t = r_2k_vec - c0_vec;  % r2k =c0+ q_t
        cor2(iter) = (h_t')*(q_t)/norm((h_t))/norm((q_t)); % correlation
        fxnW=@(z)sigma_w * z + sigma_2k * fxnAAt(z); % W = sigma_w + sigma_2k_avg* AA'
        t_cg =tic;
        CG_estimate = CG2(z, fxnW,  nx*ny, 20 , 0); % Using CG to estimate  W^{-1}z
        hist.cg_timer(iter) =toc(t_cg);
        At_CG_estimate = fxnAt(CG_estimate); % A'W^{-1}z
        c_2k_hat_vec =  r_2k_vec + sigma_2k * At_CG_estimate; 
        gamma_t = - q_t' * At_CG_estimate/sigma_2k/nx/ny;
        r_1k_vec = r_2k_vec + gamma_t^-1 * At_CG_estimate;
        
        
%%
        r_1k = vec2cell(r_1k_vec,l,scales);
        r = multiscaleRecon(r_1k); % compute r to update the true RMSE in K-space message
        h_t = r_1k_vec -c0_vec;
        cor1(iter) = (h_t)'*(q_t)/norm((h_t))/norm((q_t));            
        sigma_1k=mean(abs(r_1k_vec(:)-c0_vec(:)).^2);
        
        % saving the history of mean(||c_2k_hat-c0||^{2}) in each
        % band
        hist.true_err_c_2k_hat(iter) = sqrt(mean(abs(c_2k_hat_vec(:)-c0_vec(:)).^2));
           
    end
    
 
    % retrospectively calculate x estimate and mse 
    cor1 = abs(cor1);
    cor2 = abs(cor2);
    for iter = 1:numel(hist.timer)
        xk_tilde = multiscaleRecon(pyramidInv(hist.c_1k_hat(:,:,iter), opts.scales));
        z2 = mask.*(y - fftnc(xk_tilde));
        xk = ifftnc(fftnc(xk_tilde) + z2);
            hist.x(:,:,iter) = xk;
        hist.x_mse(iter) = immse(xk,x0);
        hist.x_nmse(iter)= 10*log10(norm(xk-x0).^2/(norm(x0).^2));%
    end
    
    hist.x_mse = hist.x_mse(hist.x_mse>0);
    
    x_hat = xk;
end




