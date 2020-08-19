function [x_hat, hist,cor1,cor2] = VDAMP_O(y, mask, prob_map, sigma_w, x0, opts)
%     The Variable Density Approximate Message Passing (VDAMP) algorithm for
%     reconstruction of a natural image from Fourier coefficients, where
%     a sparse model on Haar wavelet coefficients is assumed
%     IN:
%         y: (nx*ny) measurement data i.e. noisy undersampled Fourier coefficients...
%             unsampled coefficients have entry zero
%         mask: (nx*ny) matrix of zeros and ones that signifies of sampling location
%         prob_map: (nx*ny) matrix of sampling probabilities
%         sigma_w: variance of measurement noise 
%         x0: (nx*ny) ground truth/fully sampled image for tracking reconstruction progress
%         opts: options object with attributes
%             maxIter: maximum number of allowed iterations; default 50
%             verbose: if 1, print progress metrics; default 0
%             scales: number of wavelet decomposition scales; default 4    
%             saveHist: if 1, save detailed history of reconstruction; default 0
%             SURE: if 1, use SURE to automatically tune thresholds; 
%             default 1
%             lambda: if not using SURE, specify sparse weighting
%             denoiserDiv: if 1, use SURE to select denoising divisor Cdiv,
%             otherwise use (1-alpha) this is what we use in dissertation 
%             stopDelta: criterion that stops algorithm at iteration k when ||c1_k_hat - c1_{k-1}_hat||_2/||c1_k_hat||_2 < stopDelta. Default is 0   
%
%      OUT:
%         x_hat: VDAMP's x0 estimate
%         hist:  history object formed with attributes:
                
%               timer: (maxIter) time at every iteration
%               Cthr: (nx,ny,maxIter) post-thresholding estimate in wavelet domain
%               x_mse: (maxIter) ground truth mean-squared error of x
%               x_nmse: ground truth nmse of x
%               **the following are saved only if saveHist == 1**:
%                   r1: (nx*ny*maxIter) image representation of vector subject to thresholding
%                   true_err_r_1k: (maxIter*4*scales) band-wise ground truth RMSE of r1k, order H, V, D, A...
%                            Those scales without A have zeros
%                   belief_std_r_1k: (maxIter*4*scales) estimate of
%                   band-wise RMSE of r1k
%                   RMSEyest: (maxIter) root-mean of tau i.e. estimate of RMSE of r
%                   true_err_c_1k_hat: (maxIter*4*scales) band-wise ground truth RMSE of w (thresholded r1)
%                   belief_std_c_1k_hat: (maxIter*4*scales) estimate of band-wise RMSE of w         
%                   x: (nx*ny*maxIter) image estimate
%                   x_idx: index at minimum mean-squared error estimate i.e. index of best guess of x
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2019  Charles Millard
% Copyright (C) 2019  Siemens Healthineers     


% Notations were adapted to fit in the dissertation @ Chengpin Luo 
                
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
    
    if (isfield(opts,'saveHist') && ~isempty(opts.saveHist))
        saveHist = opts.saveHist;
    else
        saveHist = 0;
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
    
    if sureFlag  == 0 
         if (isfield(opts,'lambda') && ~isempty(opts.lambda))
             sparse_weight = opts.lambda;
         else
             error('SURE is off and no lambda has been selected: please set opts.lambda')
         end
    end
    
        
    sigma_1k = cell(scales, 1); 
    lambda = cell(scales, 1);
    err = cell(scales, 1);
    r_2k= cell(scales, 1);
    alpha_1k = cell(scales, 1);
    for s =1:scales
        sigma_1k{s} = struct('H', 0, 'V', 0, 'D', 0);
        err{s} = struct('H', 0, 'V', 0, 'D', 0);
        r_2k{s} = struct('H', [], 'V', [], 'D', []);
        alpha_1k{s} = struct('H', 0, 'V', 0, 'D', 0);       
        if sureFlag ==0
            lambda{s} = struct('H', sparse_weight, 'V', sparse_weight, 'D', sparse_weight);
        end       
    end  
    sigma_1k{scales}.A = 0;
    err{scales}.A = 0;
    r_2k{scales}.A = [];
    alpha_1k{s}.A = 0;   
    if sureFlag == 0
        lambda{scales}.A = sparse_weight/100;
    end
     
    % wavelet coefficients of original image
    c0 = multiscaleDecomp(x0, scales);
    
    
    % component-wise inverse of the probability map
    % notice that prob_map is a nx by ny matrix
    % not the diagonal nxny by nxny matrix in the paper
    inv_p = prob_map.^(-1);
    % unbaised initialisation
    m2 = inv_p.*mask.*y;
    % here we denote r as the time domain representation of r1k
    r = ifftnc(m2);
    r_1k = multiscaleDecomp(r, scales); % the input r1k to the denoiser g1
    % use oracle information for initilization of r_1k for each band
    for s = 1:scales
        subbands = fieldnames(r_1k{s});               
        for i = 1:numel(subbands)
            band_name = subbands{i};
            sigma_1k{s}.(band_name)=mean(abs(r_1k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2);
        end
    end
    % prepare history
    if saveHist
        hist.r1 = zeros(nx,ny,maxIter);
        hist.x = zeros(nx,ny,maxIter);
        hist.true_err_r_1k= zeros(maxIter, 4, scales);
        hist.belief_std_r_1k = zeros(maxIter, 4, scales);
        hist.true_err_c_1k_hat= zeros(maxIter, 4, scales);
        hist.belief_std_c_1k_hat = zeros(maxIter, 4, scales);
        hist.Cdiv = zeros(maxIter, 4, scales);
        hist.RMSEyest = zeros(1, maxIter);
        hist.x_mse = zeros(1, maxIter);
        hist.x_nmse = zeros(1, maxIter);
    end
    
    hist.C_thr = zeros(nx,ny,maxIter);
    
    
    time_init = tic;
    % Notice r_1k is a cell, which could be convert into a vector
    % representation as r_1k_vec
    r_1k_vec = cell2vec(r_1k);
    
    % h_t = r1k- c0
    h_t = r_1k_vec -cell2vec(c0);
    cor1 = zeros(maxIter,1); % input correlation in the thesis
    cor2 = zeros(maxIter,1); % output correlation in the thesis
    for iter = 1:maxIter
        
        % remember r is the time-domain representation of r1k
        if saveHist
            hist.r1(:,:,iter) = r;
        end
        
        % printing messages
        if verbose % band-wise mse of unthresholded image
            disp(['                 ITERATION ', num2str(iter)])
            fprintf( ...
                ' idwt(r1k)-x0: true RMSE = %f\n', ...
                sqrt(mean(abs(r(:)-x0(:)).^2)));
            for s=1:scales
                subbands = fieldnames(r_1k{s});
                fprintf(['\tAt scale ' int2str(s) '\n']);
                for i = 1:numel(subbands)
                    band_name = subbands{i};
                    fprintf( ...
                        '\t\tIn band %s: true RMSE = %f\tmessage std = %f\n', ...
                        band_name, ...
                        sqrt(mean(abs(r_1k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2)), ...
                        sqrt(sigma_1k{s}.(band_name)));
                end
            end
        end
        
        % thresholding, the denoising block 
        if sureFlag
            [c_1k_hat, err, df] = multiscaleSUREsoft(r_1k, sigma_1k);
        else      
            [c_1k_hat, err, df] = multiscaleComplexSoft(r_1k, sigma_1k, lambda);
        end
        
        % pyramid representation of c_1k_hat, notice c_1k_hat is a cell
        hist.c_1k_hat(:,:, iter) = pyramid(c_1k_hat);
       
        % check stopping criterion
        if stopDelta >0 && iter> 1
            c_1k_hat_old_py = pyramid(c_1k_hat_old);
            c_1k_hat_new_py = pyramid(c_1k_hat);
            dk = norm(c_1k_hat_old_py(:) - c_1k_hat_new_py(:))/norm(c_1k_hat_new_py(:));
            if dk < stopDelta
                disp('Stopping delta reached')
                break
            end
        end
            
        c_1k_hat_old = c_1k_hat;
            
        hist.timer(iter) = toc(time_init);
        if hist.timer(iter) > maxTime
            break
        end
        
        %Onsager correction in wavelet domain the VAMP way
        for s=1:scales
            subbands = fieldnames(c_1k_hat_old{s});           
            for i = 1:numel(subbands)
                band_name = subbands{i};
                
                alpha_1k{s}.(band_name) = mean(df{s}.(band_name)(:))/2;             
                r_2k{s}.(band_name) =  (c_1k_hat_old{s}.(band_name) - alpha_1k{s}.(band_name)*r_1k{s}.(band_name));
                
                if denoiserDiv == 1    
                    % VDAMP-S
                    Cdiv = r_1k{s}.(band_name)(:)'*(r_2k{s}.(band_name)(:))./norm(r_2k{s}.(band_name)(:))^2;
                else
                    % VDAMP-alpha--- we only use this in our dissertation
                    Cdiv = 1./(1 - alpha_1k{s}.(band_name));
                end
                
                r_2k{s}.(band_name) = r_2k{s}.(band_name).*Cdiv;
                
                % saving history
                % in this dissertaion, you may only care about the
                % hist.true* because we only use oracle information
                if saveHist         
                    hist.true_err_r_1k(iter,i,s) = sqrt(mean(abs(r_1k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2));
                    hist.belief_std_r_1k(iter,i,s) = sqrt(sigma_1k{s}.(band_name));
                    hist.true_err_c_1k_hat(iter,i,s) = sqrt(mean(abs(c_1k_hat{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2));
                    hist.belief_std_c_1k_hat(iter,i,s) = sqrt(sigma_1k{s}.(band_name)*alpha_1k{s}.(band_name));
                    hist.Cdiv(iter,i,s) = Cdiv;
                end        
            end
        end
 
        if verbose % band-wise progress of thresholded image 
             x = multiscaleRecon(c_1k_hat_old);
            fprintf('mean(||x-x0||^2): true RMSE = %f\n', sqrt(mean(abs(x(:)-x0(:)).^2)));
            for s=1:scales
                subbands = fieldnames(r_1k{s});
                fprintf(['\tAt scale ' int2str(s) '\n']);
                for i = 1:numel(subbands)
                    band_name = subbands{i};
                    fprintf( ...
                        '\t\tIn band %s: true RMSE = %f\tbelief std = %f\tSURE = %f\n', ...
                        band_name, ...
                        sqrt(mean(abs(c_1k_hat_old{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2)), sqrt(err{s}.(band_name)),...
                        sqrt(mean(abs(c_1k_hat_old{s}.(band_name)(:)-r_1k{s}.(band_name)(:)).^2)-sigma_1k{s}.(band_name) + 2*err{s}.(band_name)));
                end
            end
        end
        [r_2k_vec,l]= cell2vec(r_2k); % r_2k into vector form
        [c0_vec]=cell2vec(c0); %c0 into vector form
        q_t = r_2k_vec - c0_vec; % q_t =r2k-c0
        cor2(iter) = (h_t')*(q_t)/norm((h_t))/norm((q_t)); % output correlation
        % ****Reweighted gradient step****
        r_tilde = multiscaleRecon(r_2k); % r_tilde is the time domain of r_2k
        z = mask.*(y - fftnc(r_tilde));
        
        
        r = r_tilde + ifftnc(inv_p .* z);  % r is the time domain of r_1k
        r_1k = multiscaleDecomp(r, scales);
        r_1k_vec =cell2vec(r_1k);
        h_t = r_1k_vec -c0_vec; 
        cor1(iter) = (h_t')*(q_t)/norm((h_t))/norm((q_t)); %compute input correlation
    % estimate noise variance sigma_1k for next iteration
    for s = 1:scales
            subbands = fieldnames(r_2k{s});
            for i = 1:numel(subbands)
                band_name = subbands{i};
                sigma_1k{s}.(band_name)=mean(abs(r_1k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2);
            end
    end 
    end  % VDAMP MAIN ITERATION END
    % retrospectively calculate x estimate and mse 
    cor1 =abs(cor1);
    cor2 =abs(cor2);
    for iter = 1:numel(hist.timer)
        xk_tilde = multiscaleRecon(pyramidInv(hist.c_1k_hat(:,:,iter), opts.scales));
        z2 = mask.*(y - fftnc(xk_tilde));
        xk = ifftnc(fftnc(xk_tilde) + z2); % last line in the vdamp
        if saveHist
            hist.x(:,:,iter) = xk;
        end
        hist.x_mse(iter) = immse(xk,x0);%MSE
        hist.x_nmse(iter)= 10*log10(norm(xk-x0).^2/(norm(x0).^2));%NMSE
    end
    
    hist.x_mse = hist.x_mse(hist.x_mse>0);
    
    x_hat = xk; % FINAL OUTPUT OF RECONSTRUCTED IMAGE
end

