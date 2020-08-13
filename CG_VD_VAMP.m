function [x_hat, cor1, cor2,hist] = CG_VD_VAMP(y, mask, prob_map, sigma_w, x0, opts)
%  THE CG-VD-VAMP ALGORITHM PROPOSED IN THE THESIS @ 2020 CHENGPIN LUO
%  output: x_hat: the reconstructed image
%  cor1: input correlation
%  cor2: output correlation
%  hist: struct to save information about updates

%  APDATED FROM VDAMP @ Copyright (C) 2019  Charles Millard
%                     @ Copyright (C) 2019  Siemens Healthineers

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
%             otherwise use (1-alpha) 
%             stopDelta: criterion that stops algorithm at iteration k when ||C_thr_k - C_thr_{k-1}||_2/||C_thr_k||_2 < stopDelta. Default is 0   
%
%      OUT:
%         x_hat: VDAMP's x0 estimate
%         hist:  history object formed with attributes:
%               timer: (maxIter) time at every iteration
%               cg_timer =
%               c1k_hat_old: (nx,ny,maxIter) post-thresholding estimate in wavelet domain
%               x_mse: (maxIter) ground truth mean-squared error of x              
%               **the following are saved only if saveHist == 1**:
%                   r1: (nx*ny*maxIter) image representation of vector subject to thresholding
%                   true_err_r1k: (maxIter*4*scales) band-wise ground truth RMSE of r, order H, V, D, A...
%                            Those scales without A have zeros
%                   belief_std_r2k: (maxIter*4*scales) estimate of band-wise RMSE of r
%                   RMSEyest: (maxIter) root-mean of tau i.e. estimate of RMSE of r
%                   true_err_c1k_hat: (maxIter*4*scales) band-wise ground truth RMSE of w (thresholded r1)
%                   belief_std_c1k_hat: (maxIter*4*scales) estimate of band-wise RMSE of w         
%                   x: (nx*ny*maxIter) image estimate
%                   x_idx: index at minimum mean-squared error estimate i.e. index of best guess of x
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2019  Charles Millard
% Copyright (C) 2019  Siemens Healthineers

% the vector version-- all the transformation would be done in vector form
                
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
    
    if sureFlag  == 0 
         if (isfield(opts,'lambda') && ~isempty(opts.lambda))
             sparse_weight = opts.lambda;
         else
             error('SURE is off and no lambda has been selected: please set opts.lambda')
         end
    end
    
    hist.r1 = zeros(nx,ny,maxIter); % r1 is the idwt of r_1k
    hist.x = zeros(nx,ny,maxIter);
    hist.true_err_r_1k= zeros(maxIter, 4, scales); % mean(||r1k-c0||.^2)
    hist.belief_std_r_1k = zeros(maxIter, 4, scales);
    hist.true_err_c_1k_hat= zeros(maxIter, 4, scales);
    hist.belief_std_c_1k_hat = zeros(maxIter, 4, scales);% mean(||r1k_hat-c0||.^2)
    hist.true_err_r_2k= zeros(maxIter, 4, scales);
    hist.belief_std_r_2k = zeros(maxIter, 4, scales);% mean(||r2k-c0||.^2)
    hist.true_err_c_2k_hat= zeros(maxIter, 4, scales);% mean(||r2k_hat-c0||.^2)
    hist.belief_std_c_2k_hat = zeros(maxIter, 4, scales);
    hist.Cdiv = zeros(maxIter, 4, scales);
    hist.RMSEyest = zeros(1, maxIter);
    hist.x_mse = zeros(1, maxIter);   
    hist.x_nmse=zeros(1, maxIter);
    hist.cg_timer=zeros(1,maxIter);
    sigma_1k = cell(scales, 1);
    sigma_2k=cell(scales,1);
    lambda = cell(scales, 1);
    err = cell(scales, 1);
    r_2k= cell(scales, 1);
    hist.c_1k_hat = zeros(nx,ny,maxIter);
    alpha_1k = cell(scales, 1);
    for s =1:scales
        sigma_1k{s} = struct('H', 0, 'V', 0, 'D', 0);
        sigma_2k{s} = struct('H', 0, 'V', 0, 'D', 0);
        err{s} = struct('H', 0, 'V', 0, 'D', 0);
        r_2k{s} = struct('H', [], 'V', [], 'D', []);
        alpha_1k{s} = struct('H', 0, 'V', 0, 'D', 0);       
        if sureFlag ==0
            lambda{s} = struct('H', sparse_weight, 'V', sparse_weight, 'D', sparse_weight);
        end       
    end  
    sigma_1k{scales}.A = 0;
    sigma_2k{scales}.A=0;
    err{scales}.A = 0;
    r_2k{scales}.A = [];
    alpha_1k{s}.A = 0;   
    if sureFlag == 0
        lambda{scales}.A = sparse_weight/100;
    end
  
    c0 = multiscaleDecomp(x0, scales);
    
    
    inv_p = prob_map.^(-1);   %P^(-1)
    % unbaised initialisation
    m2 = inv_p.*mask.*y;
    r = ifftnc(m2);
    r_1k = multiscaleDecomp(r, scales);

    for s=1:scales
        subbands = fieldnames(r_1k{s});
        for i = 1:numel(subbands)
            band_name = subbands{i};
            % using the oracle information to initialize the sigma_1k for
            % each band
            sigma_1k{s}.(band_name)=mean(abs(r_1k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2);
        end
    end
    r_1k_vec = cell2vec(r_1k);
    h_t = r_1k_vec -cell2vec(c0); %r1k = c0+h_t
    cor1 = zeros(maxIter,1);
    cor2 = zeros(maxIter,1);
    time_init = tic;
    for iter = 1:maxIter
        hist.r1(:,:,iter) = multiscaleRecon(r_1k);
        if verbose % band-wise mse of unthresholded image
            disp(['                 ITERATION ', num2str(iter)])
            fprintf( ...
                '(iwdt(r1k)-x0): true RMSE = %f\n', ...
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
                    % VDAMP-alpha  we only use this in dissertation
                    Cdiv = 1./(1 - alpha_1k{s}.(band_name));
                end
                
                r_2k{s}.(band_name) = r_2k{s}.(band_name).*Cdiv;
                % using oracle information for estimating sigma_2k
                sigma_2k{s}.(band_name)=mean(abs(r_2k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2); 
                %sigma_2k{s}.(band_name) = sigma_1k{s}.(band_name)*alpha_1k{s}.(band_name)/(1-alpha_1k{s}.(band_name));
                % save hist
                hist.true_err_r_1k(iter,i,s) = sqrt(mean(abs(r_1k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2));
                hist.belief_std_r_1k(iter,i,s) = sqrt(sigma_1k{s}.(band_name));
                hist.true_err_c_1k_hat(iter,i,s) = sqrt(mean(abs(c_1k_hat{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2));
                hist.belief_std_c_1k_hat(iter,i,s) = sqrt(sigma_1k{s}.(band_name)*alpha_1k{s}.(band_name));
                hist.true_err_r_2k(iter,i,s)=sqrt(mean(abs(r_2k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2));
                hist.Cdiv(iter,i,s) = Cdiv;
            end
        end
        
        
        
        if verbose % band-wise progress of thresholded image 
             x_est = multiscaleRecon(c_1k_hat_old);
            fprintf('time-domain belief: true RMSE = %f\n', sqrt(mean(abs(x_est(:)-x0(:)).^2)));
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
        %% THIS IS THE CORE PART OF OUR CG-VD-VAMP ALGORITHM
        [r_2k_vec,l]= cell2vec(r_2k); % tranform the cell-from r2k into vector-form
        [c0_vec]=cell2vec(c0); % transform the oracle information c0(WAVELET OF x0) into vector form 
        [fxnA,fxnAt] = operatorsGenerate(mask,scales,prob_map,l); % generate A and A'
        z = y(:)-fxnA(r_2k_vec); % z = z-A*r2k
        % now we are going to build the matrix Sk for N(c_2;r2k,Sigma_2k)
        % Sigma_2k's diagonal is filled with sigma_2k_i bandwise
        Sigma_2k = r_2k;
        for s=1:scales
                subbands = fieldnames(r_2k{s});
                for i = 1:numel(subbands)
                    band_name = subbands{i};
                    %nums =length(r_2k{s}.band_name);
                    Sigma_2k{s}.(band_name)(:) =sigma_2k{s}.(band_name);
                end
        end
        Sigma_2k_vec=cell2vec(Sigma_2k); % Convert Sigma_2k into vector form, i.e. taking its diagonal
        fxnASAt =@(z)fxnA(Sigma_2k_vec.*fxnAt(z)); % A' Sigma_2k  A'
        q_t = r_2k_vec - c0_vec;  % r2k =c0+ q_t
        cor2(iter) = (h_t')*(q_t)/norm((h_t))/norm((q_t)); % output correlation
        fxnW=@(z)sigma_w * z + fxnASAt(z); % operator for (sigma_w I_{M} +  A' Sigma_2k_vec A' )
        t_cg =tic; % start to compute the time each 20 inner cg iteration takes
        CG_estimate = CG2(z, fxnW,  nx*ny, 20 , 0); % Using CG to estimate  W^{-1}z
        hist.cg_timer(iter) =toc(t_cg);% end toc result saved in hist.cg_timer
        At_CG_estimate = fxnAt(CG_estimate); % A'u_tilde
        % c_2k_hat =r_2k + Sigma_2k A'u_tilde line 14 in algorithm 4
        c_2k_hat_vec = r_2k_vec + Sigma_2k_vec.*At_CG_estimate;
        c_2k_hat = vec2cell(c_2k_hat_vec,l,scales); % cell form- c_2k_hat
        % update the parameters using cell-form
        % estimate r_1k for next iteration using oracle information
        for s=1:scales
                subbands = fieldnames(r_2k{s});
                for i = 1:numel(subbands)
                    band_name = subbands{i};
                    % line 16 
                    s_2k = mean(abs(c_2k_hat{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2);
                    % line 17
                    r_1k{s}.(band_name)(:) = sigma_2k{s}.(band_name)*c_2k_hat{s}.(band_name)(:)-s_2k*r_2k{s}.(band_name)(:);
                    r_1k{s}.(band_name)(:) =r_1k{s}.(band_name)(:)/(sigma_2k{s}.(band_name)-s_2k);                 
                end
        end
       
%%      
        r_1k_vec =cell2vec(r_1k);
        %r_1k = vec2cell(r_1k_vec,l,scales);
        r = multiscaleRecon(r_1k); % r is time domain of r_1k
        h_t = r_1k_vec -c0_vec;
        cor1(iter) = (h_t')*(q_t)/norm((h_t))/norm((q_t)); % input correlation
    for s = 1:scales
            subbands = fieldnames(r_2k{s});
            for i = 1:numel(subbands)
                band_name = subbands{i};
                % estimate sigma_1k for next iteration
                sigma_1k{s}.(band_name)=mean(abs(r_1k{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2);
                % saving the history of mean(||c_2k_hat-c0||^{2}) in each
                % band
                hist.true_err_c_2k_hat(iter,i,s) = sqrt(mean(abs(c_2k_hat{s}.(band_name)(:)-c0{s}.(band_name)(:)).^2));
            end
    end 
    end
    
 
    % retrospectively calculate x estimate and mse 
    cor1 =abs(cor1);
    cor2 =abs(cor2);
    for iter = 1:numel(hist.timer)
        xk_tilde = multiscaleRecon(pyramidInv(hist.c_1k_hat(:,:,iter), opts.scales));
        z2 = mask.*(y - fftnc(xk_tilde));
        xk = ifftnc(fftnc(xk_tilde) + z2);
            hist.x(:,:,iter) = xk;
        hist.x_mse(iter) = immse(xk,x0);%mse
        hist.x_nmse(iter)= 10*log10(norm(xk-x0).^2/(norm(x0).^2));%nmse
    end
    
    hist.x_mse = hist.x_mse(hist.x_mse>0);
    
    x_hat = xk; % the final output
end


