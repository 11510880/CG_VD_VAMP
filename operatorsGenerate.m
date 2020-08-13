function [fxnA,fxnAt] = operatorsGenerate(mask,scales,prob_map,l)
inv_p=prob_map.^(-1);
inv_p = inv_p(:);
fxn_Phit=@(w)idwt2v(w,l,scales);
fxn_FPhit=@(w)fftnc(reshape(fxn_Phit(w),size(mask)));
fxnA =@(w)reshape(fxn_FPhit(w).*mask,[numel(mask),1]);


%mask = mask';inv_p(:).*
fxn_invp=@(y)inv_p(:).*mask(:).*y(:); %suppose y is a vector
fxn_ifft=@(y)ifftnc(reshape(fxn_invp(y),size(mask))); 
%fxnAt=@(y)multiscaleDecomp(fxn_ifft(y),scales);
fxnAt=@(y)dwt2v(fxn_ifft(y),scales,512, 512);
%size(mask,1),size(mask,2)
end

