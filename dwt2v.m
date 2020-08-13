function [c,l] = dwt2v(img,scales,m,n)
%vector-form 2d-dwt2
% img: input image , a vector
% scales: dwt2 scales
% [m,n] size of img in matrix form
% output; c, the 2d-dwt wavelet coefficients as a vector
% l: list of number of coefficients in each band
img = reshape(img,[m,n]);
bands = multiscaleDecomp(img, scales);
%c = pyramid(bands); we should not use the pyramid structure
% instead, we should use the stack structure
scales = numel(bands);
c = bands{scales}.A(:);
l=zeros(3*scales+1,1);
l(1)=length(c);
index =1;
for s = scales:-1:1
    oldindex =index;
    index =oldindex+3;
    l(index-2)=length(bands{s}.H(:));
    l(index-1)=length(bands{s}.V(:));
    l(index)=length(bands{s}.D(:));
    c = [c;bands{s}.H(:); bands{s}.V(:);bands{s}.D(:)];
    
end

