function [C_vec,l] = cell2vec(C_cell)
%CELL2VEC cell-form coefficients to vector form coefficients
% for the vector, the start the coarsest scale [A; H; V; D; H; V; D;....]
bands = C_cell;
scales = numel(bands);
img = bands{scales}.A(:);
l=zeros(3*scales+1,1);
l(1)=length(img(:));
index =1;
for s = scales:-1:1
    oldindex =index;
    index =oldindex+3;
    l(index-2)=length(bands{s}.H(:));
    l(index-1)=length(bands{s}.V(:));
    l(index)=length(bands{s}.D(:));
    img = [img;bands{s}.H(:); bands{s}.V(:);bands{s}.D(:)];
end
C_vec = img;
end

