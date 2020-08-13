function [img] = idwt2v(c,l,scales)
%IDWT2V inverse 2d dwt for vector-form image
% vector-in-vector-out
C_cell= vec2cell(c,l,scales);
img = multiscaleRecon(C_cell);
img = img(:);
end

