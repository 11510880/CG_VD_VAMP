function [C_cell] = vec2cell(C_vec, l,scales)
%VEC2CELL convert vector-form coefficients to cell-form
% vector structure [A;H;V;D;H;V;D...]
% l: list that indicates the number of coefficients in each sub-band 
C1 =C_vec;
C = cell(scales, 1);
        l_sum = zeros(length(l),1);
        l_sum(1) = l(1);
        for i=2:length(l)
            l_sum(i)=l_sum(i-1)+ l(i);
        end
        for s = 1:scales % notice scales is at the start(top) of the vector
            index =3*(scales-s+1)+1;
            index_D = l_sum(index)-l(index)+1:l_sum(index); %index for D
            index_V = l_sum(index-1)-l(index-1)+1:l_sum(index-1); % index for V
            index_H = l_sum(index-2)-l(index-2)+1:l_sum(index-2); % index for H
            nD = sqrt(length(index_D));
            nV = sqrt(length(index_V));
            nH = sqrt(length(index_H));
            C{s} = struct( ...
                'D', reshape(C1(index_D),[nD,nD]), ...
                'V', reshape(C1(index_V),[nV,nV]), ...
                'H', reshape(C1(index_H),[nH,nH]));
        end
        C{scales}.A = reshape(C1(1:l(1)),[nD,nD]);
C_cell = C;
end

