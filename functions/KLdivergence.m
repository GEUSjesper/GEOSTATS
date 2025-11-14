function D = KLdivergence(A,B)

A_array = reshape(A,[],1);
B_array = reshape(B,[],1);

if numel(A) ~= numel(B)
    error('Probability matrices must have same size')
end

idx = A_array ~= 0 & B_array ~= 0;

D = sum(A_array(idx).*log(A_array(idx)./B_array(idx)));

    