function out = checkIdenticalRows_old(v,M)

[n_rows,n_cols] = size(M);
identical = 0*n_rows;

for i = 1:n_cols
    identical = identical + (M(:,i) == v(i));
    if max(identical) < i
        out = 0;
        return
    end
end

out = 1;