function out = checkIdenticalRows(v, M)
    out = any(all(M == v, 2));
end