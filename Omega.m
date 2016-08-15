function omega = Omega (k ,q, z , a, X)

if k == 1
    omega = 0;
else
    omega = -exp(- q *z).* (X(k-1)/X(k)).^a;
end

end