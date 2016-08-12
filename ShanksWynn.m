function E = ShanksWynn (k, S)
% Returns the array storing counterdiagonals 
% k is the transformation order
% S is the sequence of sum

E(1)=S(1);

for i = 2 : 2*k + 1
    b = 0; 
    E(i) = S(i);

    for j = i: -1 :2
        a = b;
        b = E(j-1);
        d = E(j) - b;
        E(j-1) = a + 1/d;
    end;
end


