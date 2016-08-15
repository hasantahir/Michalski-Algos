function [A,B,ratio] = LevinSidi1 (k, S, w,ratio)
% Returns the array storing counterdiagonals 
% k is the transformation order
% S is the sequence of sum

beta = 1;
term = 1/(beta + k);
A = S./w;
B = 1./w;
if k > 1
    ratio = (beta+ k)*term;
end
    for i = 0 : k-1
        fact = (k - i + beta)*term;
        A(k - i) = A(k - i + 1) - fact*A(k-i);
        B(k - i) = B(k - i + 1) - fact*B(k-i);
        term = term*ratio;
    end
end