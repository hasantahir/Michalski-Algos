function [W,e]=EpsilonAlgorithmLowStorage(x,k)
% EPSILONALGORITHMLOWSTORAGE epsilon algorithm with low storage
% [w,e]=EpsilonAlgorithmLowStorage(x,k) is an implementation of the
% epsilon algorithm of P. Wynn using only little storage. It stores
% only the diagonal of the epsilon table in W, and the last row of
% the triangular array in e in reverse order, computing 2k steps.
e(1)=x(1);
for i=2:2*k+1
    v=0; 
    e(i)=x(i);
    for j=i:-1:2
        d=e(j)-e(j-1);
        w=v+1/d;
        v=e(j-1);
        e(j-1)=w;
    end;
    W(i-1)=w;
end
