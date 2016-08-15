function [Er,E]=EpsilonAlgorithm(x,k)
% EPSILONALGORITHM epsilon algorithm of P. Wynn
% [Er,E]=EpsilonAlgorithm(x,k) is an implementation of the epsilon
% algorithm of P. Wynn. The first column of E is zero, the second
% contains the elements of the sequence x whose convergence we want
% to accelerate. 2k colums of the scheme are computed and returned
% in E, while Er is a reduced version of E, containing only the
% relevant columns (every second column).
n=2*k+1;
E=zeros(n+1); % Size increased by 1 to avoid array access errors
for i=1:n
    E(i+1,2)=x(i); % Second column has x in it
end
for i = 3 : n+1
    for j = 3 : i % Lower triangular property
        E(i,j)=E(i-1,j-2)+1/(E(i,j-1)-E(i-1,j-1));
    end
end
Er=E(:,2:2:n+1);