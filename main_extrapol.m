clear all
k=10;
% v=1;
% for j=1:2*k+1
%     y(j)=v/j; 
%     v=-v;
% end

% Example from reference paper
% jmax = 50;
% for j = 1 : jmax
%     y(j) = (-1)^(j-1)/(sqrt(j-1+1));
% end

% Example 2
jmax = 50;
for j = 1 : jmax
    y(j) = (4/5)^(j)/(j);
end

% Make a sequence of sums
S=cumsum(y);

tic
[Er,E2] = EpsilonAlgorithm(S,k);
toc
[ex,ey] = size(E2);
E2(ex,ey)
tic
[W,e] = EpsilonAlgorithmLowStorage(S,k);
toc
W(length(W))
tic
E1 = ShanksWynn(k,S);
toc
E1(1)