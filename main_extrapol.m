clear all
k=5;
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
jmax = 11;
for j = 1 : jmax
    y(j) = (4/5)^(j)/(j);
end
%
% Make a sequence of sums
S=cumsum(y);
% [S,w] = LevinInit(S,y);
% tic
% [A,B] = LevinSidi (k, S, w, y);
% toc
% ss = A(1,k)/B(1,k)

% tic
% for i = 1 : jmax   
%     [A,B,ratio] = LevinSidi1(i, S, w,ratio);   
% end
% toc
% ss = A(1)/B(1)

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
% k  = 3;
% a = zeros(k);
% for i = 1 : k
%     for j = 1 : k - i + 1
%         a(j,i) = 1;
%     end
% end
%