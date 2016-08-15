function [S,w] = LevinInit(S,y)

% t-transformation
w = zeros(size(S));
u = zeros(size(S));

for i = 1 : length (S)
    if i == 1
        u(i) = S(i);
    else 
        u(i) = S(i) - S(i-1);
    end
    w(i) = u(i); % t-transformation
%     w(i) = u(i).*y(i); % u-transformation
end
end