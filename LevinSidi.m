function [A,B] = LevinSidi (k, S, w, x)
% Returns the array storing counterdiagonals 
% k is the transformation order
% S is the sequence of sum

n = length(S);
A = zeros(k+1);
B = zeros(k+1);
reci_x = zeros(length(x)+2);
for i=1:k+1
    A(i,1) = S(i)/w(i); % Second column has x in it
    B(i,1) = 1/w(i); 
    reci_x(i+1) = 1/x(i);
end
% for i = k : -1 : 1
%     for j = 2 :  k - i + 1
%         d = 1/k - 1/(k-j);
%         A(i,j) = (A(i+1,j-1) - A(i,j-1))*d;
%         B(i,j) = (B(i+1,j-1) - B(i,j-1))*d;
%     end
% 
% end
j = 1;
for i = k+1 : -1 : 1 
            
         d = 1/k - 1/(k-j);
        A(j,i) = (A(j+1,i-1) - A(j,i-1))*d;
        B(j,i) = (B(j+1,i-1) - B(j,i-1))*d;
j = j +1;
end

end