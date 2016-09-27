%MULLERBEEF (mfile) for excecuting Muller's method.
% Edit the m file and set fa fb and fc equal to a function in terms of a b
% and c, but make the function exactly the same for fa fb and fc save the
% name of the variable. The rest of the function will excecute Muller's
% method iterating 10 times.
clear all
close
a=1;
b=2;
c=3;
load logical.mat 
% f = inline('sin(x)+cos(x)', 'x');
%Check to see if ab or c is already the root
if D(a)==0
    root=a
end
if D(b)==0
    root=b
end
if D(c)==0
    root=c
end
for k=1:10000
    temp=a;
    temp2=b;
    fa=D(a);
    fb=D(b);
    fc=D(c);
    q=(a-b)/(b-c);
    %defines q to be used later in the method for convenience
    A=(q*fa)-(q*(1+q)*fb)+(q^2*fc);
    %defines A to be used in the final quadratic formula evaluation
    B=((2*q+1)*(fa))-((1+q)^2*fb)+(q^2*fc);
    %defines B for quad formula
    C=(1+q)*fa;
    %defines c for quad formula
    root=a-(a-b)*((2*C)/(max((B+sqrt(B.^2-(4*A.*C))),(B-sqrt(B.^2-(4*A.*C))))));
    %set a tolerance
    if (B+sqrt(B.^2-(4*A.*C)) == 0)
        break
    end
    if (B-sqrt(B.^2-(4*A.*C)) == 0)
        break
    end
    %rounds d to 0 if it is smaller than our tolerance of 10^-15
    if root<1e-15
        root=0;
    end
    b=temp;
    a=root;
    c=temp2;
end
format long
root