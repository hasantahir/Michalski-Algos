function y = myfun(x) 
global pow
y = x.^(pow - 1)./FZ(x); 
end