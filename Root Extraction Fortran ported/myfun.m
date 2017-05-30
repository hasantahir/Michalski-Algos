function y = myfun(z) 
global pow
f = FZ(z);
% dd = df;
df = dFZ(z);

% save d.mat df dd
y = z.^(pow - 1).*df./f;%/FZ(z); 
end  