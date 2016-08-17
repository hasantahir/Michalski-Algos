function s = PartSum_mix_PE(fun,a,eh, e2h, n)
% This function performs the summation of the terms generated by the
% product of weights and nodes
f = fun;
ekh = eh;
s = Term_mix_PE(f, a,ekh);
for k = 2 : n
    ekh = ekh * e2h;
    s = s + Term_mix_PE(f, a,ekh);
end
end
