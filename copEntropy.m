function [S] = copEntropy(vec)

vec = vec(~isnan(vec));
counts = zeros(5,1);
for i = -2:2
    counts(i+3) = sum(vec == i);
end
p = counts / sum(counts);
avgLogP = p.*log2(p);
avgLogP(isnan(avgLogP)) = 0;
S = -sum(avgLogP);