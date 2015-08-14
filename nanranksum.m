function [p] = nanranksum(x, y)

keepX = find(~isnan(x));
keepY = find(~isnan(y));

p = ranksum(x(keepX),y(keepY));