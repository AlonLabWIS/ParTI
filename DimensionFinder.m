function [dim] = DimensionFinder(ESVdat)
%DimensionFinder looks for the correct dimensionality.
%   ESVdat is a vector of doubles

slope = (ESVdat(end) -  ESVdat(1))/(length(ESVdat)-1);
intersect =  ESVdat(1)-slope * 1;
di = zeros(length(ESVdat),1);

for i = 1:length(ESVdat)
    s = -1/slope;
    inter2 =  ESVdat(i) - s * i;
    x  = (inter2 - intersect)/  (slope - s) ; 
    y  =  s * x + inter2;
    di(i) = norm([x - i,y - ESVdat(i)]);
end
 
[~, dim] = max(di);
end

