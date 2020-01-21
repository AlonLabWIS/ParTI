function [res] = HGfuncReg(a,b,aj,bj)
% Regularized hypergeometric function used to calculate the
% probabily that bin 1 has a higher number of successes than
% bin j

temparg1 = num2cell([1+a;2+a+aj;-b],1);
temparg2 = num2cell([2+a;3+a+aj+bj],1);

%hyp = cellfun(@(x,y)hypergeom(x,y,1),temparg1,temparg2);
if exist('hypergeom') == 0
 error('The Matlab Symbolic Toolbox is needed for this function. Please install the toolbox and try again.')
end
hyp = log(cellfun(@(x,y)hypergeom(x,y,1),temparg1,temparg2));


%  t11 = gamma(2+a+b)./gamma(2+a);
%  t21 = gamma(1+bj)./gamma(1+b);
%  t31 = gamma(2+a+aj)./gamma(3+a+aj+bj);

 t = gammaln(2+a+b)- gammaln(2+a)+ gammaln(1+bj) - gammaln(1+b)+...
     gammaln(2+a+aj)- gammaln(3+a+aj+bj);


%  t1 =  prod(1:1+a+b)./prod(1:1+a);
%  t2 =  prod(1:bj)./prod(1:b);
%  t3 =  prod(1:1+a+aj)./prod(1:2+a+aj+bj);
% for i = 1:length(a);
%  t1(i) =  prod(1:1+a(i)+b(i))./prod(1:1+a(i));
%  t2(i) =  prod(1:bj(i))./prod(1:b(i));
%   t3(i) =  prod(1:1+a(i)+aj(i))./prod(1:2+a(i)+aj(i)+bj(i));
%  t31log(i) = exp((1+a(i)+aj(i))*log(1+a(i)+aj(i)) - (1+a(i)+aj(i)) - ...
%      ((2+a(i)+aj(i)+bj(i))*log(2+a(i)+aj(i)+bj(i)) - (2+a(i)+aj(i)+bj(i))));
% end
%res = 1 - hyp.*beta(1+aj,1+bj).*t1.*t2.*t31log;
 res = 1 - exp( hyp - betaln(1+aj,1+bj) + t );

% res = 1 - cellfun(@(x,y)hypergeom(x,y,1),temparg1,temparg2)... 
% .*gamma(2+a+aj).*gamma(2+a+b).*gamma(1+bj)./...
%     (gamma(2+a).*gamma(3+a+aj+bj).*beta(1+aj,1+bj).*gamma(1+b));
