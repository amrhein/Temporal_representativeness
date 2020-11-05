function out = est_error_powerlaw(tbv,tsv,t0v,beta)

% function out = integrate_powerlaw(tbv,tsv,beta)
% tbv is a vector of tau_b values
% tsv is a vector of tau_s values
% t0v is a vector of record lengths
% beta is the negative spectral slope

if max(size(t0v))==1,t0v=t0v*ones(length(tbv));end

% Get the uncertainty when alpha = 1
fun = @(x,taub,taus,b) x.^-b.*(sinc(taus.*x)-sinc(taub.*x)).^2;
% denominator of the fraction
fund = @(x,taub,b) x.^-b.*(sinc(taub.*x)).^2;
%fun3 = @(x,b) x.^-b;

out = [];

for ii = 1:length(tbv)
    taus = tsv(ii);
    taub = tbv(ii);
    tau0 = t0v(ii);
    out(ii) = integral(@(x)fun(x,taub,taus,beta),1/tau0,Inf)./...
              integral(@(x)fund(x,taub,beta),1/tau0,Inf);

end
