function out = integrate_powerlaw_vs_ideal(tbv,tsv,t0v,beta)

tbv = tbv(:);
tsv = tsv(:);
t0v = t0v(:);

% The purpose: Integrate the difference of the sinc filter minus the ideal
% filter for different taus. How big are the differences? How do they
% change with beta?

% function out = integrate_powerlaw(tbv,tsv,beta)
% tbv is a vector of tau_b values
% tsv is a vector of tau_s values
% t0v is a vector of record lengths
% beta is the negative spectral slope

if max(size(t0v))==1,t0v=t0v*ones(length(tbv));end

fun1 = @(x,taub,b) x.^-b.*(1-sinc(taub.*x)).^2;
fun2 = @(x,taub,b) x.^-b.*(sinc(taub.*x)).^2;
%fun3 = @(x,b) x.^-b;

out = [];

for ii = 1:length(tbv)
    taus = tsv(ii);
    taub = tbv(ii);
    tau0 = t0v(ii);
    out(ii) = (integral(@(x)fun1(x,taub,beta),1/tau0,1./(2*taus))+...
              integral(@(x)fun2(x,taub,beta),1./(2*taus),Inf))./...
              integral(@(x)fun2(x,taub,beta),1/tau0,Inf);
end
