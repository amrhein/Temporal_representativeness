function out = integrate_powerlaw_vs_ideal_bioturb(tbv,tsv,tdv,t0v,bet)

tbv = tbv(:);
tsv = tsv(:);
tdv = tdv(:);
t0v = t0v(:);

% The purpose: Integrate the difference of the sinc filter minus the ideal
% filter for different taus. How big are the differences? How do they
% change with beta?

% Here we also include the effect from a "bioturbative" process from
% uniform smoothing over a time duration taud

% function out = integrate_powerlaw(tbv,tsv,beta)
% tbv is a vector of tau_b values
% tsv is a vector of tau_s values
% t0v is a vector of record lengths
% beta is the negative spectral slope

if max(size(t0v))==1,t0v=t0v*ones(length(tbv));end
if max(size(tdv))==1,tdv=tdv*ones(length(tbv));end

fun1 = @(x,taub,taud,b) x.^-b.*(1-sinc(taub.*x).*sinc(taud.*x)).^2;
fun2 = @(x,taub,taud,b) x.^-b.*(sinc(taub.*x).*sinc(taud.*x)).^2;
%fun3 = @(x,b) x.^-b;

out = [];

for ii = 1:length(tbv)
    taus = tsv(ii);
    taub = tbv(ii);
    taud = tdv(ii);
    tau0 = t0v(ii);
    out(ii) = (integral(@(x)fun1(x,taub,taud,bet),1/tau0,1./(2*taus))+...
              integral(@(x)fun2(x,taub,taud,bet),1./(2*taus),Inf))./...
              integral(@(x)fun2(x,taub,taud,bet),1/tau0,Inf);
end
