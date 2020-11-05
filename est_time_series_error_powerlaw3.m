function out = est_time_series_error_powerlaw3(t0v,tav,tsv,tyv,dv,b)

% This version (3) puts the lag in the denominator!
% No mechanism for chron uncertainty right now.

% t0v is a vector of record lengths (Necessary for power law spectra)
% tav is a vector of archive smoothing time scales
% tdv is a vector of chronological uncertainty time scales
% tsv is a vector of sample spacing
% tyv is a vector observational average time scales
% dv  is a vector of time offsets between the target and observations
% b   is the negative spectral slope of the underlying climate signal. Just
%     one value at a time!

L = max([length(t0v),length(tav),length(tyv),length(dv)]);

% If any variables are scalars, make them into constant vectors
if max(size(t0v))==1,t0v=t0v*ones(L);end
if max(size(tav))==1,tav=tav*ones(L);end
if max(size(tyv))==1,tyv=tyv*ones(L);end
if max(size(dv))==1,dv=dv*ones(L);end

% Integrand in the numerator of the error fraction - Part that overlaps
% with 1 in the Heaviside ideal transfer function
funi1 =  @(taua,tauy,b,d,x)...
        x.^-b.*abs(1-exp(-2.*pi.*1i.*x.*d)...
        .*sinc(taua.*x).*sinc(tauy.*x)).^2;

% Integrand in the numerator of the error fraction - Part that overlaps
% with 0 in the Heaviside ideal transfer function
funi0 =  @(taua,tauy,b,d,x)...
        x.^-b.*abs( -exp(-2.*pi.*1i.*x.*d)...
        .*sinc(taua.*x).*sinc(tauy.*x)).^2;
    
    % Numerator of the error function with no age uncertainty (taud=0)
fun =   @(tau0,taua,taus,tauy,b,d)...
        integral(@(x)funi1(taua,tauy,b,d,x),1/tau0,1./(2*taus))+ ...
        integral(@(x)funi0(taua,tauy,b,d,x),1./(2*taus),Inf);

% Integrand in the denominator of the fraction
fundi = @(b,x)...
    x.^-b;

% Denominator of the error function
fund =  @(tau0,taus,b)...
        integral(@(x)fundi(b,x),1/tau0,1./(2*taus));

% Initialize output    
out = [];

parfor ii = 1:L
    taua = tav(ii);
    taus = tsv(ii);
    tauy = tyv(ii);
    tau0 = t0v(ii);
    d    = dv(ii);
    
    out(ii) = sqrt( fun(tau0,taua,taus,tauy,b,d)...
            ./fund(tau0,taus,b) );
end


%{
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

%}