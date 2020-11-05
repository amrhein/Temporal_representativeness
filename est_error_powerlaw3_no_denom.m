function out = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b)

% This version (3) puts the lag in the denominator!
% t0v is a vector of record lengths (Necessary for power law spectra)
% tav is a vector of archive smoothing time scales
% tdv is a vector of chronological uncertainty time scales
% txv is a vector of target average time scales
% tyv is a vector observational average time scales
% dv  is a vector of time offsets between the target and observations
% b   is the negative spectral slope of the underlying climate signal. Just
%     one value at a time!

L = max([length(t0v),length(tav),length(tdv),length(txv),length(tyv),length(dv)]);

% If any variables are scalars, make them into constant vectors
if max(size(t0v))==1,t0v=t0v*ones(L);end
if max(size(tav))==1,tav=tav*ones(L);end
if max(size(tdv))==1,tdv=tdv*ones(L);end
if max(size(txv))==1,txv=txv*ones(L);end
if max(size(tyv))==1,tyv=tyv*ones(L);end
if max(size(dv))==1,dv=dv*ones(L);end

% Integrand in the numerator of the error fraction
funi =  @(taua,taux,tauy,b,d,x)...
        x.^-b.*abs(sinc(taux.*x)-exp(-2.*pi.*1i.*x.*d)...
        .*sinc(taua.*x).*sinc(tauy.*x)).^2;

% Numerator of the error function with no age uncertainty (taud=0)
fun =   @(tau0,taua,taux,tauy,b,d)...
        integral(@(x)funi(taua,taux,tauy,b,d,x),0,Inf);

% Numerator including effects from age uncertainty
% Here d is the mean age offset, taud is the standard deviation of the age
% uncertainty, and dd is the Delta that we are integrating over.

fun_a = @(tau0,taua,taud,taux,tauy,b,d)...
    integral2(...
    @(dd,x)...
    1/((taud)*sqrt(2*pi)) * exp(-(d-dd).^2./(2*(taud+eps).^2)).*...
    x.^-b.*abs(sinc(taux.*x)-...
    exp(-2.*pi.*1i.*x.*(d-dd)).*sinc(taua.*x).*sinc(tauy.*x)).^2,...
    -Inf,Inf,1/tau0,Inf,...
    'method','iterated','AbsTol',1e-12,'RelTol',1e-4...
    );
  
% Integrand in the denominator of the fraction
fundi = @(tau0,taux,b,x)...
    x.^-b.*abs( (1-exp(-2.*pi.*1i.*x.*tau0)).*(sinc(taux.*x))).^2;

% Denominator of the error function
fund =  @(tau0,taux,b)...
        integral(@(x)fundi(tau0,taux,b,x),0,Inf);

% Initialize output    
out = [];

parfor ii = 1:L
    
% NO DENOM
    
    taua = tav(ii);
    taux = txv(ii);
    tauy = tyv(ii);
    tau0 = t0v(ii);
    d    = dv(ii);
    taud = tdv(ii);
    
    % If there is age uncertainty, use fun_a...
    if ~~taud
        ii
        out(ii) = sqrt( fun_a(tau0,taua,taud,taux,tauy,b,d)...
                 );

    % ...otherwise use fun
    else
        out(ii) = sqrt( fun(tau0,taua,taux,tauy,b,d)...
                 );
    end
end
