function out = est_error_powerlaw(t0v,tav,tdv,txv,tyv,dv,b)

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
        x.^-b.*abs(sinc(taux.*x)-exp(-2.*pi.*1i.*x.*d).*sinc(taua.*x).*sinc(tauy.*x)).^2;
%        x.^-b.*(sinc(taux.*x)-sinc(taua.*x).*sinc(tauy.*x)).^2;

% Integrand in the denominator of the fraction
% ...WHEN THE DENOMINATOR IS THE OB VARIANCE
% fundi = @(taua,tauy,b,x) x.^-b.*(sinc(taua.*x).*sinc(tauy.*x)).^2;

% Integrand in the denominator of the fraction
% ...WHEN THE DENOMINATOR IS THE OB VARIANCE
fundi = @(taux,b,x) x.^-b.*(sinc(taux.*x)).^2;

% Numerator of the error function (no age uncertainty, taud=0)
fun =   @(tau0,taua,taux,tauy,b,d)...
        integral(@(x)funi(taua,taux,tauy,b,d,x),1/tau0,Inf);

% Denominator of the error function
% ...WHEN THE DENOMINATOR IS THE OB VARIANCE
% fund =  @(tau0,taua,tauy,b)...
%        integral(@(x)fundi(taua,tauy,b,x),1/tau0,Inf);

% Denominator of the error function
fund =  @(tau0,taux,b)...
        integral(@(x)fundi(taux,b,x),1/tau0,Inf);

    % Include effects from age uncertainty
% Here d is the mean age offset, taud is the standard deviation of the age
% uncertainty, and dd is the Delta that we are integrating over.
%fun_a = @(tau0,taua,taud,taux,tauy,b,d)integral(@(dd)fun(tau0,taua,taux,tauy,b,dd),-1,1,'Arrayvalued',true);

fun_a = @(tau0,taua,taud,taux,tauy,b,d)integral(@(dd)...
        1/((taud)*sqrt(2*pi)) * exp(-(d-dd).^2./(2*(taud+eps).^2)).* ...   % Gaussian PDF
        fun(tau0,taua,taux,tauy,b,dd),-Inf,Inf,'Arrayvalued',true,'AbsTol',1e-4);

out = [];

for ii = 1:L
    taua = tav(ii);
    taux = txv(ii);
    tauy = tyv(ii);
    tau0 = t0v(ii);
    d    = dv(ii);
    taud = tdv(ii);
    
    % If there is age uncertainty, compute the expectation over it.
    if ~~taud
        out(ii) = fun_a(tau0,taua,taud,taux,tauy,b,d)...
                ./fund(tau0,taux,b);
%                ./fund(tau0,taua,tauy,b);
    
    % ...otherwise use the standard formula
    else
        out(ii) = fun(tau0,taua,taux,tauy,b,d)...
                ./fund(tau0,taux,b);
%                ./fund(tau0,taua,tauy,b);              
    end
    %ii
end
