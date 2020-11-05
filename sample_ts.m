function out = sample_ts(dat,DT,t0v,tav,txv,tyv,s)

% t0v is a vector of record lengths (Necessary for power law spectra)
% tav is a vector of archive smoothing time scales
% tdv is a vector of chronological uncertainty time scales
% txv is a vector of target average time scales
% tyv is a vector observational average time scales
% dv  is a vector of time offsets between the target and observations
% b   is the negative spectral slope of the underlying climate signal. Just
%     one value at a time!

% Samples the NGRIP d18O record

%f = tdfread('/Users/dan/Dropbox (MIT)/2017-2018/Sampling/NGRIP_oxygen_isotope_50.tab');

%dat = f.d18O_H2O;

N = length(dat);
L = max([length(tav),length(txv),length(tyv)]);

% If any variables are scalars, make them into constant vectors
if max(size(t0v))==1,t0v=t0v*ones(L);end
if max(size(tav))==1,tav=tav*ones(L);end
if max(size(txv))==1,txv=txv*ones(L);end
if max(size(tyv))==1,tyv=tyv*ones(L);end

out = nan(1,L);

%s = 100;

parfor ii = 1:L
    ta = tav(ii);
    tx = txv(ii);
    ty = tyv(ii);
    t0 = t0v(ii);

    na = round(ta/DT);
    nx = round(tx/DT);
    ny = round(ty/DT);
    wa = 1/na * ones(1,na);
    wx = 1/nx * ones(1,nx);
    wy = 1/ny * ones(1,ny);

    x  = conv(dat,wx,'same');
    
    if ~isempty(wa)
        xa = conv(dat,wa,'same');
    else
        xa = dat;
    end
    
    y  = conv(xa, wy,'same');

    xs = x(s:(end-s));
    ys = y(s:(end-s));

    % Compute lag tau0 autocovariance
    %av = 1/(N-s-d0-1)*sum((xs(d0:end)-xs(1:(end-d0+1))).^2);
    %out(ii) = 1/(N-1)*sum((xs-ys - mean(xs-ys)).^2)/av;
    d0 = round(t0/DT);
    early = xs(d0:end);
    late = xs(1:(end-d0+1));

    out(ii) = sqrt(var(xs-ys)/var(early-late));
    var(early-late)
end
