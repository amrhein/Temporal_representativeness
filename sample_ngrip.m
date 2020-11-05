function out = sample_ngrip(dat,DT,t0v,tav,txv,tyv)

% t0v is a vector of record lengths (Necessary for power law spectra)
% tav is a vector of archive smoothing time scales
% tdv is a vector of chronological uncertainty time scales
% txv is a vector of target average time scales
% tyv is a vector observational average time scales
% dv  is a vector of time offsets between the target and observations
% b   is the negative spectral slope of the underlying climate signal. Just
%     one value at a time!

% Samples the NGRIP d18O record

f = tdfread('/Users/dan/Dropbox (MIT)/2017-2018/Sampling/NGRIP_oxygen_isotope_50.tab');

time2 = 1950-(1000 * flipud(f.Age_0x5Bka_BP0x5D));
time = time2(1:2:end);

d182 = f.d18O_H2O;
d18 = d182(1:2:end);

N = length(d18);
DT = round(median(diff(time))); % 50


L = max([length(tav),length(txv),length(tyv)])

% If any variables are scalars, make them into constant vectors
if max(size(tav))==1,tav=tav*ones(L);end
if max(size(txv))==1,txv=txv*ones(L);end
if max(size(tyv))==1,tyv=tyv*ones(L);end

out = nan(1,L);

s = 100;

for ii = 1:L
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

    x  = conv(d18,wx,'same');

    if ~isempty(wa)
        xa = conv(d18,wa,'same');
    else
        xa = d18;
    end
    
    y  = conv(xa, wy,'same');

    xs = x(s:(end-s));
    ys = y(s:(end-s));

    %out(ii) = 1/(N-1)*sum((xs-ys).^2);
        d0 = round(t0/DT);
    early = xs(d0:end);
    late = xs(1:(end-d0+1));

    out(ii) = var(xs-ys)/var(early-late);
    var(early-late)

end