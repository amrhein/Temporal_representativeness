% Some constants
tyMIN = 50;
tyMAX = 4000;
txMIN = 50;
txMAX = 4000;


s1 = tyMIN:50:tyMAX;
s2 = txMIN:50:txMAX;
set1 = repmat(s1,length(s2),1);
set2 = repmat(s2,length(s1),1)';
s11 = set1(:);
s22 = set2(:);

cn = dsp.ColoredNoise(1.5,4918);
dat = cn();



% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\tau_x';

t0 = 21000;

txv = s22;
tyv = s11;
t0v = t0;
tav = 0;
tdv = 0;
dv  = 0;

out = sample_ts(dat,50,t0v,tav,txv,tyv,100);


addpath ../../export_fig
clabv = {[0.05,.2,0.4,.6,0.8,2,5,10]
         [0.01,0.05,.1,.2,.3,0.4,0.8,2,5,10]
         };

plout = nan(length(s1),length(s2));
figure()    
set(gcf,'color','w','position',[440   518   403   280])

plout(:) = out;
hold all
%pcolor(s1,s2,log10(plout+eps)),shading interp
caxis([-1.5,2])
[C,h] = contour(s1,s2,plout,[0.01,0.05,.1,.2,.3,0.4,0.8,2,5,10],'k');
colormap(flipud(hot))
ylabel([s2name ' (years)'],'fontsize',12)
xlabel([s1name ' (years)'],'fontsize',12)
grid on
axis tight
axis square
set(gca,'XTick',get(gca,'YTick'))
set(gca,'XTickLabelRotation',45,'fontsize',12)
set(gca,'YTick',get(gca,'XTick'))
%plot([1000 1000],[0,4000],'k-','LineWidth',1.5)
%clabel(C,h,'manual')
clabel(C,h)
%export_fig('-png','-r200',['Figs/ngrip_errors'])    
    
%% Sanity check

cn = dsp.ColoredNoise(1.5,90000);
dat = cn();


tx = 4000;
ty = 100;
t0 = 21000;
DT = 50;
nx = round(tx/DT);
ny = round(ty/DT);
wx = 1/nx * ones(1,nx);
wy = 1/ny * ones(1,ny);
N = length(dat);

s = 1000;

x  = conv(dat,wx,'same');

    xa = dat;

y  = conv(xa, wy,'same');

xs = x(s:(end-s));
ys = y(s:(end-s));

% Compute lag tau0 autocovariance
d0 = round(t0/DT);
early = xs(d0:end);
late = xs(1:(end-d0+1));
%av = 1/(N-s-d0-1)*sum((early-late-mean(early-late)).^2)
av = var(early-late)
%1/(N-1)*sum((xs-ys - mean(xs-ys)).^2)
var(xs-ys)
%%

sa = dsp.SpectrumAnalyzer('SpectrumType','Power density', ...
    'OverlapPercent',50,'Window','Hamming', ...
    'SpectralAverages',50,'PlotAsTwoSidedSpectrum',false, ...
    'FrequencyScale','log','YLimits',[-50 20]);

tic
while toc < 30
pink = cn();
sa(pink);
end
