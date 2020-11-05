% Some constants
t0 = 21000;
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

%L = length(s1)*length(s2);

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\tau_x';

txv = s22;
tyv = s11;
t0v = t0;
tav = 0;
tdv = 0;
dv  = 0;


f = tdfread('/Users/dan/Dropbox (MIT)/2017-2018/Sampling/NGRIP_oxygen_isotope_50.tab');

age = 1000 * f.Age_0x5Bka_BP0x5D;
d18 = f.d18O_H2O;

%out = sample_ngrip(tav,txv,tyv);
%out = sample_ngrip(dat,50,t0v,tav,txv,tyv,100);
out = sample_ts(d18,50,t0v,tav,txv,tyv,100);


%%
addpath ../../export_fig

plout = nan(length(s1),length(s2));
figure()    
set(gcf,'color','w','position',[440   518   403   280])

plout(:) = out;
hold all
pcolor(s1,s2,log10(plout+eps)),shading interp
caxis([-1.5,2])
[C,h] = contour(s1,s2,plout,[0.02,0.05,.1,.2,.3,0.4,0.8,1,2,5,10],'k');
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
export_fig('-png','-r200','Figs/lag_denom_ngrip_errors')    
    
