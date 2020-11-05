% This version uses est_error_powerlaw3 (puts lag in denominator)

% Some constants
L = 50;
t0 = 21000;
tyMIN = 10;
tyMAX = 4000;
txMIN = 10;
txMAX = 4000;

%% Case 1: vary ty and tx


s1 = linspace(tyMIN,tyMAX,L);
s2 = linspace(txMIN,txMAX,L);
set1 = repmat(s1,L,1);
set2 = repmat(s2,L,1)';
s11 = set1(:);
s22 = set2(:);

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

b = .5;
out5 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);


%%
addpath ../../export_fig
outvec = {'5' '15'};
bvec = {'0.5','1.5'};

clabv = {[.1,.2,.3,0.4,0.8,1,2,5,10]
         [0.02,0.05,.1,.2,.3,0.4,0.8,2,5,10]
         };

for ii = 1:length(outvec)
    plout = nan(L,L);
    figure()    
    set(gcf,'color','w','position',[440   518   403   280])

    plout(:) = eval(['out' outvec{ii}]);
    hold all
    pcolor(s1,s2,log10(plout+eps)),shading interp
    caxis([-1.5,2])
    [C,h] = contour(s1,s2,plout,clabv{ii},'k');
    colormap(flipud(hot))
    ylabel([s2name ' (years)'],'fontsize',12)
    xlabel([s1name ' (years)'],'fontsize',12)
    grid on
    axis tight
    axis square
    %text(400,200,['\beta = ' bvec{ii}],'fontsize',14,'fontweight','bold')
    set(gca,'XTick',get(gca,'YTick'))
    set(gca,'XTickLabelRotation',45,'fontsize',12)
    set(gca,'YTick',get(gca,'XTick'))
    %plot([1000 1000],[0,4000],'k-','LineWidth',1.5)
    clabel(C,h)%,'manual')

    s2l = s2;
    s2l(s2l<1000) = [];
    
    plot(sqrt((s2.^2-tav^2)),s2,'k--')
    export_fig('-png','-r200',['Figs/lag_denom_errors_wrt_beta_' char(outvec(ii)) '_tau0_' num2str(t0)])    
    
end


%% Case 1: vary ty and tx with ta


s1 = linspace(tyMIN,tyMAX,L);
s2 = linspace(txMIN,txMAX,L);
set1 = repmat(s1,L,1);
set2 = repmat(s2,L,1)';
s11 = set1(:);
s22 = set2(:);

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\tau_x';

txv = s22;
tyv = s11;
t0v = t0;
tav = 1000;
tdv = 0;
dv  = 0;

b = .5;
out5 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);


%%
addpath ../../export_fig
outvec = {'5' '15'};
bvec = {'0.5','1.5'};

clabv = {[.1,.2,.3,0.4,0.8,1,2,5,10]
         [0.02,0.05,.1,.2,.3,0.4,0.8,2,5,10]
         };


for ii = 1:length(outvec)
    plout = nan(L,L);
    figure()    
    set(gcf,'color','w','position',[440   518   403   280])

    plout(:) = eval(['out' outvec{ii}]);
    hold all
    pcolor(s1,s2,log10(plout+eps)),shading interp
    caxis([-1.5,2])
    [C,h] = contour(s1,s2,plout,clabv{ii},'k');
    colormap(flipud(hot))
    ylabel([s2name ' (years)'],'fontsize',12)
    xlabel([s1name ' (years)'],'fontsize',12)
    grid on
    axis tight
    axis square
    %text(400,200,['\beta = ' bvec{ii}],'fontsize',14,'fontweight','bold')
    set(gca,'XTick',get(gca,'YTick'))
    set(gca,'XTickLabelRotation',45,'fontsize',12)
    set(gca,'YTick',get(gca,'XTick'))
    plot([1000 1000],[0,4000],'k-','LineWidth',1.5)
    clabel(C,h)%,'manual')
    
    s2l = s2;
    s2l(s2l<1000) = [];
    
    plot(sqrt((s2l.^2-1000^2)),s2l,'k--')
    export_fig('-png','-r200',['Figs/lag_denom_errors_wrt_beta_taua_' char(outvec(ii)) '_tau0_' num2str(t0)])    
    
end

