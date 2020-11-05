% Some constants
L = 50;
t0 = 21000;
tyMIN = 10;
tyMAX = 4000;
tsMIN = 10;
tsMAX = 4000;

%% Case 1: vary ty and tx


s1 = linspace(tyMIN,tyMAX,L);
s2 = linspace(tsMIN,tsMAX,L);
set1 = repmat(s1,L,1);
set2 = repmat(s2,L,1)';
s11 = set1(:);
s22 = set2(:);

% Varies along horizontal plot axis
s1name = '\tau^i_y';
% Varies along vertical plot axis
s2name = '\tau^i_s';

tsv = s22;
tyv = s11;
t0v = t0;
tav = 0;
dv  = 0;

b = 0.5;
out5 = est_time_series_error_powerlaw3(t0v,tav,tsv,tyv,dv,b);
%out5 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_time_series_error_powerlaw3(t0v,tav,tsv,tyv,dv,b);
%out15 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);


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
    s2l(s2l<tav) = [];
    
    %plot(s2,.89*s2,'k--')
    %plot(sqrt(1.5)*s2,s2,'k--')
    %xlim([10,4000])
    
    export_fig('-png','-r200',['Figs/time_series_errors_' char(outvec(ii)) '_tau0_' num2str(t0)])    
    
end


%% Case 2: vary ty and tx with ta

s1 = linspace(tyMIN,tyMAX,L);
s2 = linspace(tsMIN,tsMAX,L);
set1 = repmat(s1,L,1);
set2 = repmat(s2,L,1)';
s11 = set1(:);
s22 = set2(:);

% Varies along horizontal plot axis
s1name = '\tau^i_y';
% Varies along vertical plot axis
s2name = '\tau^i_s';

tsv = s22;
tyv = s11;
t0v = t0;
tav = 1000;
dv  = 0;

b = .5;
out5 = est_time_series_error_powerlaw3(t0v,tav,tsv,tyv,dv,b);
%out5 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_time_series_error_powerlaw3(t0v,tav,tsv,tyv,dv,b);
%out15 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);


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
    
    %s2l = s2;
    %s2l(s2l<1000) = [];
    
    %plot(sqrt((1.50*s2.^2-tav^2)),s2,'k--')
    %xlim([10,4000])
    
    export_fig('-png','-r200',['Figs/time_series_errors_taua' char(outvec(ii)) '_tau0_' num2str(t0)])    
    
end

    

