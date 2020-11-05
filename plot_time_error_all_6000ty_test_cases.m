
% Some constants
L = 40;
t0 = 21000;
tyMIN = 10;
tyMAX = 20000;
txMIN = 10;
txMAX = 10000;

s1 = linspace(tyMIN,tyMAX,L);
s2 = linspace(txMIN,txMAX,L);
set1 = repmat(s1,L,1);
set2 = repmat(s2,L,1)';
s11 = set1(:);
s22 = set2(:);

%% Case 1: vary D

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\Delta';

tyv = s11;
txv = 4000;
t0v = t0;
tav = 0;
tdv = 0;
dv  = s22;

b = .1;
out5 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
%out15 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);


save('lag_denom_delt_all_10000','out5','s1','s2','s1name','s2name')


%% Case 2: vary D with different taua


% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\Delta';

tyv = s11;
txv = 4000;
t0v = t0;
tav = 1000;
tdv = 0;
dv  = s22;

b = .5;
out5 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_error_powerlaw3(t0v,tav,tdv,txv,tyv,dv,b);


%% Plot!

close all


files = {'lag_denom_delt_all_10000.mat'
         };

     
for ii = 1:length(files)

load(files{ii})

addpath ../../export_fig
outvec = {'5' '15'};
bvec = {'0.5','1.5'};
clabv = {[[0.05,.1,.2,.3,0.4:.1:2]]
         [0.01,0.05,.1,.2,.3,0.4,0.8,2,5,10]
         %[0.01,0.1:0.2:2]
         };


for jj = 1:length(outvec)
    plout = nan(L,L);
    figure()    
    set(gcf,'color','w','position',[440   475   403   323])

    plout(:) = eval(['out' outvec{jj}]);
    hold all
    pcolor(s1,s2,log10(plout+eps)),shading interp
    caxis([-1.5,2])
    [C,h] = contour(s1,s2,plout,clabv{jj},'k');
    colormap(flipud(hot(100)))
    clabel(C,h,'LabelSpacing',500)
    ylabel([s2name ' (years)'],'fontsize',12)
    xlabel([s1name ' (years)'],'fontsize',12)
    grid on
    axis tight
    axis equal
    %text(400,200,['\beta = ' bvec{ii}],'fontsize',14,'fontweight','bold')
    set(gca,'XTick',0:500:10000)
    set(gca,'XTickLabelRotation',45,'fontsize',12)
    %set(gca,'YTick',get(gca,'XTick'))
    plot([4000 4000],[0,4000],'k--','LineWidth',1.5)
    if ii ==2 || ii ==4
        plot([1000 1000],[0,4000],'k-','LineWidth',1.5)
    end
    %if ii == 1
    %    plot(s1,abs(4000-s1)/2,'k-','LineWidth',1.5)
    %end
    %export_fig('-png','-r200',['Figs/lag_denom_errors_wrt_beta_6000_' bvec{jj} '_' files{ii} '_tau0_' num2str(t0)])    
    
end

end

