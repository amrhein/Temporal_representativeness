
% Some constants
L = 30;
t0 = 21000;
tyMIN = 10;
tyMAX = 4000;
txMIN = 10;
txMAX = 4000;

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

b = .5;
out5 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);


save('delt_all_1','out5','out15','s1','s2','s1name','s2name')


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
out5 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);


save('delt_all_1a','out5','out15','s1','s2','s1name','s2name')


%% Case 3: vary tdv

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\sigma_\Delta';

tyv = s11;
txv = 4000;
t0v = t0;
tav = 0;
tdv = s22;
dv  = 0;

b = .5;
out5 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

save('sigdelt_all_1','out5','out15','s1','s2','s1name','s2name')
%% Case 3: vary tdv with different taua

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\sigma_\Delta';

tyv = s11;
txv = 4000;
t0v = t0;
tav = 1000;
tdv = s22;
dv  = 0;

b = .5;
out5 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

b = 1.5;
out15 = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

save('sigdelt_all_1a','out5','out15','s1','s2','s1name','s2name')

%% Plot!

close all


files = {'delt_all_1.mat'
         'delt_all_1a.mat'
         'sigdelt_all_1.mat'
         'sigdelt_all_1a.mat'
         }

     
for ii = 1:length(files)

load(files{ii})

addpath ../../export_fig
outvec = {'5' '15'};
bvec = {'0.5','1.5'};
clabv = {[0.05,.2,0.4,.6,0.8,2,5,10]
         [0.01,0.05,.1,.2,.3,0.4,0.8,2,5,10]
         %[0.01,0.1:0.2:2]
         };


for jj = 1:length(outvec)
    plout = nan(L,L);
    figure()    
    set(gcf,'color','w','position',[440   518   403   280])

    plout(:) = eval(['out' outvec{jj}]);
    hold all
    pcolor(s1,s2,log10(plout+eps)),shading interp
    caxis([-1.5,2])
    [C,h] = contour(s1,s2,plout,clabv{jj},'k');
    colormap(flipud(hot))
    clabel(C,h,'LabelSpacing',500)
    ylabel([s2name ' (years)'],'fontsize',12)
    xlabel([s1name ' (years)'],'fontsize',12)
    grid on
    axis tight
    axis square
    %text(400,200,['\beta = ' bvec{ii}],'fontsize',14,'fontweight','bold')
    set(gca,'XTick',get(gca,'YTick'))
    set(gca,'XTickLabelRotation',45,'fontsize',12)
    set(gca,'YTick',get(gca,'XTick'))
    export_fig('-png','-r200',['Figs/errors_wrt_beta_' bvec{jj} '_' files{ii} '_tau0_' num2str(t0)])    
    
end

end
