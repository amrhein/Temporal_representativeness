
% Some constants
L = 50;
t0 = 21000;
tyMIN = 10;
tyMAX = 4000;
txMIN = 10;
txMAX = 4000;

%% Case 1: vary D

s1 = linspace(tyMIN,tyMAX,L);

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\tau_x';

tyvec = [500,1000,4000];

out5 = [];
out15 = [];

for ii = 1:length(tyvec)

    tyv = tyvec(ii)
    txv = 4000;
    t0v = t0;
    tav = 0;
    tdv = 0;
    dv  = s1;

    b = .5;
    out5(:,ii) = est_error_powerlaw(t0v,tav,tdv,txv,tyv,dv,b);

    b = 1.5;
    out15(:,ii) = est_error_powerlaw(t0v,tav,tdv,txv,tyv,dv,b);
end

save('delt_1','out5','out15','s1')



%% Case 2: vary D with different taua
s1 = linspace(tyMIN,tyMAX,L);

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\tau_x';

tyvec = [500,1000,4000];

out5 = [];
out15 = [];

for ii = 1:length(tyvec)

    tyv = tyvec(ii)
    txv = 4000;
    t0v = t0;
    tav = 1000;
    tdv = 0;
    dv  = s1;

    b = .5;
    out5(:,ii) = est_error_powerlaw(t0v,tav,tdv,txv,tyv,dv,b);

    b = 1.5;
    out15(:,ii) = est_error_powerlaw(t0v,tav,tdv,txv,tyv,dv,b);
end

save('delt_1a','out5','out15','s1')
%% Case 3: vary tdv

s1 = linspace(tyMIN,tyMAX,L);

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\tau_x';

tyvec = [500,1000,4000];

out5 = [];
out15 = [];

for ii = 1:length(tyvec)

    tyv = tyvec(ii)
    txv = 4000;
    t0v = t0;
    tav = 0;
    tdv = s1;
    dv  = 0;

    b = .5;
    out5(:,ii) = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

    b = 1.5;
    out15(:,ii) = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);
end

save('sigdelt_2','out5','out15','s1')

%{
%%
addpath ../../export_fig

figure()
set(gcf,'color','w')
plot(s1,out5)
xlabel('\sigma_\Delta (years)','fontsize',12)
legend('500','1000','4000')
export_fig('-pdf','Figs/sigdelt5')

figure()
set(gcf,'color','w')
plot(s1,out15)
xlabel('\sigma_\Delta (years)','fontsize',12)
legend('500','1000','4000')
export_fig('-pdf','Figs/sigdelt15')
%}

%% Case 3: vary tdv with different taua


s1 = linspace(tyMIN,tyMAX,L);

% Varies along horizontal plot axis
s1name = '\tau_y';
% Varies along vertical plot axis
s2name = '\tau_x';

tyvec = [500,1000,4000];

out5 = [];
out15 = [];

for ii = 1:length(tyvec)

    tyv = tyvec(ii)
    txv = 4000;
    t0v = t0;
    tav = 1000;
    tdv = s1;
    dv  = 0;

    b = .5;
    out5(:,ii) = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);

    b = 1.5;
    out15(:,ii) = est_error_powerlaw2(t0v,tav,tdv,txv,tyv,dv,b);
end

save('sigdelt_2a','out5','out15','s1')


%% Plot!

close all

load delt_1
addpath ../../export_fig

figure()
set(gcf,'color','w')
plot(s1,out5)
xlabel('\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('Delta .5')
export_fig('-pdf','Figs/delt5')

figure()
set(gcf,'color','w')
plot(s1,out15)
xlabel('\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('Delta 1.5')
export_fig('-pdf','Figs/delt15')

load delt_1a
addpath ../../export_fig

figure()
set(gcf,'color','w')
plot(s1,out5)
xlabel('\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('Delta .5 a')
export_fig('-pdf','Figs/delt5a')

figure()
set(gcf,'color','w')
plot(s1,out15)
xlabel('\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('Delta 1.5 a')
export_fig('-pdf','Figs/delt15a')

load sigdelt_2
addpath ../../export_fig

figure()
set(gcf,'color','w')
plot(s1,out5)
xlabel('\sigma_\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('SigDelta .5')
export_fig('-pdf','Figs/sigdelt5')

figure()
set(gcf,'color','w')
plot(s1,out15)
xlabel('\sigma_\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('SigDelta 1.5')
export_fig('-pdf','Figs/sigdelt15')


load sigdelt_2a
addpath ../../export_fig

figure()
set(gcf,'color','w')
plot(s1,out5)
xlabel('\sigma_\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('SigDelta .5 a')
export_fig('-pdf','Figs/sigdelt5a')

figure()
set(gcf,'color','w')
plot(s1,out15)
xlabel('\sigma_\Delta (years)','fontsize',12)
legend('500','1000','4000')
title('SigDelta 1.5a')
export_fig('-pdf','Figs/sigdelt15a')

