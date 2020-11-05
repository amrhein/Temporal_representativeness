% Fix pmtm plot
% Make a squared misfit plot?

clear

addpath('/Users/dan/Dropbox (MIT)/2011-2012/820/PS3/')

f = tdfread('/Users/dan/Dropbox (MIT)/2017-2018/Sampling/NGRIP_oxygen_isotope_50.tab');

time2 = 1950-(1000 * flipud(f.Age_0x5Bka_BP0x5D));
time = time2(1:2:end);

d182 = flipud(f.d18O_H2O);
d18 = d182(1:2:end);

N = length(d18);
dt = round(median(diff(time))); % 50

n = length(d18);

%% Time series smoothing stuff
s = 40;

ta = 0;
tx = 4000;
ty = 1000;

na = round(ta/dt);
nx = round(tx/dt);
ny = round(ty/dt);
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

if s>0
    xs = x(s:(end-s));
    ys = y(s:(end-s));
    time_s = time(s:(end-s));
else
    xs = x;
    yx = y;
    time_s = time;
end


%% Wavelet stuff


pad = 1;      % pad the time series with zeroes (recommended)
dj = 0.1;    % this will do 4 sub-octaves per octave
s0 = 4*dt;    % this says start at a scale of 6 months
j1 = 7/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72;  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(d18,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Wavelet transform:
%[wave,period,scale,coi] = wavelet(d18,dt);%,pad,dj,s0,j1,mother);
%power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Global wavelet spectrum & significance levels:
%global_ws = variance*(sum(power')/n);   % time-average over all times
%dof = n - scale;  % the -scale corrects for padding at edges

freq = 1./period;

xl = [time(1), time(end)];  % plotting range

%% Plotting

xTickL = {'120','100','80','60','40','20','0'};

%--- Plot time series

close all
set(gcf,'color','w','position',[322         171        1071         627])
subplot('position',[0.280    0.7571    0.6500    0.2000])

yyaxis left
hold on
plot(time,d18,'color',.7*[1 1 1])
plot(time_s,xs,'k-')
plot(time_s,ys,'r-')
set(gca,'XLim',xl(:))
xlabel('Time (kya)')
ylabel('\delta^{18}O (VSMOW)')
set(gca,'XTicklabel',xTickL)
set(gca,'Ylim',[-50 -32])
set(gca,'YTick',-44:2:34)
set(gca,'YColor','k')
ylabh = get(gca,'yLabel');
set(ylabh,'Position',get(ylabh,'Position') + [0 .2 0])


yyaxis right
hold all
plot(time_s,(ys-xs),'k-','linewidth',1)
plot([min(time),max(time)],[0 0],'color',.8*[1 1 1])

set(gca,'XLim',xl(:))
set(gca,'Ylim',[-4 22])
set(gca,'YTick',[-2 0 2])
set(gca,'YColor','k')
ylabel('\theta')


columnlegend(4,{'NGRIP \delta^{18}O of water','4000-yr running mean','1000-yr running mean','\theta (difference)'},'location','north','box','off')


%%  Squared misfit subplot


te = 1/2379*sum((xs-ys).^2);

subplot('position',[0.280    0.6316    0.6500    0.0558])

hold on
plot(time_s,(xs-ys).^2)
plot([min(time) max(time)],[te te],'r-')
%area(time_s,(xs-ys).^2,'EdgeColor','none','FaceColor',[0 0.4470 0.7410])
%legend('\theta^2','<\theta^2>','box','off','location','northwest')
columnlegend(2,{'\theta^2','<\theta^2>'},'box','off','location','northwest')
%columnlegend(2,{'t','t'},'box','off','location','northwest')
set(gca,'XLim',xl(:))
xlabel('Time (kya)')
ylabel('(\delta^{18}O)^2')



%% Wavelet subplot
hold off

%--- Contour plot wavelet power spectrum
subplot('position',[0.280    0.3317    0.6500    0.2265])
levels = [eps,.125,.25,.5,1,2,4,8,16,32,64,128,256,1000] ;
%Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
Yticks = [500 1000 2000 4000 8000 16000];
ytl = {'1/500' '1/1000' '1/2000' '1/4000' '1/8000' '1/16000'};

contourf(time,log2(period),log2(power),log2(levels),'edgecolor','none');  %*** or use 'contourfill'
%imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
colormap((parula))
xlabel('Time (kya)')
ylabel('Frequency (yrs^{-1})')
set(gca,'XLim',xl(:))
set(gca,'YLim',log2([min(period),max(period)]), ...
	'YDir','reverse', ...
    'XTicklabel',xTickL,...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',ytl)
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
%contour(time,log2(period),sig95,[-99,1],'k');
hold on
% cone-of-influence, anything "below" is dubious
plot(time,log2(coi),'k')
h = colorbar('eastoutside')
set(get(h,'ylabel'),'string','Power, (\delta^{18}O)^2')
set(h,'YTick',log2([1 10 100 1000 10000 100000]),...
    'YTicklabel',{'10^0','10^1','10^2','10^3','10^4','10^5'})

set(gca,'position',[0.280    0.3317    0.6500    0.2265])
plot([min(time),max(time)],log2([2257,2257]),'color',0.5*[1 1 1],'LineWidth',1)
plot([min(time),max(time)],log2([5298,5298]),'color',0.5*[1 1 1],'LineWidth',1)
hold off
caxis([-3,8])

%% pmtm subplot

[p,f] = pmtm(d18,3,length(d18),1/50);

subplot('position',[0.0663    0.3288    0.1232    0.2263])
plot((log2(p)),(log2(1./f)),'r')
ax1 = gca;
set(gca,'YLim',log2([min(period),max(period)]), ...
    'Xtick',log2([1 10 100 1000 10000 100000]),...
    'XTicklabel',{'10^0','10^1','10^2','10^3','10^4','10^5'},...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',ytl,...
    'ydir','r',...
    'xdir','r')
hold on
% Plot spectral slope fit
plot((log2(0.01*(1./f).^1.53)),(log2(1./f)))

% Plot sinc function
ylabel('Frequency (yrs^{-1})')
xlabel('Power, (\delta^{18}O)^2')
set(gca,'xColor','r')

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on
plot(((sinc(1000*f)-sinc(4000*f)).^2),log2(1./f),'b','Parent',ax2)
plot(log2([min(p),max(p)]),log2([5298,5298]),'color',0.5*[1 1 1],'LineWidth',1)
plot(log2([min(p),max(p)]),log2([2257,2257]),'color',0.5*[1 1 1],'LineWidth',1)
xlim([0,1.5])
xlabel('Power transfer function')

set(gca,'YLim',log2([min(period),max(period)]), ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',{},...
    'ydir','r',...
    'xdir','r',...
    'xcolor','b')


%%
export_fig('-pdf','Figs/NGRIP')
