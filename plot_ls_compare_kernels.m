t = 200;
rtp = 10:10:500;
[out] = ls_diff_sinc_heavyside(t,rtp);

%% Plot
figure(2)
plot(rtp,out,'k-o')
grid
set(gcf,'position',[ 440   608   568   190],'color','w')
ylabel({'RMSE between ideal' 'and sinc^2(\tau_b\nu) transfer functions'})
xlabel('\tau_b')

export_fig('-pdf','Figs/kernel_misfit')
