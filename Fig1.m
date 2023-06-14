clc;clear all;close all;
%% signal
N  = 1024;
t  = (0:N-1)/N;
fs = 0:N-1;
ft = 0:N-1;
bt = 1:N;
phi1 = 260*t+100*t.^2;phi2 = 220*t+100*t.^2;
s = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));
Nfft = N; 

sigma = 0.03;
[h, Lh] = create_gaussian_window(N, Nfft, sigma);
[STFT, ~,~,Q] = FM_operators(s,Nfft,h, Lh, sigma);
[Cs,ind,~,~,~,~] = R_RD_multi(STFT,Lh,Q,0);
xn=0.4;xN=0.5;yn1=290;yN1=390;

figure;
imagesc(t,fs,abs(STFT).^2);
set(gca,'ydir','normal');
xlabel('time','FontSize',30);
ylabel('frequency','FontSize',30);
title('spectrogram','FontSize',30);
ax = gca;
set(gca,'TickLength',[0 0])
set(gca,'Yticklabel',[]) 
set(gca,'Xticklabel',[])
%ax.FontSize = 20;
xlim([xn xN]);
ylim([yn1 yN1]);
hold on;
plot(t(ind{1}(:)+Lh-1),Cs{1}(:)-1,'Linewidth',2);
plot(t(ind{2}(:)+Lh-1),Cs{2}(:)-1,'Linewidth',2);

[maxi, mini] = get_points_min_max(abs(STFT));

[n,m] = size(STFT);
[mi1, mi2] = ind2sub([n, m],find(mini));
[ma1, ma2] = ind2sub([n, m],find(maxi));
plot(t(mi2),ft(mi1),'bo','Linewidth',2,'MarkerSize',10);
plot(t(ma2),ft(ma1),'ro','Linewidth',2,'MarkerSize',10);
legend('upper ridge','lower ridge','zeros','local maxima','FontSize',25);

phi1 = 260*t+100*t.^2;phi2 = 240*t+100*t.^2;
s = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));
[STFT, ~,~,Q] = FM_operators(s,Nfft,h, Lh, sigma);
[Cs,ind,jmax,~,~,~] = R_RD_multi(STFT,Lh,Q,0);

figure;
imagesc(t,fs,abs(STFT).^2);
set(gca,'ydir','normal');
xlabel('time','FontSize',30);
ylabel('frequency','FontSize',30);
title('spectrogram','FontSize',30);
ax = gca;
set(gca,'TickLength',[0 0])
set(gca,'Yticklabel',[]) 
set(gca,'Xticklabel',[])
%ax.FontSize = 20;
xlim([xn xN]);
ylim([yn1 yN1]);
hold on;
[maxi, mini] = get_points_min_max(abs(STFT));
[n,m] = size(STFT);
[mi1, mi2] = ind2sub([n, m],find(mini));
[ma1, ma2] = ind2sub([n, m],find(maxi));
plot(t(mi2),ft(mi1),'bo','Linewidth',2,'MarkerSize',10);
plot(t(ma2),ft(ma1),'ro','Linewidth',2,'MarkerSize',10);
for k = 1:jmax
 plot(t(ind{k}(:)+Lh-1),Cs{k}(:)-1,'Linewidth',2);
 legend off;
end
legend('zeros','local maxima','FontSize',25);
