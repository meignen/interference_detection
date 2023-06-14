clc;clear all;close all;
%% signal
N  = 1024;
t  = (0:N-1)/N;
fs = 0:N-1;
ft = 0:N-1;
bt = 1:N;

phi2 = 260*t;phi1 = 240*t;
s = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));
phi1_p = 240;
phi2_p = 260;
 
sigma = 0.02;
%computation of the singularities, here we have an analytic expression
alpha = sqrt(pi/2)*sigma*(phi2_p-phi1_p);
eta_star = (phi1_p+phi2_p)/2;
A = acos(-1+2*alpha^2)/(2*pi);
k = 0;
cur = 0;
T = [];
while cur < 1
 t_k1 = (k-A)/(phi2_p-phi1_p);
 t_k2 = (k+A)/(phi2_p-phi1_p);
 T = [T t_k1 t_k2];
 cur = t_k1;
 k = k+1;
end 
T = T((T>0)&(T<1));
singularities = zeros(2,length(T));
singularities(1,:) = T;
singularities(2,:) = eta_star*ones(1,length(T));    

Nfft = N;
[h, Lh] = create_gaussian_window(N, Nfft, sigma);
[STFT, ~,~,Q] = FM_operators(s,Nfft,h, Lh, sigma);
[Cs,ind,jmax,~,~,~] = R_RD_multi(STFT,Lh,Q,0);
xn=0.5;xN=0.6;yn1=210;yN1=280;
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

plot(singularities(1,:),singularities(2,:),'*','Linewidth',2,'MarkerSize',10,'Color','r');
for q = 1:jmax
  plot(t(ind{q}(:)+Lh-1),Cs{q}(:)-1,'Linewidth',2)   
  legend off;
end
legend('zeros','local maxima','singularities of $\widehat{q_f}$','Interpreter','latex','FontSize',25);

hold off

phi1 = 240*t;phi2 = 260*t;
s = 1.3*exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));

sigma = 0.03;
alpha = sqrt(pi/2)*sigma*(phi2_p-phi1_p);
eta_star = (phi1_p+phi2_p)/2+log(1.3)/(2*pi*sigma^2*(phi2_p-phi1_p));
A = acos(-1+2*alpha^2)/(2*pi);
k = 0;
cur = 0;
T = [];
while cur < 1
 t_k1 = (k-A)/(phi2_p-phi1_p);
 t_k2 = (k+A)/(phi2_p-phi1_p);
 T = [T t_k1 t_k2];
 cur = t_k1;
 k = k+1;
end 
T = T((T>0)&(T<1));
singularities = zeros(2,length(T));
singularities(1,:) = T;
singularities(2,:) = eta_star*ones(1,length(T));    

Nfft = N;
[h, Lh] = create_gaussian_window(N, Nfft, sigma);
%q_denome  = Calcul_singularity_Q(s, Nfft, h, Lh, sigma);
[STFT, ~,~,Q] = FM_operators(s,Nfft,h, Lh, sigma);
[Cs,ind,jmax,~,~,~] = R_RD_multi(STFT,Lh,Q,0);
xn=0.5;xN=0.6;yn1=210;yN1=280;
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

% [~, mini] = get_points_min_max(abs(q_denome));
% [n,m] = size(q_denome);
% [mi1, mi2] = ind2sub([n, m],find(mini));
% plot(t(mi2),ft(mi1),'*','Linewidth',2,'MarkerSize',10,'Color','r');
plot(singularities(1,:),singularities(2,:),'*','Linewidth',2,'MarkerSize',10,'Color','r');
for q = 1:jmax
  plot(t(ind{q}(:)+Lh-1),Cs{q}(:)-1,'Linewidth',2)   
  legend off;
end
legend('zeros','local maxima','singularities of $\widehat{q_f}$','Interpreter','latex','FontSize',25);

hold off