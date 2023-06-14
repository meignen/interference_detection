clear; close all

%definition of the signal
N = 1024;
X(:,1) = (fmconst(N, 0.15));
X(:,2) = (fmlin(N,0.08,0.25));
X(:,3) = (fmsin(N,0.25,0.40,320,1,0.3,+1));
s = sum(X,2);
t  = (0:N-1)/N;
fs = 0:N/2-1;
ft = 1:N;
Nfft = N;


sigma = [0.03 0.035 0.04 0.045 0.05 0.055];

for k=1:length(sigma)   
 [h, Lh]       = create_gaussian_window(N,Nfft,sigma(k));
 %we compute the STFT and the modulation operator
 [STFT, omega,omega2,Q] = FM_operators(s,Nfft,h, Lh, sigma(k));
 
 [Cs,ind,jmax,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,0);
 points = calcul_points_bubbles(STFT,Lh,Cs,Tx_ridge,Ap_ridge,Pos_ridge,0);

 [aux,~]=size(points);
 cantidad_puntos_no_noise(k) = aux;
 entropies_no_noise(k) = renyi_entropy(abs(STFT(1:N/2,:)).^2,2);
end 

SNR = [20 10 0];
J = 30;
entropies = zeros(3,J,length(sigma));
cantidad_puntos = zeros(3,J,length(sigma));

for k=1:length(sigma)
 k   
 [h, Lh] = create_gaussian_window(N,Nfft,sigma(k));
 for j=1:J
  j
  n = randn(N,1)+1i*randn(N,1);
  for kk=1:3
   [sn] = sigmerge(s(:),n,SNR(kk));
   [STFT, ~,~,Q] = FM_operators(sn,Nfft,h, Lh, sigma(k));
   %we compute the ridge until a certain amount of the energy is contained in
   %the vicinity of the ridges 
   [Cs,~,~,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,1);
   points = calcul_points_bubbles(STFT,Lh,Cs,Tx_ridge,Ap_ridge,Pos_ridge,0);
   [aux,~]=size(points);
   cantidad_puntos(kk,j,k) = aux;
   entropies(kk,j,k) = renyi_entropy(abs(STFT(1:N/2,:)).^2,2);
  end
 end
end

figure; 
plot(sigma,[entropies_no_noise;squeeze(mean(entropies,2))],'o-','Linewidth',2,'MarkerSize',10)
ylabel('Rényi entropy')
xlabel('$\sigma$','Interpreter','Latex')
xlim([sigma(1) sigma(end)]);
ylim([16 18]);
legend('Noise free','SNR in = 20 dB','SNR in = 10 dB','SNR in = 0 dB')
set(gca,'fontsize',30)
figure
plot(sigma,[cantidad_puntos_no_noise;squeeze(mean(cantidad_puntos,2))],'o-','Linewidth',2,'MarkerSize',10)
xlim([sigma(1) sigma(end-1)]);
ylim([0 30]);
ylabel('number of TFBs points')
xlabel('$\sigma$','Interpreter','Latex')
legend('Noise free','SNR in = 20 dB','SNR in = 10 dB','SNR in = 0 dB')
set(gca,'fontsize',30)
