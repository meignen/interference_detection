clear; close all

%definition of the signal
N = 1024;
t  = (0:N-1)/N;
fs = 0:N/2-1;
ft = 1:N;
Nfft = N;

%phip = 256 + 100*cos(2*pi*(2*t+8*t.^2));%chirp rate = 6
phip = 256 + 100*cos(2*pi*(2*t+4*t.^2));%chirp rate = 6
phip_prim = -100*(2*pi*(2+ 8*t)).*sin(2*pi*(2*t+4*t.^2));
phip_prim = 356+ cumsum(phip_prim)/N;
phi = cumsum(phip)/N;
s = exp(1i*2*pi*phi)+exp(1i*2*pi*(400*t+15*t.^2));
sigma = [0.012 0.014 0.016 0.018 0.020 0.022 0.024];

for k=1:length(sigma)   
 [h, Lh]       = create_gaussian_window(N,Nfft,sigma(k));
 %we compute the STFT and the modulation operator
 [STFT, omega,omega2,Q] = FM_operators(s,Nfft,h, Lh, sigma(k));
 
 [Cs,~,~,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,0);
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
xlim([sigma(1) sigma(end)]);
ylim([17 20]);
ylabel('Rényi entropy')
xlabel('$\sigma$','Interpreter','Latex')
legend('Noise free','SNR in = 20 dB','SNR in = 10 dB','SNR in = 0 dB')
set(gca,'fontsize',30)
figure
plot(sigma,[cantidad_puntos_no_noise;squeeze(mean(cantidad_puntos,2))],'o-','Linewidth',2,'MarkerSize',10)
xlim([sigma(1) sigma(end)]);
ylabel('number of TFB points')
xlabel('$\sigma$','Interpreter','Latex')
legend('Noise free','SNR in = 20 dB','SNR in = 10 dB','SNR in = 0 dB')
set(gca,'fontsize',30)
