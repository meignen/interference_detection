clear; close all

x = importdata('62-a_lhl.wav');
fs = x.fs;
x = x.data;
x = resample(x,1,10);
% x = [x; zeros(2*8192-length(x),1)];

Fs = fs/10;
Lx = length(x);

C_Lx = Lx/Fs;

%% filter design

%sigma_w = 0.005; % 0.05 para N = 8192, 0.025 para N = 4096 - revisar bien los valores de sigma

%% Time frequency

Nfft = 4096;
Nx = 1:Lx;
t = (Nx - 1)/Fs;

max_f_norm = Nfft/2;

% snr_in = 10;
noise_r = randn(1, Lx);

sigma_w = [0.003 0.005 0.007 0.009 0.011 0.013 0.015];% original es 0.01 
entropies_no_noise = zeros(1,length(sigma_w));
cantidad_puntos_no_noise = zeros(1,length(sigma_w));


i = 0; 
for k=1:length(sigma_w)
 k
 [h, Lh]            = create_gaussian_window(Lx,Nfft,sigma_w(k));
 %we compute the STFT and the modulation operator
 [STFT, ~,omega2,Q] = FM_operators(x,Nfft,h, Lh, sigma_w(k));
 [Cs,ind,jmax,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi_voice(STFT(1:Nfft/2,:),Lh,Q,0);
 %we compute the points associated with two ridges
 points=[];
 [points(:,1),points(:,2)]=find((Tx_ridge ==2));
 % we also consider as singularities the points corresponding to an
 % interruption of the ridge
 for j = 1:jmax 
  if (ind{j}(1) > 1)
   points = [points ; Cs{j}(1),ind{j}(1)] ; 
  end
  if (ind{j}(end) < Lx-2*Lh+1)
   points = [points ; Cs{j}(end),ind{j}(end)] ; 
  end
 end
 [aux,~]=size(points);
 cantidad_puntos_no_noise(k) = aux;
 entropies_no_noise(k) = renyi_entropy(abs(STFT(1:Nfft/2,Lh:Lx-Lh)).^2,2);  
end 

SNR = [20 10];
nbreal = 10;
entropies = zeros(length(SNR),nbreal,length(sigma_w));
cantidad_puntos = zeros(length(SNR),nbreal,length(sigma_w));

for k=1:length(sigma_w)
 k
 [h, Lh]       = create_gaussian_window(Lx,Nfft,sigma_w(k));
 for q = 1:length(SNR)
  q
  for p = 1:nbreal  
   p
   n = randn(Lx,1)+1i*randn(Lx,1);
   [sn] = sigmerge(x(:),n,SNR(q));
   %we compute the STFT and the modulation operator
   [STFT, ~,omega2,Q] = FM_operators(sn,Nfft,h, Lh, sigma_w(k));
   [Cs,ind,jmax,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT(1:Nfft/2,:),Lh,Q,1);
   %we compute the points associated with two ridges
   points = [];
   [points(:,1),points(:,2)] = find((Tx_ridge ==2));
   % we also consider as singularities the points corresponding to an
   % interruption of the ridge
   for j = 1:jmax 
    if (ind{j}(1) > 1)
     points = [points ; Cs{j}(1),ind{j}(1)] ; 
    end
    if (ind{j}(end) < Lx-2*Lh+1)
     points = [points ; Cs{j}(end),ind{j}(end)] ; 
    end
   end
   [aux,~]=size(points);
   cantidad_puntos(q,p,k) = aux;
   entropies(q,p,k) = renyi_entropy(abs(STFT(1:Nfft/2,Lh:Lx-Lh)).^2,2);  
  end
 end
end

figure;
plot(sigma_w,[entropies_no_noise;squeeze(mean(entropies,2))],'o-','Linewidth',2,'MarkerSize',10)
ylabel('Rényi entropy')
xlabel('$\sigma$','Interpreter','Latex')
legend('Noise free','SNR in = 20 dB','SNR in = 10 dB')
set(gca,'fontsize',30)
figure
plot(sigma_w,[cantidad_puntos_no_noise;squeeze(mean(cantidad_puntos,2))],'o-','Linewidth',2,'MarkerSize',10)
ylabel('number of singularities')
xlabel('$\sigma$','Interpreter','Latex')
legend('Noise free','SNR in = 20 dB','SNR in = 10 dB')
set(gca,'fontsize',30)
