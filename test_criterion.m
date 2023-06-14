function [True_val,R,R_noise,sigma] = test_criterion

 N  = 1024;
 t  = (0:N-1)/N;
 Nfft = N;
 fs = 0:N/2-1;
 phi2 = 260*t; 
 phi1 = 230*t;
 s = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));
 sigma = 0.024:0.001:0.034;
 % we give an illustration when there are no bubbles
 sigma_test = 0.03;
 [h, Lh]    = create_gaussian_window(N,Nfft,sigma_test);
 %we compute the STFT and the modulation operator
 [STFT, ~,~,Q] = FM_operators(s,Nfft,h, Lh, sigma_test);
 [Cs,ind,jmax,~,~,~] = R_RD_multi(STFT,Lh,Q,0);
 
 figure 
 imagesc(t(Lh:N-Lh),fs,abs(STFT(1:Nfft/2,Lh:N-Lh)));
 hold on;
 xn=0.4;xN=0.5;yn1=210;yN1=280;
 xlim([xn xN]);
 ylim([yn1 yN1]);
 set(gca,'ydir','normal');
 xlabel('time','FontSize',30);
 ylabel('frequency','FontSize',30);
 title('spectrogram','FontSize',30);
 ax = gca;
 set(gca,'TickLength',[0 0])
 set(gca,'Yticklabel',[]) 
 set(gca,'Xticklabel',[])
 for k = 1:jmax
  plot(t(ind{k}(:)+Lh-1),Cs{k}(:)-1,'Linewidth',2);
  legend off;
 end
 hold off;
 %only one equation to determine the existence of TFBs
 True_val = (sqrt(pi/2) *sigma*(30) <=  1);
 
 %we consider a function that tells us whether a intersection of ridges has
 %been found, for us intersection of ridges equal TFBs
 R = zeros(1,length(sigma));
 %R is supposed to be equal to one when a bubble has been detected
 for k = 1:length(sigma)  
  [h, Lh] = create_gaussian_window(N,Nfft,sigma(k));
  %we compute the STFT and the modulation operator
  [STFT, ~,~,Q] = FM_operators(s,Nfft,h, Lh, sigma(k));
  %we detect the ridges and the intersections
  [Cs,~,~,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,0);
  points = calcul_points_bubbles(STFT,Lh,Cs,Tx_ridge,Ap_ridge,Pos_ridge,1);
  A = isempty(points); % no intersection detected 
  R(k)= 1-A;
 end
 %we perform the same computation but with some noise superimposed. 
 SNR = [20 10 5]; 
 nb_real = 30;
 R_noise_int = zeros(nb_real,length(SNR),length(sigma));
 R_noise = zeros(length(SNR),length(sigma));
 
 for r = 1:length(sigma)
  r   
  [h, Lh] = create_gaussian_window(N,Nfft,sigma(r));
  for k = 1:length(SNR)  
   k
   for p = 1:nb_real  
     n = randn(N,1)+1i*randn(N,1);
     [sn] = sigmerge(s(:),n,SNR(k));
     [STFT, ~,~,Q] = FM_operators(sn,Nfft,h, Lh, sigma(r));
     [Cs,~,~,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,1);
     points = calcul_points_bubbles(STFT,Lh,Cs,Tx_ridge,Ap_ridge,Pos_ridge,1);
     A = isempty(points);
     R_noise_int(p,k,r) = 1-A;
   end
  end
 end 
 R_noise(:,:) = mean(R_noise_int);
end
