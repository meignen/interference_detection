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

sigma_w = [0.003 0.005 0.007 0.009 0.011 0.013 0.015];

i = 0; 
for k=1:length(sigma_w)
 k
 [h, Lh]       = create_gaussian_window(Lx,Nfft,sigma_w(k));
 x1 = 8330-Lh+1;
 x2 = 8830-Lh+1;
 y1 = 820;
 y2 = 1150;
 %we compute the STFT and the modulation operator
 [STFT, ~,omega2,Q] = FM_operators(x,Nfft,h, Lh, sigma_w(k));
 SST = synchro(STFT,omega2*Nfft/Lx,0);
 
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
 
 if (k==1)||(k == 3)||(k == 6)
  
  subplot(3,3,i+1)
  imagesc(abs(STFT(1:Nfft/2,Lh:Lx-Lh)));
  set(gca,'ydir','normal');
  text(400,1850,['$\sigma =$ ' num2str(sigma_w(k))],'fontsize',20,'color',[1 1 1],'Interpreter','Latex')
  title('spectrogram');
  set(gca,'fontsize',10);
  set(gca,'TickLength',[0 0])
  set(gca,'Yticklabel',[]) 
  set(gca,'Xticklabel',[])
  hold on
  for q = 1:jmax
   plot(ind{q}(:),Cs{q}(:)-1,'Linewidth',2)   
  end
  if isempty(points) == 0
   plot(points(:,2),points(:,1)-1,'*','Linewidth',2,'Markersize',10,'Color','r');
  end
  rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','green','Linewidth',1.5);
  hold off

  subplot(3,3,i+2)
  imagesc(abs(STFT(1:Nfft/2,Lh:Lx-Lh)));
  xlim([x1 x2])
  ylim([y1 y2])
  set(gca,'ydir','normal');
  title('zoomed spectrogram');
  set(gca,'TickLength',[0 0])
  set(gca,'Yticklabel',[]) 
  set(gca,'Xticklabel',[])
  set(gca,'fontsize',10);
  
  subplot(3,3,i+3)
  imagesc(abs(SST(y1:y2,x1:x2)));
  set(gca,'ydir','normal');
  title('zoomed SST2');
  set(gca,'TickLength',[0 0])
  set(gca,'Yticklabel',[]) 
  set(gca,'Xticklabel',[])
  set(gca,'fontsize',10);
  i = i+3;
 end
end 
