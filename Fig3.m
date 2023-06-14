clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = './';
%% required paths 
addpath(folder);

%definition of the signal
N = 1024;
X(:,1) = (fmconst(N, 0.15));
X(:,2) = (fmlin(N,0.08,0.25));
X(:,3) = (fmsin(N,0.25,0.40,320,1,0.3,+1));
s = sum(X,2);
X = transpose(X);
t  = (0:N-1)/N;
fs = 0:N/2-1;
ft = 1:N;
Nfft = N;

sigma = [0.03 0.04 0.05];

for k=1:length(sigma)
 %we compute the filter
 [h, Lh]       = create_gaussian_window(N,Nfft,sigma(k));
 %we compute the STFT and the modulation operator
 [STFT, ~,~,Q] = FM_operators(s,Nfft,h, Lh, sigma(k));
 %we compute the ridge until a certain amount of the energy is contained in
 %the vicinity of the ridges
 [Cs,ind,jmax,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,0);
 points = calcul_points_bubbles(STFT,Lh,Cs,Tx_ridge,Ap_ridge,Pos_ridge,0);
 
 figure 
 imagesc(t(Lh:N-Lh),fs,abs(STFT(1:Nfft/2,Lh:N-Lh)));
 set(gca,'ydir','normal');
 xlabel('time','FontSize',30);
 ylabel('frequency','FontSize',30);
 title('spectrogram','FontSize',30);
 ax = gca;
 set(gca,'TickLength',[0 0])
 set(gca,'Yticklabel',[]) 
 set(gca,'Xticklabel',[])
 hold on;

 for q = 1:jmax
  plot(t(ind{q}(:)+Lh-1),Cs{q}(:)-1,'Linewidth',2)   
 end
 plot(t(points(:,2)+Lh-1),points(:,1)-1,'*','Linewidth',2,'Markersize',10,'Color','r');
 hold off;
 
end

