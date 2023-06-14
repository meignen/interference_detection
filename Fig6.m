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
s = exp(1i*2*pi*phi) + exp(1i*2*pi*(400*t+15*t.^2));

% s = s + 0.1*std(s)*(randn(size(s))+1i*randn(size(s)));

sigma = [0.012 0.018 0.024];
entropies = zeros(1,length(sigma));
cantidad_puntos = zeros(1,length(sigma));
i = 1;
for k=1:length(sigma)
 k   
 [h, Lh]  = create_gaussian_window(N,Nfft,sigma(k));
 %we compute the STFT and the modulation operator
 [STFT, omega,omega2,Q] = FM_operators(s,Nfft,h, Lh, sigma(k));
 SST2 = synchro(STFT,omega2,0);
 
 %we compute the ridges
 [Cs,ind,jmax,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,0);
 points = calcul_points_bubbles(STFT,Lh,Cs,Tx_ridge,Ap_ridge,Pos_ridge,0);
 figure
  imagesc(abs(STFT(1:N/2,Lh:N-Lh)))
  set(gca,'ydir','normal');
  text(50,450,['$\sigma$ = ' num2str(sigma(k))],'fontsize',30,'color',[1 1 1],'Interpreter','Latex')
  xlabel('time','FontSize',30);
  ylabel('frequency','FontSize',30);
  title('spectrogram','FontSize',30);
  ax = gca;
  set(gca,'TickLength',[0 0])
  set(gca,'Yticklabel',[]) 
  set(gca,'Xticklabel',[])
  hold on
  for q = 1:jmax
   plot(ind{q}(:),Cs{q}(:),'Linewidth',2)   
  end
  if isempty(points) == 0
   plot(points(:,2),points(:,1),'*','Linewidth',2,'Markersize',10,'Color','r');
  end 
  hold off
end 
