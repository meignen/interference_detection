 clear all; close all;
 [True_val,R,R_noise,sigma] = test_criterion;
 plot(sigma,True_val,sigma,R,'-*',sigma,R_noise(1,:),'-d',sigma,R_noise(2,:),'-o',sigma,R_noise(3,:),'-s','Linewidth',2,'MarkerSize',20);
 legend('ground truth','R','R,SNR = 20 dB','R,SNR = 10 dB','R,SNR = 5 dB');
 xlabel('\sigma');
 ylabel('TFB detection');
 set(gca,'fontsize',30)
 hold off;

 