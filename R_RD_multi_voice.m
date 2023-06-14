function [Cs,ind,jmax,Tx_ridge,Ap_ridge,Pos_ridge] = R_RD_multi(STFT,Lh,Q,cas)
 % Function R_RD_multi : extracts the ridges of a multicomponent signal
 %                       using R_RD technique 
 % Inputs:
 %   STFT       : STFT of s
 %   Q          : modulation operator
 %   cas        : 0 if no noise, 1 if noise
  
 % Outputs:
 %   Cs  :  A structure that contains the frequencies of the computed ridges. 
 %   ind :  The time indices associated with this construction
 %   jmax : The number of ridges computed
 %   Tx_ridge : indicates the number of ridges a point belongs to
 %   Ap_ridge : indicates the ridge number the points belongs to
 %   Pos_ridge : indicates the position on the associated ridge
 %   local_min : determines the local minima
 
 Tx       = abs(STFT); 
 [Nfft,N] = size(Tx);

 nb_max = 20; % we consider a maximum of 20 ridge portions
 Cs  = cell(1,nb_max); %the ridges 
 ind = cell(1,nb_max); %the corresponding indices
 %we compute the set of the local maxima
 
 A = zeros(Nfft,N);
 for k = 1:N
  A(:,k)   = islocalmax(Tx(:,k));
 end
 
 if cas == 0
  Tx = Tx.*A; %we consider the values at local maxima
  thresh = 0; % no threshold put on the coefficients
 else
  % where there is noise we only consider the extrema maxima that above the noise level 
  Y2 = real(STFT);
  gamma_estime = median(abs(Y2(:)))/0.6745;
  thresh = 3*gamma_estime;
  Tx = Tx.*A; %we consider the values at local maxima
 end

 %To initialize we consider the whole set of local maxima
 %we restrict ourselves to inner points.
 Tx = Tx(:,Lh:N-Lh);
 Tx_init = Tx; 
 sz = size(Tx);
 Tx_ridge  = zeros(size(Tx));
 Ap_ridge  = zeros(sz(1),sz(2),2);
 Pos_ridge = zeros(sz(1),sz(2),2); 
 
 % we compute the ridges until the length of the ridges is below 10 
 B0 = max(max(Tx));
 B = B0;
 
 j=1;
 while (B > B0/20)&&(j <= nb_max)
  %extraction of the ridges one at a time   
  %indcur corresponds to the time indices associated with the extracted ridge
  [Cs_cur,indcur] = R_RD(Tx_init,Tx,Lh,Q(:,Lh:N-Lh),thresh);
  if (isempty(indcur) == 0)
   if (length(indcur) >= 300)
    for b=1:length(indcur)
     k = Cs_cur(b);
     n = indcur(b);
     Tx_init(k,n) = 0; % this local maximum can no longer be considered for initialization
     Tx_ridge(k,n) = Tx_ridge(k,n)+1;%we remark that (k,n) is a ridge point
    
     if (Tx_ridge(k,n) == 1)
      Ap_ridge(k,n,1)  = j; % (n,k) belongs to the ridge j
      Pos_ridge(k,n,1) = b; % with position b on that ridge
     elseif (Tx_ridge(k,n) == 2)
      Ap_ridge(k,n,2) = j; % (n,k) belongs to the ridge j
      Pos_ridge(k,n,2) = b; % with position b on that ridge
     end    
    end
  
    Cs{j}  = Cs_cur;
    ind{j} = indcur; 
    j = j+1;
    B = max(max(Tx_init));
   else
   % we remove the detected ridge portion to avoid the alogirthl to get trapped   
    for b=1:length(indcur)
     k = Cs_cur(b);
     n = indcur(b);
     Tx_init(k,n) = 0;
    end
    B = max(max(Tx_init));
   end
  else
   break;   
  end
 end
 jmax =j-1;
end