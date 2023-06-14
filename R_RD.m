function [c,ind] = R_RD(Tx_init,Tx,Lh,Q,thresh)
 % R_RD : extracts  the ridge curve by maximizing some energy 
 %   The result is an approximation computd by a greedy algorithm.
 %   The algorithm uses several random initializations and then
 %   forward/backward search.
 %
 % INPUTS:    
 %   Tx_init  : local maxima of STFT modulus used for initialization 
 %   Tx       : local maxima of STFT modulus 
 %   Q        : the modulation operator    

 % OUTPUTS:  
 %   c   : vector containing the indices of the ridge. 
 %   ind : indice corresponding to the detected ridge

 [~,N]=size(Tx_init);
 %we find the global maximum of the spectrogram 
 val_max = max(max(Tx_init));
 if (val_max > thresh) 
  [index_row,index_column] = find(Tx_init == val_max);
  %we implement the bilateral relation
  %we consider the first maximum in case there are several
  [ridge] = novel_partial_RD(Tx,index_column(1),index_row(1),Q,thresh);
  ind = ridge{1};
  c   = ridge{2};
 else
  ind =[];
  c   =[];
 end
end
  
function [ridge] = novel_partial_RD(Tx,m,km,Q,thresh)
  
 %we contruct the ridge starting from [m,km]
 ridge = cell(1,2); %we put the abscissae and the ordinate
 rq_n0 = real(Q(km,m));
 abscisse_ridge = m;
 ordinate_ridge = km;
 
 [~,N]=size(Tx);
 %by construction the selected point should not belong to an existing ridge 
 
 % dir = 1 : forward iteration
 k_cur = km;
 rq = rq_n0;
 dir = 1;
 for n=(m + dir):dir:N   
  [k_next] = ridge_step(Tx,rq,k_cur,n,dir); 
  if (k_next ~= 0)
   if (Tx(k_next,n) > thresh) %we only consider prolongating the ridge if it is above the noise threshold   
    rq = real(Q(k_next,n));
    [k_next_1] = ridge_step(Tx,rq,k_next,n-1,-dir); % we check the stability by moving backward   
    if (k_next_1 == k_cur)
     k_cur = k_next;
     abscisse_ridge = [abscisse_ridge n];
     ordinate_ridge = [ordinate_ridge k_next];
    else 
     k_cur = k_next;
     abscisse_ridge = [abscisse_ridge n];
     ordinate_ridge = [ordinate_ridge k_next];
     break;
    end
   else
    break;   
   end
  else 
   break;   
  end
 end
 
 % dir = -1 : backward iteration
   
 k_cur = km;
 rq = rq_n0;
 dir = -1;
   
 for n = (m + dir):dir:1 
  [k_next] = ridge_step(Tx,rq,k_cur,n,dir); 
  if (k_next ~= 0)
   if (Tx(k_next,n) > thresh) %we only consider prolongating the ridge if it is above the noise threshold   
    rq=real(Q(k_next,n));
    [k_next_1] = ridge_step(Tx,rq,k_next,n+1,-dir); % we check the stability by moving backward   
    if (k_next_1 == k_cur) 
     k_cur = k_next;
     abscisse_ridge = [n abscisse_ridge];
     ordinate_ridge = [k_next ordinate_ridge];
    else 
     k_cur = k_next;
     abscisse_ridge = [n abscisse_ridge];
     ordinate_ridge = [k_next ordinate_ridge];
     break;
    end
   else
    break;   
   end
  else
   break;   
  end
 end
 
 ridge{1} = abscisse_ridge;
 ridge{2} = ordinate_ridge;
end

function [k_next_1] = ridge_step(Tx,rq,k_next,n,dir)
  
 %Tx contains the set of local maxima in which to look for the next point
 %k_next_1 : non zero and positive if a value is found, 0 if not
 
 [Nfft, N] = size(Tx);
 rq_grid = round(rq*Nfft/(N^2));
 k_rq    = max(1, min(Nfft, k_next+dir*rq_grid));
 
 %we compute the closest local maxima in the direction given by k_rq 
 V = 0;
 while 1
  a = max(1,k_rq - V);
  b = min(Nfft, k_rq + V);
  if Tx(b, n) > 0
   k_next_1 = b;
   if Tx(a, n) > 0
    if Tx(b, n) < Tx(a, n)
     k_next_1 = a;
    end
   end
   break;
  elseif Tx(a, n) > 0
   k_next_1 = a;
   break;
  elseif a == 1 && b == Nfft
   % no local maximum in at time index n
   k_next_1 = 0;
   break;
  end
  V = V + 1;
 end
end 