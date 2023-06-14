function points = calcul_points_bubbles(STFT,Lh,Cs,Tx_ridge,Ap_ridge,Pos_ridge,cas) 
 %input:
 % STFT     : STFT of the signal 
 % Cs       : frequency corresponding to the ridges
 % Tx_ridge :  indicates the number of ridge a point belongs to
 % Ap_ridge :  indicates the ridge index associated with this point
 % Pos_ridge: indicates the position on that ridge  
 % cas      : if 1 stops as soon as a bubble is found 
 %output:
 % points: the points associated with different ridges where interference
 %         are detected. 
 
 [~,N] = size(STFT);
 STFT = STFT(:,Lh:N-Lh);%we consider only inner points
 Tx = abs(STFT);
 
 %we compute the points associated with two ridges
 [points(:,1),points(:,2)]=find((Tx_ridge ==2)); %row index first, colum index second
 %We now identify the loops, to remove those that are not associated with
 %bubbles, a bubbles must contain a single local min
 if (isempty(points) == 0)
  a = size(points);
  sz = size(Tx);
  A = zeros(sz);
  B = zeros(sz);
  for k = 1:sz(1)
   A(k,:) = islocalmin(Tx(k,:));
  end
  for k = 1:sz(2)
   B(:,k) = islocalmin(Tx(:,k));
  end
  C = A.*B;
  [points_min(:,1),points_min(:,2)] = find((C == 1));
  b = size(points_min);
  mini = zeros(sz(1),sz(2));
  for k = 1:b(1) 
    k1 = points_min(k,1);
    k2 = points_min(k,2);
    [~,I] = min([Tx(k1,k2) Tx(k1-1,k2-1) Tx(k1+1,k2+1) Tx(k1-1,k2+1) Tx(k1+1,k2-1)]);
    if (I(1) == 1)
     mini(k1,k2) = 1;  
    end
 end  
 points_new = [];

 for k = 1:a(1)   
  for k1= k+1:a(1)
   if Ap_ridge(points(k,1),points(k,2),:) == Ap_ridge(points(k1,1),points(k1,2),:)
    
    A = min(points(k,2), points(k1,2));
    B = max(points(k,2), points(k1,2));
    
 
    p1 = Pos_ridge(points(k,1),points(k,2),1);
    p2 = Pos_ridge(points(k1,1),points(k1,2),1);
    freq_r1  = Cs{Ap_ridge(points(k,1),points(k,2),1)}(min(p1:p2):max(p1:p2));
    
    p3 = Pos_ridge(points(k,1),points(k,2),2);
    p4 = Pos_ridge(points(k1,1),points(k1,2),2);
    freq_r2  = Cs{Ap_ridge(points(k,1),points(k,2),2)}(min(p3:p4):max(p3:p4));
    
    freq =[freq_r1,freq_r2];
    C = min(freq);
    D = max(freq);
    
    if (sum(sum(mini(C:D,A:B)))==1)
     if cas == 1
      points_new = [points(k,1),points(k,2);points(k1,1),points(k1,2)]; %just to say it is non empty;
      break;   
     else
      points_new = [points_new;points(k,1),points(k,2);points(k1,1),points(k1,2)];
     end
    end    
   end
  end
 end
 points = points_new;
 end
end