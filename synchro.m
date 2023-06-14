function T = synchro(STFT,omega,gamma)

[neta,nb] = size(STFT);
T = zeros(neta,nb);

% %% Reassignment step
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))>gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%SST1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            k = 1+round(omega(eta,b));
            if k>=1 && k<=neta
                T(k,b)  = T(k,b) + STFT(eta,b);
            end

% %%%%%%%%%%%%%%%%%%%%SST2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             k = 1+round(omega2(eta,b));
%              if k>=1 && k<=neta
%                SST2(k,b)  = SST2(k,b) + STFT(eta,b);
%              end
% 
% %%%%%%%%%%%%%%%%%%%%%SST3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              k = 1+floor(omega3(eta,b));
%              if k>=1 && k<=neta
%                SST3(k,b)  = SST3(k,b) + STFT(eta,b);
%              end
% 
% %%%%%%%%%%%%%%%%%%%%%SST4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              k = 1+floor(omega4(eta,b));
%              if k>=1 && k<=neta
%                SST4(k,b)  = SST4(k,b) + STFT(eta,b);
%              end

        end
    end
end

end