% Compute the self similarity join of time series
% Author information omitted for ICDM review.
% For details of the AAMP_Optimized algorithm, see:
% "Efficient Matrix Profile Algorithms for Normalized and Non-Normalized Distances", submitted to KDE 2021.
% Usage:
% [mindist, minind] = AAMP_Optimized(X, m)
% Output:
%     mindist: matrix porfile of the self-join (vector)
%     minind: matrix porfile index of the self-join (vector)

% Input:
%     X: input time series (vector)
%     m: interested subsequence length (scalar)
%%
function [mindist, minind] = AAMP(X,m)
    exc_zone = round(m / 2);
    [~, Nb]=size(X); 
    s = Nb-m;
    Dmin = realmax*ones(1,s+1);
    minind = ones(1,s+1);
    
    matchFlag = false;
    for k = 1:1:s
        
        D = sum((X(1:m) - X(k+1:m+k)).^2); 
        
  %      exc_ed_first = min(pro_len, 1+exc_zone);
        if(k > exc_zone)
            matchFlag = true;
        end
        if ( (D < Dmin(1)) && (matchFlag)  )
            Dmin(1) = D;
            minind(1) = k+1;
        end
        
        if ( (D < Dmin(k+1)) && (matchFlag)  )
            Dmin(k+1) = D;
            minind(k+1) = 1;
        end
        
        for i = 2:1:s-k+1
            kplusi = k+i;
            D = D - (X(i-1) - X(kplusi-1))^2 +(X(i+m-1)-X(kplusi+m-1))^2;   % X(i-1) is the fist element of previous query sub-sequence
                                                                            % X(kplusi-1) is the fist element of previous target sub-sequence
                                                                            % X(i+m-1) is the last element of new query sub-sequence
                                                                            % X(kplusi+m-1) is the last element of new target sub-sequence
            if ((Dmin(i) >  D) && (matchFlag))
                minind(i) = kplusi;
                Dmin(i) = D;
            end

            if ((Dmin(kplusi) >  D) && (matchFlag))
                minind(kplusi) = i;
                Dmin(kplusi) = D;
            end
        end
    end
    mindist = sqrt(Dmin);
end