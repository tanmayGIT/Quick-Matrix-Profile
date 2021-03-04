% Compute the self similarity join of time series
% Author information omitted for ICDM review.
% For details of the ACMP_Optimized algorithm, see:
% "Efficient Matrix Profile Algorithms for Normalized and Non-Normalized Distances", submitted to KDE 2021.
% Usage:
% [mindist, minind] = ACAMP_Optimized(data, sub_len)
% Output:
%     mindist: matrix porfile of the self-join (vector)
%     minind: matrix porfile index of the self-join (vector)

% Input:
%     data: input time series (vector)
%     sub_len: interested subsequence length (scalar)
%%

function [mindist, minind] = ACAMP_Optimized(data, sub_len)

exc_zone = round(sub_len / 2);
[Nb, ~]=size(data);
s = Nb-sub_len;
Dmin = realmax*ones(1,s+1);
minind = ones(1,s+1);

matchFlag = false;
pro_len = s+1;
%initialization
data_mu = zeros(pro_len, 1);
data_sig = zeros(pro_len, 1);

[data_mu(:, 1), data_sig(:, 1), data(:,1)] = mass_pre(data(:, 1), Nb, sub_len);

for k = 1:1:s
    if(  (data_sig(k+1, 1) ~= 0))
        [D, product_me] = dist_Z_Norm(data, k+1, 1,sub_len, data_mu, data_sig);
        
        D = abs(D);
        D = real(D);
        
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
            product_me = product_me - data(i-1)*data(kplusi-1) + data(i+sub_len-1)*data(kplusi+sub_len-1);
            D = 2 * (sub_len - ( (product_me - (sub_len*data_mu(i)* ...
                data_mu(kplusi))    )   /    (data_sig(i) * data_sig(kplusi)) ) );
            
            D = abs(D);
            D = real(D);
            
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
end

mindist = sqrt(Dmin);

end





function [data_mu, data_sig, data] = mass_pre(data, data_len, sub_len)

% This particular function is taken from the following code :
% https://github.com/mcyeh/mstamp/blob/master/MATLAB/mstamp.m
% The following two functions are modified from the code provided in the following URL
% http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html


nanIndexesTar = (isnan(data(:) ));
infIndexesTar = (isinf(data(:) ));
rowNonInf = (~isfinite(data));

allBadCell = bitor((bitor(nanIndexesTar,infIndexesTar)),rowNonInf);
tAllIndex = 1:numel(data(:));
data(allBadCell) = interp1(tAllIndex(~allBadCell), data(~allBadCell), tAllIndex(allBadCell));

data_cum = cumsum(data);
data2_cum =  cumsum(data.^2);


data_sum = data_cum(sub_len:data_len) - ...
    [0; data_cum(1:data_len-sub_len)];
data2_sum = data2_cum(sub_len:data_len) - ...
    [0; data2_cum(1:data_len-sub_len)];


data_mu = data_sum./sub_len;

data_sig2 = (data2_sum./sub_len)-(data_mu.^2);
data_sig2 = real(data_sig2);
data_sig2 = max(data_sig2, 0);

data_sig = sqrt(data_sig2);

end

function [dist_pro, product_me] = dist_Z_Norm(data, targetIndx, queryIndx,sub_len, ...
    data_mu, data_sig)


query = data(queryIndx:1:((queryIndx+sub_len)-1));
target = data(targetIndx:1:((targetIndx+sub_len)-1));

% compute the product and sum
product_me = sum(target.*query);


% compute the distance profile
dist_pro = 2 * (sub_len - ...
    ( (product_me - (sub_len*data_mu(queryIndx)*data_mu(targetIndx)))    /  (data_sig(queryIndx) * data_sig(targetIndx))  )     );
end


function [dist_pro, last_prod] = mass(data_freq, query, ...
    data_len, sub_len, data_mu, data_sig, query_mu, query_sig)

% This particular function is taken from the following code :
% https://github.com/mcyeh/mstamp/blob/master/MATLAB/mstamp.m


% pre-process query for fft
query = query(end:-1:1);
query(sub_len+1:(sub_len+data_len),1) = 0;

% compute the product
query_freq = fft(query);
product_freq = data_freq.*query_freq;
product = ifft(product_freq);

% compute the distance profile
dist_pro = 2 * (sub_len - ...
    (product(sub_len:data_len) - sub_len*data_mu*query_mu)./...
    (data_sig * query_sig));
last_prod = real(product(sub_len:data_len));
end