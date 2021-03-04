
% We have written this code by using the code from the following link :
% https://github.com/mcyeh/mstamp/blob/master/MATLAB/mstamp.m


% For details of the algorithm, please refer to paper:
% Yan Zhu, Zachary Zimmerman, Nader Shakibay Senobari, 
% Chin-Chia Michael Yeh, Gareth Funning, Abdullah Mueen, Philip Brisk and Eamonn Keogh
% Matrix Profile II: Exploiting a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and Joins

% Usage:
% [pro_mul, pro_idx] = STOMP(data, sub_len)
% Output:
%     pro_mul: matrix porfile of the self-join (vector)
%     pro_idx: matrix porfile index of the self-join (vector)

% Input:
%     data: input time series (vector)
%     sub_len: interested subsequence length (scalar)


function [pro_mul, pro_idx] = STOMP(data, sub_len)
%% get various length
exc_zone = round(sub_len / 2); % I think this is for trivial matches
data_len = size(data, 1);
pro_len = data_len - sub_len + 1;

%% check input
% if sub_len > data_len / 2
%     error(['Error: Time series is too short relative ', ...
%         'to desired subsequence length']);
% end
% if sub_len < 4
%     error('Error: Subsequence length must be at least 4');
% end


%% initialization
data_freq = zeros((sub_len + data_len), 1);
data_mu = zeros(pro_len, 1);
data_sig = zeros(pro_len, 1);
first_prod = zeros(pro_len, 1);

[data_freq(:, 1), data_mu(:, 1), data_sig(:, 1), data(:,1)] = ...
    mass_pre(data(:, 1), data_len, sub_len);

[~, first_prod(:, 1)] = mass(...
    data_freq(:, 1), data(1:sub_len, 1), data_len, ...
    sub_len, data_mu(:, 1), data_sig(:, 1), ...
    data_mu(1, 1), data_sig(1, 1));

%% compute the matrix profile
pro_mul = Inf(pro_len, 1);
pro_mul_refined = Inf(1, 1);

pro_idx = zeros(pro_len, 1);
pro_idx_refined = zeros(1, 2);

dist_pro = zeros(pro_len, 1);
last_prod = zeros(pro_len, 1);
drop_val = zeros(1, 1);

goodCnt = 1;
for i = 1:pro_len
    
    query = data(i:i+sub_len-1, :);
    if(data_sig(i, 1) ~= 0) % otherwise all the denominator will become 0, hence all last_prod elements become -Inf
        if i==1
            [dist_pro(:, 1), last_prod(:, 1)] = mass(data_freq(:, 1), query(:, 1), ...
                data_len, sub_len, data_mu(:, 1), data_sig(:, 1), data_mu(i, 1), data_sig(i, 1));
            
        else
            last_prod(2:data_len - sub_len + 1, :) = last_prod(1:data_len - sub_len, :) - data(1:data_len - sub_len, :) ...
                .* repmat(drop_val, pro_len - 1, 1) + data(sub_len + 1:data_len, :) .* repmat(query(sub_len, :), pro_len - 1, 1); % (sub_len+1):data_len; this starts from the last
            % element of the 2nd sub-sequence upto the last element of the last sub-sequence possible
            
            last_prod(1, :) = first_prod(i, :); % see the previous line and the last_prod (2:pp, :) means calculating value from 2 to pp for all dimention
            % so that's why 1th row the value is fulfilled by
            dist_pro = 2 * (sub_len - (last_prod ...
                - sub_len * data_mu .* repmat(data_mu(i, :), pro_len, 1)) ...
                ./ (data_sig .* repmat(data_sig(i, :), pro_len, 1)));
        end
        dist_pro = abs(dist_pro);
        dist_pro = real(dist_pro);
        drop_val(:) = query(1, :);
        
        
%         dist_pro(isnan(dist_pro(:,1)),1) = intmax('int32');  % if find nan value then replace by high value
         dist_pro(isinf(dist_pro(:,1)), 1) = intmax('int32'); % if find inf value then replace by high value
%         dist_pro(~isreal(dist_pro(:,1)), 1) = intmax('int32'); % if find complex vallue then replace by high value
%         dist_pro(~isfinite(dist_pro(:,1)), 1) = intmax('int32');
         dist_pro((dist_pro(:,1)< 0), 1) = intmax('int32');   % if the value if negative then replace it by a bigh value
        
        
        % apply exclusion zone
        exc_st = max(1, i - exc_zone);
        exc_ed = min(pro_len, i+exc_zone);
        dist_pro(exc_st:exc_ed, :) = inf;
        
        dist_pro(data_sig < eps) = inf;
        
        [min_val, min_idx] = min(dist_pro);
        pro_mul(i, 1) = min_val;
        pro_idx(i, 1) = min_idx;
        
        pro_mul_refined(goodCnt, 1) = min_val; % the minimum dist
        pro_idx_refined(goodCnt, 1) = i; % the real index of the sub-sequence
        pro_idx_refined(goodCnt, 2) = min_idx; % matched with which index 
        goodCnt = goodCnt+1;
    else
        pro_idx(i, 1) = i;
    end
    
end

pro_mul = sqrt(pro_mul);
end



% This particular function is taken from the following code : 
% https://github.com/mcyeh/mstamp/blob/master/MATLAB/mstamp.m

% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [data_freq, data_mu, data_sig, cleanData] = mass_pre(data, data_len, sub_len)

nanIndexesTar = (isnan(data(:) ));
infIndexesTar = (isinf(data(:) ));
rowNonInf = (~isfinite(data));

allBadCell = bitor((bitor(nanIndexesTar,infIndexesTar)),rowNonInf);
tAllIndex = 1:numel(data(:));
data(allBadCell) = interp1(tAllIndex(~allBadCell), data(~allBadCell), tAllIndex(allBadCell));

cleanData = data;

data(data_len+1:(sub_len+data_len)) = 0;
data_freq = fft(data);

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


% This particular function is taken from the following code : 
% https://github.com/mcyeh/mstamp/blob/master/MATLAB/mstamp.m
function [dist_pro, last_prod] = mass(data_freq, query, ...
    data_len, sub_len, data_mu, data_sig, query_mu, query_sig)
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