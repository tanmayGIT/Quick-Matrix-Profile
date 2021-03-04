
% This code is taken from the following link :
% https://sites.google.com/site/scrimpplusplus/

% You can get the original version of the code by clicking on "GUI tools"
% in this line 

% "Interactive Tool: To demonstrate the interactivity of SCRIMP++, we have 
% created a GUI tool as a courtesy to the community"

% We have just removed the portions, related to GUI and the remaining
% things are same as the original code





%%
% The prototype for interactive matrix profile calculation using SCRIMP++
%
% Yan Zhu, Chin-Chia Michael Yeh and Eamonn Keogh
%
% For details of the algorithm, please refer to our paper:
% Yan Zhu, Chin-Chia Michael Yeh, Zachary Zimmerman, Kaveh Kamgar and
% Eamonn Keogh, Matrix Proivle XI: SCRIMP++: Time Series Motif Discovery
% at Interactive Speeds, ICDM 2018.
%
% Usage:
% [matrixProfile, profileIndex] = SCRIMP_Codde_Test(data, subLen, whichAnyTime);
% Output:
%     matrixProfile: matrix porfile of the self-join (vector)
%     profileIndex: matrix porfile index of the self-join (vector)
%     motifIndex: index of the first, second, and third motifs and their 
%                 associated neighbors when stopped (3x2 cell)
%                 +-----------------------+------------------------+
%                 | indices for 1st motif | neighbors of 1st motif |
%                 +-----------------------+------------------------+
%                 | indices for 2nd motif | neighbors of 2nd motif |
%                 +-----------------------+------------------------+
%                 | indices for 3rd motif | neighbors of 3rd motif |
%                 +-----------------------+------------------------+
%     discordIndex: index of discords when stopped (vector)
% Input:
%     data: input time series (vector)
%     subLen: subsequence length (scalar)
%     whichAnyTime : 2 if you want to run SCRIMP
%                  : 3 if you want to run SCRIMP++ 
%%


function [matrixProfile, profileIndex] = ...
    SCRIMP_Codde_Test(data, subLen, whichAnyTime)
%% options for the algorithm
excZoneLen = round(subLen );
radius = 2.0;
updatePeriod = 1; % in second
anytimeMode = whichAnyTime; % 1: original with mass O(n^2 log n); 
                 % 2: SCRIMP
                 % 3: SCRIMP++
preStep = floor(subLen * 0.25); % s in the paper, which means sampling rate

%% check input
dataLen = length(data);
if subLen > dataLen / 2
    error(['Error: Time series is too short ', ...
        'relative to desired subsequence length']);
end
if subLen < 4
    error('Error: Subsequence length must be at least 4');
end
if subLen > dataLen / 5
    error('Error: subsequenceLength > dataLength / 20')
end
if dataLen == size(data, 2)
    data = data';
end


%% locate nan and inf
proLen = dataLen - subLen + 1;
isSkip = false(proLen, 1);
for i = 1:proLen
    if any(isnan(data(i:i + subLen - 1))) || ...
            any(isinf(data(i:i + subLen - 1)))
        isSkip(i) = true;
    end
end
data(isnan(data) | isinf(data)) = 0;

%% preprocess for matrix profile
[dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen);  % ######  line 2 in of the pseudo code in paper
matrixProfile = inf(proLen, 1);
profileIndex = zeros(proLen, 1);
if anytimeMode == 1
    idxOrder = randperm(proLen);
elseif anytimeMode == 2 
    idxOrder = excZoneLen + 1:proLen;
    idxOrder = idxOrder(randperm(length(idxOrder)));
elseif anytimeMode == 3
    Preidx = 1 : preStep : proLen;
    PreidxOrder = Preidx(randperm(length(Preidx))); % line 3 in of the pseudo code in paper
    
    idxOrder = excZoneLen + 1:proLen; % this is to note all the first element of each subsequence, starting from the 2nd sub-sequence
    idxOrder = idxOrder(randperm(length(idxOrder))); % doing the random combination of these indexes 
end

if anytimeMode ==1 || anytimeMode ==2
    iterationEnd = length(idxOrder);
elseif anytimeMode==3
    iterationEnd = length(PreidxOrder) + length(idxOrder);  % getting the summation or total of number of these indexes 
    dotproduct = zeros(proLen,1);
    refine_distance = inf(proLen,1);
    orig_index = 1:proLen;
end

for i = 1:iterationEnd
    if anytimeMode==1 || anytimeMode==2
        idx = idxOrder(i);
    elseif anytimeMode==3
        if i<= length(PreidxOrder)
            curStage=1;
            idx = PreidxOrder(i);
        else
            curStage=2;
            idx = idxOrder(i-length(PreidxOrder));
        end
    end
        
    if isSkip(idx)
        continue
    end
    
    % compute the distance profile
    query = data(idx:idx+subLen-1);
    if anytimeMode == 1
        distProfile = mass(dataFreq, query, dataLen, subLen, ...
            dataMu, dataSig, dataMu(idx), dataSig(idx));
        distProfile = abs(distProfile);
        distProfile = sqrt(distProfile);
        
        distProfile(isSkip) = inf;
        excZoneStart = max(1, idx - excZoneLen);
        excZoneEnd = min(proLen, idx + excZoneLen);
        distProfile(excZoneStart:excZoneEnd) = inf;
        
        updatePos = distProfile < matrixProfile;
        profileIndex(updatePos) = idx;
        matrixProfile(updatePos) = distProfile(updatePos);
        [matrixProfile(idx), profileIndex(idx)] = min(distProfile);
    elseif anytimeMode == 2 || (anytimeMode == 3 && curStage == 2)
        distProfile = diagonalDist(...
            data, idx, dataLen, subLen, proLen, dataMu, dataSig);
        distProfile = abs(distProfile);
        distProfile = sqrt(distProfile);
        
        pos1 = idx:proLen;
        pos2 = 1:proLen - idx + 1;
        
        updatePos = matrixProfile(pos1) > distProfile;
        profileIndex(pos1(updatePos)) = pos2(updatePos);
        matrixProfile(pos1(updatePos)) = distProfile(updatePos);
        
        
        updatePos = matrixProfile(pos2) > distProfile;
        profileIndex(pos2(updatePos)) = pos1(updatePos);
        matrixProfile(pos2(updatePos)) = distProfile(updatePos);
        
        matrixProfile(isSkip) = inf;
        profileIndex(isSkip) = 0;
        
    elseif (anytimeMode==3 && curStage ==1) %PreSCRIMP
        % sampling
        distProfile = mass(dataFreq, query, dataLen, subLen, ...
            dataMu, dataSig, dataMu(idx), dataSig(idx));
        distProfile = abs(distProfile);
        distProfile = sqrt(distProfile);
        
        distProfile(isSkip) = inf;
        excZoneStart = max(1, idx - excZoneLen);
        excZoneEnd = min(proLen, idx + excZoneLen);
        distProfile(excZoneStart:excZoneEnd) = inf;
        
        updatePos = distProfile < matrixProfile;
        profileIndex(updatePos) = idx;
        matrixProfile(updatePos) = distProfile(updatePos);
        [matrixProfile(idx), profileIndex(idx)] = min(distProfile);
        
        % refinement
        idx_nn = profileIndex(idx);
        idx_diff = idx_nn-idx;
        dotproduct(idx) = (subLen - matrixProfile(idx)^2/2)  *  dataSig(idx)  *  dataSig(idx_nn) + subLen   *    dataMu(idx)   *   dataMu(idx_nn);
    
        endidx = min([proLen, idx+preStep-1, proLen-idx_diff]);
        
        
        dotproduct(idx+1:endidx) =  dotproduct(idx)    +     cumsum(data(idx+subLen:endidx+subLen-1).*data(idx_nn+subLen:endidx+subLen-1+idx_diff)     -      data(idx:endidx-1).*data(idx_nn:endidx-1+idx_diff));
        
        
        refine_distance(idx+1:endidx) = sqrt(abs(2*(subLen-(dotproduct(idx+1:endidx)-subLen*dataMu(idx+1:endidx).*dataMu(idx_nn+1:endidx+idx_diff))./ ...
                                                                                    (dataSig(idx+1:endidx).*dataSig(idx_nn+1:endidx+idx_diff)))));
    
        beginidx = max([1, idx-preStep+1, 1-idx_diff]);
        dotproduct(idx-1:-1:beginidx) = dotproduct(idx)  +  cumsum(data(idx-1:-1:beginidx).*data(idx_nn-1:-1:beginidx+idx_diff)  -  data(idx-1+subLen:-1:beginidx+subLen).*data(idx_nn-1+subLen:-1:beginidx+idx_diff+subLen));
        
        
        refine_distance(beginidx:idx-1) = sqrt(abs(2*(subLen-(dotproduct(beginidx:idx-1)-subLen*dataMu(beginidx:idx-1).*dataMu(beginidx+idx_diff:idx_nn-1)) ...
                                                                                    ./(dataSig(beginidx:idx-1).*dataSig(beginidx+idx_diff:idx_nn-1)))));
    
                                                                                
        update_pos1 = find(refine_distance(beginidx:endidx) < matrixProfile(beginidx:endidx)); 
        matrixProfile(update_pos1+beginidx-1) = refine_distance(update_pos1+beginidx-1);
        profileIndex(update_pos1+beginidx-1) = orig_index(update_pos1+beginidx-1)+idx_diff;
    
        
        update_pos2 = find(refine_distance(beginidx:endidx) < matrixProfile(beginidx+idx_diff:endidx+idx_diff));
        matrixProfile(update_pos2+beginidx+idx_diff-1) = refine_distance(update_pos2+beginidx-1);
        profileIndex(update_pos2+beginidx+idx_diff-1) = orig_index(update_pos2+beginidx+idx_diff-1)-idx_diff;
    end
    
end
end


% The following two functions are modified from the code provided in the 
% following URL
% http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [dataFreq, dataMu, dataSig] = massPre(data, dataLen, subLen)
data(dataLen + 1:(subLen + dataLen)) = 0;
dataFreq = fft(data);
dataCumsum = cumsum(data);
data2Cumsum =  cumsum(data .^ 2);
data2Sum = data2Cumsum(subLen:dataLen) - ...
    [0; data2Cumsum(1:dataLen - subLen)];
dataSum = dataCumsum(subLen:dataLen) - ...
    [0; dataCumsum(1:dataLen - subLen)];
dataMu = dataSum ./ subLen;
data2Sig = (data2Sum ./ subLen) - (dataMu .^ 2);
dataSig = sqrt(data2Sig);
end


function distProfile = mass(dataFreq, query, ...
    dataLen, subLen, dataMu, dataSig, queryMu, querySig)
query = query(end:-1:1);
query(subLen+1:(subLen+dataLen)) = 0;
queryFreq = fft(query);
productFreq = dataFreq .* queryFreq;
product = ifft(productFreq);
distProfile = 2 * (subLen - ...
    (product(subLen:dataLen) - subLen * dataMu * queryMu) ./ ...
    (dataSig * querySig));

end

function distProfile = diagonalDist(...
    data, idx, dataLen, subLen, proLen, dataMu, dataSig)
xTerm = ones(proLen - idx + 1, 1) * ...
    (data(idx:idx + subLen - 1)' * data(1:subLen)); % idx is holding the random sub-sequence index, we cut this particular sub-sequence
                                                    % and then get the dot product with very 1st sub-sequence. Then as we need to calculate
                                                    % the distance with the remaining subsequent sub-sequences, so we multiply it a vector
                                                    % of 1's of size (proLen - idx + 1) which means the starting from idx index until the last sub-sequence index 
mTerm = data(idx:proLen - 1) .* ...    % this is to get t_(j-1) for all the sub-sequences after the index idx  : remember this vector should contain the first element of the previous sub-sequence
    data(1:proLen - idx);              % this is to get t_(i-1) for all the sub-sequences from 1 to (proLen - idx)  : remember this vector should contain the first element of the previous sub-sequence


aTerm = data(idx + subLen:end) .* ... % this is the last element of all the sub-sequences, starting from (idx+1)th sub-sequence
    data(subLen + 1:dataLen - idx + 1); % this is the last element of all the sub-sequences, starting from 2nd sub-sequence untill (dataLen - idx + 1)th sub-sequence
if proLen ~= idx
    xTerm(2:end) = xTerm(2:end) - cumsum(mTerm) + cumsum(aTerm);
end

distProfile = (xTerm - ...
    subLen .* dataMu(idx:end) .* dataMu(1:proLen - idx + 1)) ./ ...
    (subLen .* dataSig(idx:end) .* dataSig(1:proLen - idx + 1));
distProfile = 2 * subLen * (1 - distProfile);
end

function x = zeroOneNorm(x)
x = x - min(x(~isinf(x) & ~isnan(x)));
x = x / max(x(~isinf(x) & ~isnan(x)));
end
