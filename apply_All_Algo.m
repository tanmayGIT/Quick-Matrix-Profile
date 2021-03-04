function apply_All_Algo(keepAllTargetTogether, subSeqLen)


tic
[AAMP_pro_mul, AAMP_pro_idx] = AAMP(keepAllTargetTogether', subSeqLen);
wholetimeAAMP = toc;
fprintf('The time taken by AAMP : %d \n', wholetimeAAMP);


tic
[STOMP_pro_mul, STOMP_pro_idx] = STOMP(keepAllTargetTogether, subSeqLen);
wholetimeSTOMP = toc;
fprintf('The time taken by STOMP : %d \n', wholetimeSTOMP);

tic
[SCRIMP_pro_mul, SCRIMP_pro_idx] = SCRIMP_Codde_Test(keepAllTargetTogether, subSeqLen, 2); % SCRIMP algo
wholetimeSCRIMP = toc;
fprintf('The time taken by SCRIMP : %d \n', wholetimeSCRIMP);


tic
[SCRIMP_Plus_pro_mul, SCRIMP_Plus_pro_idx] = SCRIMP_Codde_Test(keepAllTargetTogether, subSeqLen, 3); % SCRIMP++ algo
wholetimeSCRIMP_Plus = toc;
fprintf('The time taken by SCRIMP Plus Plus : %d \n', wholetimeSCRIMP_Plus);

tic
[ACAMP_mul, ACAMP_pro_idx] = ACAMP_Optimized(keepAllTargetTogether, subSeqLen); % ACAMP algo
wholetimeACAMP = toc;
fprintf('The time taken by ACAMP : %d \n', wholetimeACAMP);


tic
[ACAMP_Colez_pro_mul, ACAMP_Colez_pro_idx] = ACAMP_1(keepAllTargetTogether, subSeqLen); % ACAMP Colez algo and formula
wholetimeACAMP_Colez = toc;
fprintf('The time taken by ACAMP_Colez : %d \n\n\n', wholetimeACAMP_Colez);



% [maxAAMPDist, maxAAMPDistIndx] = max(AAMP_pro_mul);
% aampMatchIndx = AAMP_pro_idx(maxAAMPDistIndx);
% stompMatchIndx = SCRIMP_Plus_pro_idx(maxAAMPDistIndx);
% plotTheGraph(keepAllTargetTogether, maxAAMPDistIndx, aampMatchIndx, stompMatchIndx, subSeqLen);


return;
end





function plotTheGraph(Varib2, origOutlierIndx, aampMatchIndx, stompMatchIndx, winLen)

hFig1 = figure();
plot(1:length(Varib2), Varib2, 'b-', 'LineWidth',1);
xlim([0 length(Varib2)]);

hold on;

plot(origOutlierIndx:(origOutlierIndx+winLen-1), Varib2(origOutlierIndx:(origOutlierIndx+winLen-1)),'Color','m', 'LineWidth',1);  
plot(aampMatchIndx:(aampMatchIndx+winLen-1), Varib2(aampMatchIndx:(aampMatchIndx+winLen-1)),'Color','r', 'LineWidth',1); 
plot(stompMatchIndx:(stompMatchIndx+winLen-1), Varib2(stompMatchIndx:(stompMatchIndx+winLen-1)),'Color','k', 'LineWidth',1);

hold off;
end




function [keepAllTargetTogether, keepDataFileInfo] = ConcatenateAllSeries(subFoldersTarget )
% Get all the target sequences together and merged
keepAllTargetTogether = zeros(1,1);
fullPtCnt = 1;
keepDataFileInfo = cell(1,1);

getGoodFileCnt = 1;
for lTarget = 1:1:size(subFoldersTarget,1) % get the target files
    C1Target = subFoldersTarget(lTarget,:);
    getLengthTarget = length(C1Target);
    
    C1TargetArr = zeros(getLengthTarget,1);
    
    C1TargetArr(:,1) = C1Target(1, :);
    
    clearvars C1Target
    
    if (getGoodFileCnt == 1)
        keepAllTargetTogether(1:(getLengthTarget),1) = C1TargetArr(1:end);
        
        keepDataFileInfo{getGoodFileCnt,1}.FileNum = lTarget;
        keepDataFileInfo{getGoodFileCnt,1}.DataStart = 1;
        keepDataFileInfo{getGoodFileCnt,1}.DataEnd = (getLengthTarget);
        fullPtCnt = fullPtCnt + (getLengthTarget);
    else
        keepAllTargetBackup = keepAllTargetTogether;
        keepAllTargetTogether = zeros ((size(keepAllTargetBackup,1)+size(C1TargetArr,1)),1);
        keepAllTargetTogether(1:size(keepAllTargetBackup,1),:) = keepAllTargetBackup(:,:);
        
        clearvars keepAllTargetBackup
        keepAllTargetTogether(fullPtCnt:(fullPtCnt+(getLengthTarget-1)),1) = C1TargetArr(:,1);
        
        keepDataFileInfo{getGoodFileCnt,1}.FileNum = lTarget;
        keepDataFileInfo{getGoodFileCnt,1}.DataStart = fullPtCnt;
        keepDataFileInfo{getGoodFileCnt,1}.DataEnd = (fullPtCnt+(getLengthTarget-1));
        
        if(( keepDataFileInfo{getGoodFileCnt,1}.DataEnd - ...
                keepDataFileInfo{getGoodFileCnt,1}.DataStart) > getLengthTarget)
            error('need to check it, there could be some problem here');
        end
        fullPtCnt = fullPtCnt + (getLengthTarget);
    end
    getGoodFileCnt = getGoodFileCnt + 1;
end

return
end

% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [data_freq, data_mu, data_sig,data] = mass_pre(data, data_len, sub_len)

% nanIndexesTar = find(isnan(data(:) ));
% if(~isempty(nanIndexesTar) )
%     data(nanIndexesTar,1) = 0;
% end
% infIndexesTar = find(isinf(data(:) ));
% if(~isempty(infIndexesTar) )
%     data(infIndexesTar,1) = 0;
% end
% [rowNonInf,~,~] = find(~isfinite(data));
% data(rowNonInf) = 0; % making NAN and INF to zero to avaoid the problems

nanIndexesTar = (isnan(data(:) ));
infIndexesTar = (isinf(data(:) ));
rowNonInf = (~isfinite(data));

allBadCell = bitor((bitor(nanIndexesTar,infIndexesTar)),rowNonInf);
tAllIndex = 1:numel(data(:));
data(allBadCell) = interp1(tAllIndex(~allBadCell), data(~allBadCell), tAllIndex(allBadCell));

data(data_len+1:(sub_len+data_len)) = 0;
data_freq = fft(data);
data_cum = cumsum(data);
data2_cum =  cumsum(data.^2);
data2_sum = data2_cum(sub_len:data_len) - ...
    [0; data2_cum(1:data_len-sub_len)];
data_sum = data_cum(sub_len:data_len) - ...
    [0; data_cum(1:data_len-sub_len)];
data_mu = data_sum./sub_len;
data_sig2 = (data2_sum./sub_len)-(data_mu.^2);
data_sig2 = real(data_sig2);
data_sig2 = max(data_sig2, 0);
data_sig = sqrt(data_sig2);
end


function [dist_pro, last_prod] = mass(data_freq, query, ...
    data_len, sub_len, data_mu, data_sig, query_mu, query_sig, lnNum)
% pre-process query for fft
query = query(end:-1:1);
query(sub_len+1:(sub_len+data_len)) = 0;

% compute the product
query_freq = fft(query);
product_freq = data_freq.*query_freq;
product = ifft(product_freq);

% compute the distance profile
dist_pro = 2 * (sub_len - ...
    (product(sub_len:data_len) - sub_len*data_mu*query_mu)./...
    (data_sig * query_sig));
if(~isfinite(dist_pro))
    fprintf('The bad line number %d \n', lnNum);
    fprintf('The data_freq is \n');
    disp(data_freq);
    
    fprintf('The query is \n');
    disp(query);
    
    fprintf('The data_mu %d data_sig %d query_mu %d query_sig %d \n', data_mu, data_sig, query_mu, query_sig);
    disp('need checking');
end
last_prod = real(product(sub_len:data_len));
end

