% To show the problem of sigma = 0 in original STOMP algorithm and how the
% proposed algorithm "AAMP" is capable to overcome this issue

% To run this code, you need to download the dataset from link : 
% https://drive.google.com/drive/folders/10WHOK5qEaUiZREW5Tf9dReRJoAua4dMy

% Author information omitted for KDE review.
% For details of the code, see:
% "Efficient Matrix Profile Algorithms for Normalized and Non-Normalized Distances", submitted to KDE 2021.

clear
close all
clc


targetFilePath = 'Seismic Dataset/0194.dat';

getTable = readtable(targetFilePath);

noOfOutliers = 10;
winLen = 15;

Varib2 = getTable.Var2;
Varib3 = getTable.Var3;
Varib4 = getTable.Var4;


[~,fileNamOnly,~] = fileparts(targetFilePath);



% The first column
[pro_mul_1_STOMP, pro_idx_1_STOMP] = STOMP(Varib2, winLen);
[pro_mul_1_AAMP, pro_idx_1_AAMP] = AAMP(Varib2', winLen);
[pro_mul_1_ACAMP, pro_idx_1_ACAMP] = ACAMP_Optimized(Varib2, winLen); % ACAMP algo

pro_mul_1_AAMP(isinf(pro_mul_1_AAMP(:,1)), 1) = 0; % if there are some inf values then replace them by zeros so that they couldn't occur in outliers 
[sortVal1, sortIndx1] = sort(pro_mul_1_AAMP, 'descend');



% The second column
[pro_mul_2_AAMP_STOMP, pro_idx_2_AAMP_STOMP] = STOMP(Varib3, winLen);
[pro_mul_2_AAMP, pro_idx_2_AAMP] = AAMP(Varib3', winLen);
[pro_mul_2_ACAMP, pro_idx_2_ACAMP] = ACAMP_Optimized(Varib3, winLen); % ACAMP algo

pro_mul_2_AAMP(isinf(pro_mul_2_AAMP(:,1)), 1) = 0; % if there are some inf values then replace them by zeros so that they couldn't occur in outliers 
[sortVal2, sortIndx2] = sort(pro_mul_2_AAMP, 'descend');



% The third column
[pro_mul_3_AAMP_STOMP, pro_idx_3_AAMP_STOMP] = STOMP(Varib4, winLen);
[pro_mul_3_AAMP, pro_idx_3_AAMP] = AAMP(Varib4', winLen);
[pro_mul_3_ACAMP, pro_idx_3_ACAMP] = ACAMP_Optimized(Varib4, winLen); % ACAMP algo

pro_mul_3_AAMP(isinf(pro_mul_3_AAMP(:,1)), 1) = 0; % if there are some inf values then replace them by zeros so that they couldn't occur in outliers 
[sortVal3, sortIndx3] = sort(pro_mul_3_AAMP, 'descend');



plotTheGraph(Varib2, sortIndx1, noOfOutliers, pro_mul_1_AAMP, pro_mul_1_STOMP, pro_mul_1_ACAMP, winLen, 'Longitude');
% plotTheGraph(Varib3, sortIndx2, noOfOutliers, pro_mul_2_AAMP, winLen, 'Latitude');
plotTheGraph(Varib4, sortIndx3, noOfOutliers, pro_mul_3_AAMP, pro_mul_3_AAMP_STOMP, pro_mul_3_ACAMP,  winLen, 'Height');




function plotTheGraph(Varib2, sortIndx1, noOfOutliers, pro_mul_1_AAMP, pro_mul_1_STOMP, pro_mul_1_ACAMP, winLen, str)

hFig = figure();
subplot(2,2,1);
    
plot(1:length(Varib2), Varib2, 'k-', 'LineWidth',1);
hold on;
for iOut = 1:1:noOfOutliers
    outlierIndx = (sortIndx1(iOut)); % here we get the paticular sub-sequence and then we repaint it accordingly
    plot(outlierIndx:(outlierIndx+winLen-1), Varib2(outlierIndx:(outlierIndx+winLen-1)),'Color','r', 'LineWidth',2);  
end
hold off;

strTitle = strcat('Original Data-', str);
title(strTitle);


subplot(2,2,2);
plot(pro_mul_1_AAMP, 'r-', 'LineWidth',1);
strTitle = strcat('Matrix Profile AAMP-', str);
title(strTitle);

subplot(2,2,3);
plot(pro_mul_1_STOMP, 'b-', 'LineWidth',1);
strTitle = strcat('Matrix Profile STOMP-', str);
title(strTitle);

subplot(2,2,4);
plot(pro_mul_1_ACAMP, 'b-', 'LineWidth',1);
strTitle = strcat('Matrix Profile ACAMP-', str);
title(strTitle);

end