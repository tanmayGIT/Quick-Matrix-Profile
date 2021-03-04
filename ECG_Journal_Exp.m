% To show the better outlier detection capabilities of AAMP and ACAMP over
% STOMP and SCRIMP++ algorithm

% To run this code, you need to download the dataset from link : 
% https://drive.google.com/drive/folders/10WHOK5qEaUiZREW5Tf9dReRJoAua4dMy

% Author information omitted for KDE review.
% For details of the code, see:
% "Efficient Matrix Profile Algorithms for Normalized and Non-Normalized Distances", submitted to KDE 2021.

clear
close all
clc


targetFilePath_2 = 'ECG_Data/chfdb_chf13_45590.txt';

targetFilePath_7 = 'ECG_Data/xmitdb_x108_0.txt';


getTable_2 = readtable(targetFilePath_2);
getTable_7 = readtable(targetFilePath_7);

noOfOutliers = 10;
winLen = 15;



Varib2_1 = getTable_2.Var1;
Varib2_2 = getTable_2.Var2;
Varib2_3 = getTable_2.Var3;

Varib7_1 = getTable_7.Var1;
Varib7_2 = getTable_7.Var2;
Varib7_3 = getTable_7.Var3;

[AAMP_pro_mul, AAMP_pro_idx] = AAMP(Varib2_2', 50);
[SCRIMP_Plus_pro_mul, SCRIMP_Plus_pro_idx] = SCRIMP_Codde_Test(Varib2_2, 50, 3); 


[AAMP_pro_mul_7_2, AAMP_pro_idx_7_2] = AAMP(Varib7_2', 50);
[SCRIMP_Plus_pro_mul_7_2, SCRIMP_Plus_pro_idx_7_2] = SCRIMP_Codde_Test(Varib7_2, 50, 3);


plotTheGraph(Varib2_2, AAMP_pro_mul, SCRIMP_Plus_pro_mul);
plotTheGraph(Varib7_2, AAMP_pro_mul_7_2, SCRIMP_Plus_pro_mul_7_2);


function plotTheGraph(Varib2, pro_mul_1_AAMP, pro_mul_1_SCRIMP)

hFig = figure();
subplot(1,3,1);
    
plot(1:length(Varib2), Varib2, 'k-', 'LineWidth',1);
strTitle = 'Original Data';
title(strTitle);


subplot(1,3,2);
plot(pro_mul_1_AAMP, 'r-', 'LineWidth',1);
strTitle = 'Matrix Profile AAMP';
title(strTitle);

subplot(1,3,3);
plot(pro_mul_1_SCRIMP, 'b-', 'LineWidth',1);
strTitle = 'Matrix Profile Scrimp';
title(strTitle);

end
