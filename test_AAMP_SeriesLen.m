% Compute the change in comutational time with the increase of time series
% length i.e. n
% Author information omitted for KDE review.
% For details of the code, see:
% "Efficient Matrix Profile Algorithms for Normalized and Non-Normalized Distances", submitted to KDE 2021.

clear
clc
close all;

subSequenceExpSheepData();
subSequenceExpProteinData();



function subSequenceExpProteinData()
clear

subSeqLen = 50;  % The length of the sub-sequence
load ('ProteinData.mat');

[keepAllTargetTogether, ~] = ConcatenateAllSeries(realData(1:100,:)); % concatenate all the time series

takenRws = 20;
while (takenRws < 150)
    fprintf('Protien Data : The numer of taken rows are  : %d \n', takenRws);
    
    apply_All_Algo(keepAllTargetTogether, subSeqLen);   % working for the 2nd column here
    takenRws = takenRws + 20;
end
end


function subSequenceExpSheepData()
clear

% Loading the sheep data
load('SheepDataFull.mat')


[keepAllTargetTogether, ~] = ConcatenateAllSeries(keepAllData(1:100,:)); % concatenate all the time series

subSeqLen = 50;  % The length of the sub-sequence
takenRws = 20;
while (takenRws < 150)
    fprintf('Sheep Data : The numer of taken rows are %d \n', takenRws);
    apply_All_Algo(keepAllTargetTogether, subSeqLen);   % working for the 2nd column here
    
    takenRws = takenRws + 20;
end
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