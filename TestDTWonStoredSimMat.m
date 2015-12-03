clear all;
close all;
clc;

% load('LagMatrices.mat');
factor=6;
% SDM=LagC;%
SDM=imread('TestImage.png');
SDM=rgb2gray(SDM);
SDM=double(SDM);

SDM=resample(SDM,1,factor);
SDM=resample(SDM',1,factor);
SDM=SDM';
[rw,cl]=size(SDM);


level = graythresh(SDM);
SDM_erode = im2bw(SDM,level);
imagesc(SDM_erode); axis xy;

% BLyrics=LagLyricsResize; 
BLyrics=imread('TestImageShifted.png');
BLyrics=rgb2gray(BLyrics);
BLyrics=double(BLyrics);

BLyrics=resample(BLyrics,1,factor);
BLyrics=resample(BLyrics',1,factor);
BLyrics=BLyrics';

level = graythresh(BLyrics);
BW = im2bw(BLyrics,level);
figure; imagesc(BW); axis xy;

LyricsResize=double(BW);
C=double(SDM_erode);
D=zeros(length(SDM_erode),length(SDM_erode));

for i=1:length(SDM_erode)
    for j=1:length(SDM_erode)
%         if i>=j
        CostMat=pdist2(LyricsResize(:,i),C(:,j));
%         [path,vals]=SimpleDtw(CostMat);
        [path, vals] = ToolSimpleDtw(CostMat);
        D(i,j)=vals(end,end);
%         end
    end
end
figure;imagesc(D);axis xy; xlabel('Lyrics'); ylabel('Audio');

[path, vals] = ToolSimpleDtw(D);