clear all;
close all;
clc;

SDM=imread('vertLines.png');
[rw,cl]=size(SDM);
factor=4;
SDM=imresize(SDM,[round(rw/factor),round(cl/factor)]);
% SDM=imrotate(SDM,45);

level = graythresh(SDM);
SDM_erode = im2bw(SDM',level);
imshow(SDM_erode);

BLyrics=imread('vertLines.png');

BLyrics_resize=imresize(BLyrics,[round(rw/factor),round(cl/factor)]);
BLyrics_resize=imrotate(BLyrics_resize,45);
level = graythresh(BLyrics_resize);
BW = im2bw(BLyrics_resize',level);
figure; imshow(BW);

LyricsResize=double(BW);
C=double(SDM_erode);
D=zeros(length(SDM_erode),length(SDM_erode));

for i=1:length(SDM_erode)
    for j=1:length(SDM_erode)
        CostMat=pdist2(LyricsResize(:,i),C(:,j));
%         [path,vals]=SimpleDtw(CostMat);
        [path, vals] = ToolSimpleDtw(CostMat);
        D(i,j)=vals(end,end);
    end
end
figure;imagesc(D);axis xy; 

figure; imagesc(LagC); axis xy;
figure; imagesc(LagLyricsResize); axis xy;
