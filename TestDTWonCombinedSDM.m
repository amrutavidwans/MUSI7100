clear all;
fclose all;
close all;
clc;

oldpath=('/Users/Amruta/Documents/MS GTCMT/Sem1/Research Project 7100/Audios_16kHz/');
% path='.\Audios_16kHz\';
file='Christina Perri - A Thousand Years.lrc';%'Alice Cooper - I''m Eighteen.lrc';
[~,filename,~]=fileparts(file);
[Time_final,OffsetSec,Lyrics]= ReadLrc(oldpath,file);

TimeLyrics=Time_final+(ones(1,length(Time_final))*OffsetSec);
% find the locations of silences in the Lrc
% emptyCells = cellfun(@isempty,Lyrics);
% Idx=find(emptyCells==1);

%% similarity of Lyrics
SimMatLyrics = LyricsSimilarity(Lyrics);
figure; imagesc(SimMatLyrics); axis xy;


%% similarity using SMtoolbox
addpath(genpath('/Users/Amruta/Documents/MS GTCMT/Sem1/Research Project 7100/MATLAB_SM-Toolbox_1.0/'));

oldpath=('/Users/Amruta/Documents/MS GTCMT/Sem1/Research Project 7100/Audios_16kHz/');

fname = [filename '.wav'];
[f_audio,sideinfo] = wav_to_audio('',oldpath, fname);
paramPitch.winLenSTMSP = 4410;
% paramPitch.visualize = 1;      % visualise pitch based feature inbuit in SM toolbox
[f_pitch] = audio_to_pitch_via_FB(f_audio,paramPitch);
paramCENS.winLenSmooth = 11;
paramCENS.downsampSmooth = 5;
% paramCENS.visualize = 1;      % visualize the chroma feature iin SMtoolbox
% paramCENS.featureRate = 2;
[f_CENS] = pitch_to_CENS(f_pitch,paramCENS);

S = features_to_SM(f_CENS,f_CENS);
paramVis.colormapPreset = 2;
paramVis.featureRate=2;
visualizeSM(S,paramVis);

paramVis.colormapPreset = 2;
paramSM.smoothLenSM = 20;
paramSM.tempoRelMin = 0.5;
paramSM.tempoRelMax = 2;
paramSM.tempoNum = 7;
paramSM.forwardBackward = 1;
paramSM.circShift = [0:11];
[S,I,parameter_ftosm] = features_to_SM(f_CENS,f_CENS,paramSM);

paramThres.threshTechnique = 1;
paramThres.threshValue = 0.85;
paramThres.applyBinarize = 1;
[S_thres,paramThres] = threshSM(S,paramThres); %(IdxStrt:IdxEnd,IdxStrt:IdxEnd)
visualizeSM(S_thres,paramVis);

TimeStamps= (0:(length(S_thres))-1)/paramVis.featureRate;

%% erode the SDM across diagonals
se = strel('line', 10, -45);
SDM_erode = imerode(S_thres,se); %(IdxStrt:IdxEnd,IdxStrt:IdxEnd) matrix considering IdxStrt and IdxEnd
figure;imagesc(SDM_erode); axis xy;
% hold on;
MarkGTonSDM(SDM_erode,TimeStamps,Lyrics,TimeLyrics); %TimeStamps(IdxStrt:IdxEnd)
hold off;
%% find the boundaries of connected components
[B,L] = bwboundaries(SDM_erode,'noholes');
% figure;imagesc(TimeStamps,TimeStamps,label2rgb(L, @jet, [.5 .5 .5]));  axis xy;
figure;imagesc(label2rgb(L, @jet, [.5 .5 .5]));  axis xy;
hold on
for k = 1:length(B)
boundary = B{k};
plot(boundary(:,1), boundary(:,2), 'w', 'LineWidth', 2)
end

%% filter out the edges with v low pixels inside it

for i=1:length(B)
chk_p(i)=length(B{i});
end
maxchkB=max(chk_p);
clear chk_p;

hardThreshBnd=round(0.05*maxchkB);
NewBnds={};
for i=1:length(B)
   if length(B{i})>hardThreshBnd && length(B{i})<maxchkB && B{i}(2,1)<B{i}(2,2)
      NewBnds=[NewBnds B{i}];
   end
end

NewBnds=NewBnds';

%% plot new vals
figure;imagesc(label2rgb(L, @jet, [.5 .5 .5]));  axis xy;
hold on
for k = 1:length(NewBnds)
boundary = NewBnds{k};
plot(boundary(:,1), boundary(:,2), 'w', 'LineWidth', 2)
end

%% take pixels start points

NewBndsOneSide=zeros(length(NewBnds),2);
for i=1:length(NewBnds)
    NewBndsOneSide(i,:)=NewBnds{i}(1,:);
end
% NewBndsOneSide=NewBndsOneSide';

%% pick the y co-ordinates of the start points with x co-ordinate same

ThresRep=2;
TempRep=[];
RepChk={};
CopyNewBndsOS=sortrows(NewBndsOneSide,1);
chk1=-1*ThresRep:ThresRep;
cnt1=0;
cnt2=1;
lenpst=length(CopyNewBndsOS);
lenprsnt=length(CopyNewBndsOS)-1;

while isempty(CopyNewBndsOS)~=1 && lenprsnt<lenpst && lenprsnt~=1
    
%     cnt2=0;
    
%     for i=1:length(CopyNewBndsOS)
        cnt3=1;
        tempstr=CopyNewBndsOS(1,1);
%         RepChk{cnt1}(cnt2,:)=NewBndsOneSide(i,:);
        TempRep=[];
        TempRep=[TempRep ; CopyNewBndsOS(1,:)];
%         CopyNewBndsOS(cnt1,:)=[];
        indx=[];
        indx(cnt3)=1;
        
        for j=2:length(CopyNewBndsOS)
            
%             if j~=i %&& abs(NewBndsOneSide(j,1)-2) < ThresRep 
                
%                 chk2=repmat(tempstr,1,length(chk1));
                chk=abs(tempstr-CopyNewBndsOS(j,1));
                if sum(chk<=ThresRep)>=1
                    cnt3=cnt3+1;
                    TempRep=[TempRep ; CopyNewBndsOS(j,:)];
                    indx(cnt3)=j;
                end
%             end
        end
        [rwTR,clTR]=size(TempRep);
        if isempty(TempRep)~=1 && rwTR>1
            cnt1=cnt1+1;
            RepChk{cnt1}=TempRep;
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenpst=rwcp;
            CopyNewBndsOS(indx,:)=[];
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenprsnt=rwcp;
        else 
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenpst=rwcp;
            CopyNewBndsOS(indx,:)=[];
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenprsnt=rwcp;
        end
end

% hold on
% for k = 1:length(RepChk)
% boundary = RepChk{k};
% plot(boundary(:,1), boundary(:,2), 'w', 'LineWidth', 2)
% end

%% keep only diagonals in the lyrics similarity
se = strel('line', 5, -45);
SDM_edLyrics = imerode(SimMatLyrics,se);
figure; imagesc(SDM_edLyrics); axis xy;

%% resize the Lyrics SDM to the size of feature SDM
[rw,cl]=size(SDM_erode);
BLyrics_resize = imresize(SDM_edLyrics,[rw cl]);
% figure; imagesc(BLyrics_resize); axis xy;

% convert to black and white
level = graythresh(BLyrics_resize);
BLyrics_resize=im2bw(BLyrics_resize,level);
se = strel('line', 10, -45);
BLyrics_resize = imerode(BLyrics_resize,se);
figure; imagesc(BLyrics_resize); axis xy;

% boundaries of the diagonals in Lyrics similarity matrix
[BLyrics,Lchk] = bwboundaries(BLyrics_resize,'noholes');


%% take boundaries from one side
chkLen=zeros(length(BLyrics),1);
for i=1:length(BLyrics)
    chkLen(i)=length(BLyrics{i});
end
maxchkB=max(chkLen);
clear chkLen;

NewBndsLyrics={};
hardThreshBnd=round(0.1*maxchkB);
for i=1:length(BLyrics)
   if length(BLyrics{i})>hardThreshBnd && length(BLyrics{i})<maxchkB && BLyrics{i}(2,1)<BLyrics{i}(2,2)
      NewBndsLyrics=[NewBndsLyrics BLyrics{i}];
   end
end
NewBndsLyrics=NewBndsLyrics';

figure;imagesc(label2rgb(Lchk, @jet, [.5 .5 .5]));  axis xy;
hold on
for k = 1:length(NewBndsLyrics)
boundary = NewBndsLyrics{k};
plot(boundary(:,1), boundary(:,2), 'w', 'LineWidth', 2)
end

%% experiment with direct multiplication of feature SDM and Lyrics SDM

% A=BLyrics_resize.*SDM_erode;
% figure;imagesc(A); axis xy;
% 
% [Achk,ALchk] = bwboundaries(A,'noholes');
% figure;imagesc(label2rgb(ALchk, @jet, [.5 .5 .5])); axis xy;
% hold on
% for k = 1:length(Achk)
% boundary = Achk{k};
% plot(boundary(:,1), boundary(:,2), 'w', 'LineWidth', 2)
% end

%% find out the centroids of the connected components in both Lyrics SDM and feature SDM

CCLyrics=zeros(length(NewBndsLyrics),2);
for i=1:length(NewBndsLyrics)
    CCLyrics(i,1)=mean(NewBndsLyrics{i}(:,1));
    CCLyrics(i,2)=mean(NewBndsLyrics{i}(:,2));
end
% hold on
% plot(CCLyrics(:,1), CCLyrics(:,2), 'b*')

CCfeat=zeros(length(NewBnds),2);
for i=1:length(NewBnds)
    CCfeat(i,1)=mean(NewBnds{i}(:,1));
    CCfeat(i,2)=mean(NewBnds{i}(:,2));
end
% hold on
% plot(CCfeat(:,1), CCfeat(:,2), 'b*')

% vector to centroids from the center of the SDMs
centFull=ones(size(CCLyrics)).*length(BLyrics_resize)/2;
CCLyricsVec=CCLyrics-centFull;
lenLyricsCC=length(centFull);
CCLyricsVecNew=zeros(size(CCLyricsVec));%zeros(length(CCLyricsVec),1)];

centFull=ones(size(CCfeat)).*rw/2;
CCfeatVec=CCfeat-centFull;
lenfeatCC=length(CCfeatVec);
CCfeatVecNew=zeros(size(CCfeatVec));%,zeros(length(CCfeatVec),1)];

% find magnitude and angle of all the vectors
CCLyricsVecNew(:,1)=sqrt(CCLyricsVec(:,2).^2+CCLyricsVec(:,1).^2);
CCLyricsVecNew(:,2)=atan2(CCLyricsVec(:,2),CCLyricsVec(:,1));

CCfeatVecNew(:,1)=sqrt(CCfeatVec(:,2).^2+CCfeatVec(:,1).^2);
CCfeatVecNew(:,2)=atan2(CCfeatVec(:,2),CCfeatVec(:,1));

% min cross is match
VecAngleDiff=zeros((size(CCLyricsVec,1)),(size(CCfeatVec,1)));
for i=1:(size(CCLyricsVec,1))
    for j=1:(size(CCfeatVec,1))
        VecAngleDiff(i,j)=1*(abs(CCLyricsVecNew(i,2)-CCfeatVecNew(j,2)))+0*(abs(CCLyricsVecNew(i,1)-CCfeatVecNew(j,1)));
    end
end

for i=1:(size(CCLyricsVec,1))
    [~,VecMatchLoc(i)]=min(VecAngleDiff(i,:));
end

%% keep centroids and vectors of the matched diagonals in SDM
MatchedCCVec=zeros(length(VecMatchLoc),2);
MatchedVec={};%zeros(length(VecMatchLoc),2);
MatchedNewBndsStrt=zeros(length(VecMatchLoc),2);

for i=1:(length(VecMatchLoc))
    MatchedCCVec(i,:)=CCfeat(VecMatchLoc(i),:);
    MatchedVec{i}=NewBnds{VecMatchLoc(i)};
    MatchedNewBndsStrt(i,:)=NewBndsOneSide(VecMatchLoc(i),:);
end

figure;imagesc(label2rgb(L, @jet, [.5 .5 .5]));  axis xy;
hold on
for k = 1:length(B)
boundary = B{k};
plot(boundary(:,1), boundary(:,2), 'w', 'LineWidth', 2)
end
hold on
plot(MatchedCCVec(:,1), MatchedCCVec(:,2), 'r*')

se = strel('line', 10, 45);
BLyrics_resize = imdilate(BLyrics_resize,se);
figure; imagesc(BLyrics_resize); axis xy;

 A=BLyrics_resize.*SDM_erode;
 figure;imagesc(A); axis xy;
 
 %% DTW on SDMs
 % B, BLyrics_resize SDM of SDM
% C=double(BLyrics_resize);
% E=SDM_erode.^2;
% D=zeros(576,576);
% for i=1:length(SDM_erode)
% for j=1:length(SDM_erode)
% D(i,j)=(sum((E(:,i))-(C(:,j))));
% end
% end
% figure;imagesc(D);axis xy;

LyricsResize=double(BLyrics_resize);
se = strel('line', 2, 45);
C=imdilate(SDM_erode,se);
C=double(C);
C=resample(C,1,6);
C=resample(C',1,6);
C=C';
LyricsResize=resample(LyricsResize,1,6);
LyricsResize=resample(LyricsResize',1,6);
LyricsResize=LyricsResize';

D=zeros(size(C));
LagC=computeLagDistMatrix(C);
LagLyricsResize=computeLagDistMatrix(LyricsResize);

figure; imagesc(LagC); axis xy;
figure; imagesc(LagLyricsResize); axis xy;

for i=1:length(LagC)
    for j=1:length(LagC)
        CostMat=pdist2(LagLyricsResize(:,i),LagC(:,j));
%         [path,vals]=SimpleDtw(CostMat);
        [path, vals] = ToolSimpleDtw(CostMat);
        D(i,j)=vals(end,end);
    end
end
figure;imagesc(D);axis xy; 

 %  D=pdist2(B,BLyrics_resize);
%  figure; imagesc(D); axis xy;
%  [p,C]=SimpleDtw(D);