close all;
fclose all;
clear all;
clc;

tic
%%% Extract Pitch Contour using Melodia of the wav files in folder
%%% .\Audios_16kHz and save them as .csv files in the same folder
% ExtractPitchUsingMelodia

%%% Extract the Spectrograms of the wav files in folder .\Audios_16kHz and
%%% save them as mat files with same name as the wav file in the same folder
% ExtractAndSaveSpec

% Read the Lrc files to get the ground truth
oldpath=('C:\Users\amrut_000\Documents\MS GTCMT\Sem1\Research Project 7100\Audios_16kHz\');
% path='.\Audios_16kHz\';
file='Blank Space.lrc';
[~,filename,~]=fileparts(file);
[Time_final,OffsetSec,Lyrics]= ReadLrc(oldpath,file);

TimeLyrics=Time_final+(ones(1,length(Time_final))*OffsetSec);
% find the locations of silences in the Lrc
emptyCells = cellfun(@isempty,Lyrics);
Idx=find(emptyCells==1);

% load([oldpath filename '.mat']);

% using melodia SVD output:

% TotalLines=length(Lyrics)-length(Idx);
TotalLines=length(Lyrics);

% MFCCs
hop=0.05; % hop in s
frm=0.2; % frm in s
AudioFileName=[oldpath filename '.wav'];
[AudioFile,fs]=audioread(AudioFileName);
addpath('.\MFCC');
coeff=melfcc(AudioFile,fs, 'maxfreq', 8000, 'numcep', 13, 'nbands', 40, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1, 'wintime',frm, 'hoptime', hop, 'preemph', 0, 'dither', 1);

% normalize the coeff
[r,c]=size(coeff);
mn=repmat(mean(coeff,2),1,c);
sdev=repmat(std(coeff,0,2),1,c);
newcoeff= (coeff-mn)./sdev;

sim_mat_coeff=pdist2(newcoeff',newcoeff','Euclidean');
TimeLen=length(AudioFile)/16000;
TimeInter=0:hop:TimeLen;

figure;
imagesc(TimeInter,TimeInter,sim_mat_coeff); axis xy; %colormap gray;

figure; imagesc(newcoeff); axis xy;

k = [0 1 2; -1 0 1; -2 -1 0];
H = conv2(double(sim_mat_coeff),k, 'same');
V = conv2(double(sim_mat_coeff),k','same');

E = sqrt(H.*H + V.*V);
edgeImage = uint8((E > threshold) * 255);

figure; imagesc(TimeInter,TimeInter,edgeImage); colormap gray; axis xy;

% subplot(2,1,2);
% y1arr=ones(length(LyricsApproxTiming),1)*y1;
% for itr=1:length(LyricsApproxTiming)
%     hold on; line([LyricsApproxTiming(itr) LyricsApproxTiming(itr)],y1, 'LineWidth',2,'Color','c');
% end
% for itr=1:length(TimeLyrics)
%     hold on; line([TimeLyrics(itr) TimeLyrics(itr)],y1, 'LineWidth',2);
% end
% % linkaxes('xy')
% 
% figure; title(filename);
% for itr=1:length(TimeLyrics)
%     hold on; line([TimeLyrics(itr) TimeLyrics(itr)],y1, 'LineWidth',2);
% end
% hold on;
% for itr=1:length(Idx)
%     hold on; line([TimeLyrics(Idx(itr)) TimeLyrics(Idx(itr))],y1, 'LineWidth',2,'Color','m');
% end
% 
% % error calculation
% thresh=0.5;
% [prec, rec]= PrecRec(TimeLyrics, LyricsApproxTiming, thresh);
% fmeasure=(2*prec*rec)/(prec+rec);
% 
% 
% idx_lbl=label_segments(LyricsApproxTiming,PitchMelodia(:,1));
% gt_lbl=label_segments(TimeLyrics,PitchMelodia(:,1));
% [r_e,acp,r_a,asp,K]=clust_purity(idx_lbl,gt_lbl);

%display outputs
display(fmeasure);display(acp);display(asp);

toc
% using low level features:
