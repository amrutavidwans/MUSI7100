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
file='Christina Perri - A Thousand Years.lrc';
[~,filename,~]=fileparts(file);
[Time_final,OffsetSec,Lyrics]= ReadLrc(oldpath,file);

TimeLyrics=Time_final+(ones(1,length(Time_final))*OffsetSec);
% find the locations of silences in the Lrc
emptyCells = cellfun(@isempty,Lyrics);
Idx=find(emptyCells==1);

% load([oldpath filename '.mat']);

AudioFileName=[oldpath filename '.wav'];
[AudioFile,fs]=audioread(AudioFileName);
addpath('.\Chroma');
% ChromA
fftlen=1024;
C = chromagram_IF(AudioFile,fs,fftlen);
% The frame advance is always one quarter of the FFT length.  Thus,
% the columns  of C are at timebase of fftlen/4/sr
tt = [1:size(C,2)]*fftlen/4/fs;

imagesc(tt,(1:12),20*log10(C+eps));
axis xy

[r,c]=size(C);
mn=repmat(mean(C,2),1,c);
sdev=repmat(std(C,0,2),1,c);
newcoeff= (C-mn)./sdev;

sim_mat_coeff=pdist2(newcoeff',newcoeff','Euclidean');

figure;
imagesc(tt,tt,sim_mat_coeff); axis xy; %colormap gray;

k = [0 1 2; -1 0 1; -2 -1 0];
H = conv2(double(sim_mat_coeff),k, 'same');
V = conv2(double(sim_mat_coeff),k','same');

E = sqrt(H.*H + V.*V);
threshold=0.5;
edgeImage = uint8((E > threshold) * 255);

figure; imagesc(tt,tt,1-edgeImage); colormap gray; axis xy;
% linkaxes('xy')

% figure; title(filename);
% % subplot(2,1,1);
% imagesc(xTime,yFreq,1-20*log10(abs(SpecVals))); axis xy; colormap gray; %ylim([0 4000]); xlabel('Time (sec)'); ylabel('Freq (Hz)');
% hold on;
% plot(SmoothPitchMelodia(:,1),SmoothPitchMelodia(:,2),'.r','LineWidth',4); 
% y1=get(gca,'ylim');
% hold on;
% for itr=1:length(TimeLyrics)
%     hold on; line([TimeLyrics(itr) TimeLyrics(itr)],y1, 'LineWidth',2);
% end
% hold on;
% for itr=1:length(Idx)
%     hold on; line([TimeLyrics(Idx(itr)) TimeLyrics(Idx(itr))],y1, 'LineWidth',2,'Color','m');
% end

% error calculation
thresh=0.5;
[prec, rec]= PrecRec(TimeLyrics, LyricsApproxTiming, thresh);
fmeasure=(2*prec*rec)/(prec+rec);


idx_lbl=label_segments(LyricsApproxTiming,SmoothPitchMelodia(:,1));
gt_lbl=label_segments(TimeLyrics,SmoothPitchMelodia(:,1));
[r_e,acp,r_a,asp,K]=clust_purity(idx_lbl,gt_lbl);

%display outputs
display(fmeasure);display(acp);display(asp);

toc
% using low level features:

% Nz_MP_vals=find(SmoothPitchMelodia(:,2)~=0);
% StartLyrics=SmoothPitchMelodia(Nz_MP_vals(1),1);
% EndLyrics=SmoothPitchMelodia(Nz_MP_vals(end),1);
% 
% Nz_PM_ds=find(SmoothPitchMelodia(:,2)~=0);
% NzPitchMelodia(:,1)=SmoothPitchMelodia(Nz_PM_ds,1);
% NzPitchMelodia(:,2)=SmoothPitchMelodia(Nz_PM_ds,2);
% 
% 
% NzPitchMelodiaMidi=69+12*log2(NzPitchMelodia(:,2)/440);
% plot(NzPitchMelodia(:,1),NzPitchMelodiaMidi);
% 
% [sim_mat,nov_score]=SDM_nov(5,'Euclidean',NzPitchMelodiaMidi);
% 
% 
% Nz_nov_score=zeros(length(SmoothPitchMelodia),1);
% Nz_nov_score(Nz_PM_ds)=nov_score;
% figure;
% imagesc(SmoothPitchMelodia(:,1),SmoothPitchMelodia(:,1),sim_mat); axis xy; axis square;
% figure;
% plot(SmoothPitchMelodia(:,1),nov_score); axis square;
% 
% figure;
% subplot(211);imagesc(NzPitchMelodia(:,1),NzPitchMelodia(:,1),sim_mat);
% subplot(212); plot(SmoothPitchMelodia(:,1),Nz_nov_score); 


% ApproxLocIntrvl=(EndLyrics-StartLyrics)/(TotalLines-1);
% ApproxLocs=StartLyrics:ApproxLocIntrvl:EndLyrics;
% 
% for iter=1:length(ApproxLocs)
%     [~,chk_loc]=min(abs(PitchMelodia(Nz_MP_vals,1)-ApproxLocs(iter)));
%     LyricsApproxTiming(iter)=PitchMelodia(Nz_MP_vals(chk_loc),1);
% end