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
oldpath=('/Users/Amruta/Documents/MS GTCMT/Sem1/Research Project 7100/Audios_16kHz/');
% path='.\Audios_16kHz\';
file='Christina Perri - A Thousand Years.lrc';
[~,filename,~]=fileparts(file);
[Time_final,OffsetSec,Lyrics]= ReadLrc(oldpath,file);

TimeLyrics=Time_final+(ones(1,length(Time_final))*OffsetSec);
% find the locations of silences in the Lrc
emptyCells = cellfun(@isempty,Lyrics);
Idx=find(emptyCells==1);

load([oldpath filename '.mat']);

% using melodia SVD output:
MelodiaFile=[oldpath filename '_vamp_mtg-melodia_melodia_melody.csv'];
PitchMelodia=load(MelodiaFile);
PitchMelodia(PitchMelodia(:,2)<0,2)=0;

SmoothPitchMelodia=PitchMelodia;
windowSmooth=15;
ForCentr=floor(windowSmooth/2);
% smoothing of melody 
for i=ForCentr+1:length(PitchMelodia(:,2))-ForCentr
    temp=PitchMelodia(i-ForCentr:i+ForCentr,2);
    if all(temp)~=0 || sum(temp==0)>round(length(temp))*0.5;
        SmoothPitchMelodia(i,2)=median(temp);
    end
end

% TotalLines=length(Lyrics)-length(Idx);
TotalLines=length(Lyrics);

Nz_MP_vals=find(PitchMelodia(:,2)~=0);
StartLyrics=PitchMelodia(Nz_MP_vals(1),1);
EndLyrics=PitchMelodia(Nz_MP_vals(end),1);

[PM_ds]=MelodiaDownSample(MelodiaFile,0.01);

Nz_PM_ds=find(PM_ds(:,2)~=0);
NzPitchMelodia(:,1)=PM_ds(Nz_PM_ds,1);
NzPitchMelodia(:,2)=PM_ds(Nz_PM_ds,2);


NzPitchMelodiaCents=69+12*log2(NzPitchMelodia(:,2)/440);
plot(NzPitchMelodia(:,1),NzPitchMelodiaCents);

[sim_mat,nov_score]=SDM_nov(5,'Euclidean',NzPitchMelodiaCents);


Nz_nov_score=zeros(length(PM_ds),1);
Nz_nov_score(Nz_PM_ds)=nov_score;

figure;
imagesc(NzPitchMelodia(:,1),NzPitchMelodia(:,1),sim_mat); colormap gray; axis xy; axis square;
figure;
plot(NzPitchMelodia(:,1),nov_score); axis square;

figure;
subplot(211);imagesc(NzPitchMelodia(:,1),NzPitchMelodia(:,1),sim_mat);
subplot(212); plot(PM_ds(:,1),Nz_nov_score); 

figure;imagesc(1-(sim_mat/max(max(sim_mat)).^2));

% MFCCs
% hop=0.01; % hop in s
% frm=0.04; % frm in s
% AudioFileName=[oldpath filename '.wav'];
% AudioFile=audioread(AudioFileName);
% addpath('.\MFCC');
% coeff=melfcc(AudioFile,fs, 'maxfreq', 8000, 'numcep', 13, 'nbands', 40, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1, 'wintime',frm, 'hoptime', hop, 'preemph', 0, 'dither', 1);
% [sim_mat_coeff,nov_score_coeff]=SDM_nov(5,'Euclidean',coeff');
% figure;
% imagesc(PM_ds(:,1),PM_ds(:,1),sim_mat_coeff); %colormap gray;

% plot spectrogram, original lyrics and estimated lyrics
figure; title(filename);
% subplot(2,1,1);
imagesc(xTime,yFreq,1-20*log10(abs(SpecVals))); axis xy; colormap gray; %ylim([0 4000]); xlabel('Time (sec)'); ylabel('Freq (Hz)');
hold on;
plot(PitchMelodia(:,1),PitchMelodia(:,2),'.r','LineWidth',4); 
y1=get(gca,'ylim');

% subplot(2,1,2);
hold on;
% y1arr=ones(length(LyricsApproxTiming),1)*y1;
for itr=1:length(LyricsApproxTiming)
    hold on; line([LyricsApproxTiming(itr) LyricsApproxTiming(itr)],y1, 'LineWidth',2,'Color','c');
end
for itr=1:length(TimeLyrics)
    hold on; line([TimeLyrics(itr) TimeLyrics(itr)],y1, 'LineWidth',2);
end
% linkaxes('xy')

figure; title(filename);
% subplot(2,1,1);
imagesc(xTime,yFreq,1-20*log10(abs(SpecVals))); axis xy; colormap gray; %ylim([0 4000]); xlabel('Time (sec)'); ylabel('Freq (Hz)');
hold on;
plot(PitchMelodia(:,1),PitchMelodia(:,2),'.r','LineWidth',4); 
y1=get(gca,'ylim');
hold on;
for itr=1:length(TimeLyrics)
    hold on; line([TimeLyrics(itr) TimeLyrics(itr)],y1, 'LineWidth',2);
end
hold on;
for itr=1:length(Idx)
    hold on; line([TimeLyrics(Idx(itr)) TimeLyrics(Idx(itr))],y1, 'LineWidth',2,'Color','m');
end

% error calculation
thresh=0.5;
[prec, rec]= PrecRec(TimeLyrics, LyricsApproxTiming, thresh);
fmeasure=(2*prec*rec)/(prec+rec);


idx_lbl=label_segments(LyricsApproxTiming,PitchMelodia(:,1));
gt_lbl=label_segments(TimeLyrics,PitchMelodia(:,1));
[r_e,acp,r_a,asp,K]=clust_purity(idx_lbl,gt_lbl);

%display outputs
display(fmeasure);display(acp);display(asp);

toc
% using low level features:


% ApproxLocIntrvl=(EndLyrics-StartLyrics)/(TotalLines-1);
% ApproxLocs=StartLyrics:ApproxLocIntrvl:EndLyrics;
% 
% for iter=1:length(ApproxLocs)
%     [~,chk_loc]=min(abs(PitchMelodia(Nz_MP_vals,1)-ApproxLocs(iter)));
%     LyricsApproxTiming(iter)=PitchMelodia(Nz_MP_vals(chk_loc),1);
% end