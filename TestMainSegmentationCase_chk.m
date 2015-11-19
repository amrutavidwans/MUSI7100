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

% load([oldpath filename '.mat']);

% using melodia SVD output:
MelodiaFile=[oldpath filename '_vamp_mtg-melodia_melodia_melody.csv'];
% PitchMelodia=load(MelodiaFile);
% PitchMelodia(PitchMelodia(:,2)<0,2)=0;

reqd_dur=0.04;
[PitchMelodia]=MelodiaDownSample(MelodiaFile,reqd_dur);

SmoothPitchMelodia=PitchMelodia;
windowSmooth=7;
ForCentr=floor(windowSmooth/2);
% smoothing of melody 
for i=ForCentr+1:length(PitchMelodia(:,2))-ForCentr
    temp=PitchMelodia(i-ForCentr:i+ForCentr,2);
    if all(temp)~=0 || sum(temp~=0)>round(length(temp))*0.5;
        SmoothPitchMelodia(i,2)=mean(temp);
    end
end

% TotalLines=length(Lyrics)-length(Idx);
TotalLines=length(Lyrics);

SmoothPitchMelodiaMidi(:,1)=SmoothPitchMelodia(:,1);
SmoothPitchMelodiaMidi(:,2)=69+12*log2(SmoothPitchMelodia(:,2)/440);

sim_mat=pdist2(SmoothPitchMelodia(:,2),SmoothPitchMelodia(:,2),'Euclidean');

figure;
imagesc(SmoothPitchMelodia(:,1),SmoothPitchMelodia(:,1),sim_mat); axis xy; axis square;

figure;imagesc(1-(sim_mat/max(max(sim_mat)).^2));

Nz_PM_ds=find(SmoothPitchMelodiaMidi(:,2)~=-Inf);
NzPitchMelodia(:,1)=SmoothPitchMelodiaMidi(Nz_PM_ds,1);
NzPitchMelodia(:,2)=SmoothPitchMelodiaMidi(Nz_PM_ds,2);
sim_mat=pdist2(NzPitchMelodia(:,2),NzPitchMelodia(:,2),'Euclidean');

figure;
imagesc(sim_mat); axis xy; axis square;


% plot spectrogram, original lyrics and estimated lyrics
% figure; title(filename);
% imagesc(xTime,yFreq,1-20*log10(abs(SpecVals))); axis xy; colormap gray; %ylim([0 4000]); xlabel('Time (sec)'); ylabel('Freq (Hz)');
% hold on;
% plot(SmoothPitchMelodia(:,1),SmoothPitchMelodia(:,2),'.r','LineWidth',4); 
% y1=get(gca,'ylim');

% subplot(2,1,2);
% hold on;
% y1arr=ones(length(LyricsApproxTiming),1)*y1;

% for itr=1:length(LyricsApproxTiming)
%     hold on; line([LyricsApproxTiming(itr) LyricsApproxTiming(itr)],y1, 'LineWidth',2,'Color','c');
% end
% for itr=1:length(TimeLyrics)
%     hold on; line([TimeLyrics(itr) TimeLyrics(itr)],y1, 'LineWidth',2);
% end
% % linkaxes('xy')
% 
% % error calculation
% thresh=0.5;
% [prec, rec]= PrecRec(TimeLyrics, LyricsApproxTiming, thresh);
% fmeasure=(2*prec*rec)/(prec+rec);
% 
% 
% idx_lbl=label_segments(LyricsApproxTiming,SmoothPitchMelodia(:,1));
% gt_lbl=label_segments(TimeLyrics,SmoothPitchMelodia(:,1));
% [r_e,acp,r_a,asp,K]=clust_purity(idx_lbl,gt_lbl);

%display outputs
% display(fmeasure);display(acp);display(asp);

toc


% ApproxLocIntrvl=(EndLyrics-StartLyrics)/(TotalLines-1);
% ApproxLocs=StartLyrics:ApproxLocIntrvl:EndLyrics;
% 
% for iter=1:length(ApproxLocs)
%     [~,chk_loc]=min(abs(PitchMelodia(Nz_MP_vals,1)-ApproxLocs(iter)));
%     LyricsApproxTiming(iter)=PitchMelodia(Nz_MP_vals(chk_loc),1);
% end