close all;
fclose all;
clear all;
clc;

%%% Extract Pitch Contour using Melodia of the wav files in folder
%%% .\Audios_16kHz and save them as .csv files in the same folder
% ExtractPitchUsingMelodia

%%% Extract the Spectrograms of the wav files in folder .\Audios_16kHz and
%%% save them as mat files with same name as the wav file in the same folder
% ExtractAndSaveSpec

% Read the Lrc files to get the ground truth
path='.\Audios_16kHz\';
file='Way Back Into Love.lrc';
[~,filename,~]=fileparts(file);
[Time_final,OffsetSec,Lyrics]= ReadLrc(path,file);

TimeLyrics=Time_final-(ones(1,length(Time_final))*OffsetSec);
% find the locations of silences in the Lrc
emptyCells = cellfun(@isempty,Lyrics);
Idx=find(emptyCells==1);

load([path filename '.mat']);

% using melodia SVD output:
PitchMelodia=load([path filename '_vamp_mtg-melodia_melodia_melody.csv']);
PitchMelodia(PitchMelodia(:,2)<0,2)=0;

% TotalLines=length(Lyrics)-length(Idx);
TotalLines=length(Lyrics);

Nz_MP_vals=find(PitchMelodia(:,2)~=0);
StartLyrics=PitchMelodia(Nz_MP_vals(1),1);
EndLyrics=PitchMelodia(Nz_MP_vals(end),1);

ApproxLocIntrvl=(EndLyrics-StartLyrics)/TotalLines;
ApproxLocs=StartLyrics:ApproxLocIntrvl:EndLyrics;

for iter=1:length(ApproxLocs)
    [~,chk_loc]=min(abs(PitchMelodia(Nz_MP_vals,1)-ApproxLocs(iter)));
    LyricsApproxTiming(iter)=PitchMelodia(Nz_MP_vals(chk_loc),1);
end

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
% using low level features:
