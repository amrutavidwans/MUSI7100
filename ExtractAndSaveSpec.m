clear all;
close all;
clc;

oldpath=('C:\Users\amrut_000\Documents\MS GTCMT\Sem1\Research Project 7100\Audios_16kHz\');
% path='.\Audios_16kHz\';

fnames = dir([oldpath '*.wav']);
numfids = length(fnames);


for iter=1:numfids
    wavFile=fnames(iter).name;
    display(wavFile);
    [~,wavName,~]=fileparts([oldpath wavFile]);
    [fwavVals,fs] = wavread([oldpath wavFile]);
    window=0.04; noverlap=0.75; nfft=(window*fs*noverlap*5);
    [SpecVals,yFreq,xTime]=spectrogram(fwavVals,(window*fs),(window*fs*noverlap),nfft,fs);
    save([oldpath wavName '.mat'],'wavFile','fwavVals','window','noverlap','nfft','SpecVals','yFreq','xTime','fs');
    clear 'wavName' 'fwavVals' 'window' 'noverlap' 'nfft' 'SpecVals' 'yFreq' 'xTime' 'fs';
end

% FileName_Pitch='Coldplay - Paradise_16k_vamp_mtg-melodia_melodia_melody.csv';
% 
% PitchMelodia = load([path FileName_Pitch],'r');  
% PitchMelodia(PitchMelodia(:,2)<0,2)=0;
% 
% load([path 'Coldplay - Paradise_16k.mat']);
% 
% %---------------------------For plotting PredomF0 on spectrogram
% figure;imagesc(xTime,yFreq,1-20*log10(abs(SpecVals))); axis xy; colormap gray; %ylim([0 4000]); xlabel('Time (sec)'); ylabel('Freq (Hz)');
% hold on;
% plot(PitchMelodia(:,1),PitchMelodia(:,2),'.r','LineWidth',4); 
% [xInput,yInput,button] = ginput(2);

fclose all;

