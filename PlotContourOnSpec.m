clear all;
close all;
clc;

oldpath=('E:\Users\Admin\Documents\MS GTCMT\Sem1\Research Project 7100\Audios_16kHz\');
% path='.\Audios_16kHz\';

fnames = dir([oldpath '*.mat']);
numfids = length(fnames);


for iter=1:numfids
    
    load([oldpath fnames(iter).name]);
    display(wavFile);
    [~,wavName,~]=fileparts([oldpath wavFile]);
    FileName_Pitch=[oldpath wavName '_vamp_mtg-melodia_melodia_melody.csv'];
    PitchMelodia = load(FileName_Pitch,'r');  
    PitchMelodia(PitchMelodia(:,2)<0,2)=0;
    
    %---------------------------For plotting PredomF0 on spectrogram
    figure; title(wavName);
    imagesc(xTime,yFreq,1-20*log10(abs(SpecVals))); axis xy; colormap gray; %ylim([0 4000]); xlabel('Time (sec)'); ylabel('Freq (Hz)');
    hold on;
    plot(PitchMelodia(:,1),PitchMelodia(:,2),'.r','LineWidth',4); 
    clearvars -except path fnames numfids iter;
    
end
    
fclose all;



% [xInput,yInput,button] = ginput(2);
    
% % create the plot of audio samples
% figure; hold on;
% plot(y, 'b'); % plot audio data
% title('Audio Data');
% xlabel(strcat('Sample Number (fs = ', num2str(Fs), ')'));
% ylabel('Amplitude');
% ylimits = get(gca, 'YLim'); % get the y-axis limits
% plotdata = [ylimits(1):0.1:ylimits(2)];
% hline = plot(repmat(0, size(plotdata)), plotdata, 'r'); % plot the marker
% 
% % instantiate the audioplayer object
% player = audioplayer(y, Fs);
% 
% % setup the timer for the audioplayer object
% player.TimerFcn = {@plotMarker, player, gcf, plotdata}; % timer callback function (defined below)
% player.TimerPeriod = 0.01; % period of the timer in seconds
% 
% % start playing the audio
% % this will move the marker over the audio plot at intervals of 0.01 s
% play(player);