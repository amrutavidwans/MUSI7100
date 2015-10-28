clc;
close all;
fclose all;
clear all;

%C:\Program Files (x86)\Sonic Visualiser

% to run for single file
%sonic-annotator.exe -d vamp:mtg-melodia:melodia:melody C:\audio\song.wav -w csv

%For processing all files in a directory (recursively):
%sonic-annotator.exe -d vamp:mtg-melodia:melodia:melody -r C:\audio -w csv

% command='sonic-annotator.exe -d vamp:mtg-melodia:melodia:melody -r "E:\Users\Admin\Documents\MS GTCMT\Sem1\Research Project 7100\Audios_16kHz" -w csv';
% status = system(command);
% oldpath=('C:\Users\amrut_000\Documents\MS GTCMT\Sem1\Research Project 7100\Audios_16kHz\');
% path='.\Audios_16kHz\';
path='C:\Users\amrut_000\Documents\MS GTCMT\Sem1\Research Project 7100\Audios_16kHz\';

fnames = dir([path '*.wav']);
numfids = length(fnames);


for iter=1:numfids
    wavFile=fnames(iter).name;
    display(wavFile);
    command=['sonic-annotator.exe -d vamp:mtg-melodia:melodia:melody -r "' path wavFile '" -w csv'];
    status = system(command);
end
%sonic-annotator.exe -d vamp:mtg-melodia:melodia:melody "E:\Users\Admin\Documents\MS GTCMT\Sem1\Research Project 7100\Coldplay - Paradise_16k.wav" -w csv