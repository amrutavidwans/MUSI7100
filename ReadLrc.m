% clear all;
% close all;
% clc;

% Inputs:
% path='.\Audios_16kHz\';
%filename='Way Back Into Love.lrc'

% Outputs:
% Time_final=Time Stamps in sec of the lyrics
% OffsetSec= Offset to be added to the lyrics to align properly
% Lyrics= Cell array of lyrics corresponding to time stamps in Time_final

 
function [Time_final,OffsetSec,Lyrics]= ReadLrc(path,filename)

fid_lrc=fopen([path filename]);

tline_lrc = fgetl(fid_lrc);
Key   = 'offset';
Index = strfind(tline_lrc, Key);
Offset = sscanf(tline_lrc(Index(1) + length(Key)+1:end), '%g', 1);
OffsetSec=Offset/1000;

tline_lrc = fgetl(fid_lrc);
tline_lrc = fgetl(fid_lrc);
tline_lrc = fgetl(fid_lrc);
tline_lrc = fgetl(fid_lrc);

count=1;
while ischar(tline_lrc)
    % extract lyrics
    Key   = ']';
    Index = strfind(tline_lrc, Key);
    Lyrics{count} = tline_lrc(Index(1) + length(Key):end);
    
    %extract timing
    Key   = '[';
    Index = strfind(tline_lrc, Key);
    Timing_min= sscanf(tline_lrc(Index(1) + length(Key):end), '%f', 1);
    
    Key   = ':';
    Index = strfind(tline_lrc, Key);
    Timing_sec=sscanf(tline_lrc(Index(1) + length(Key):end), '%f', 1);
    Time_final(count)=Timing_min*60+Timing_sec;
    
    count=count+1;
    tline_lrc=fgetl(fid_lrc);
    
end

fclose(fid_lrc);

end