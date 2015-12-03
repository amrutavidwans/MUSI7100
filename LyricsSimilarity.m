% input: Lyrics: cell array of all the lines in the lyrics .lrc file
% output: SimMatLyrics: Similarity matrix of lyrics with each line compared
% with all other lines

function SimMatLyrics = LyricsSimilarity(Lyrics)
    
SimMatLyrics=zeros(length(Lyrics),length(Lyrics));

for i=1:length(Lyrics)
%     if isempty(Lyrics{i})~=1    
        for j=1:length(Lyrics)
            SimMatLyrics(i,j) = strcmp(Lyrics{i},Lyrics{j});
        end
%     end
end


end