function MarkGTonSDM(SDM,TimeStamps,Lyrics,TimeLyrics)

figure;imagesc(TimeStamps,TimeStamps,SDM); axis xy;
hold on;
TimeLyrics=[TimeLyrics TimeStamps(end)];
TimeLyrics=round(TimeLyrics*10)/10;

% Lyrics{end+1}=[];

% for i=1:length(Lyrics)
%     if isempty(Lyrics{i})~=1    
%         for j=1:length(Lyrics)
%             if strcmp(Lyrics{i},Lyrics{j})==1
%                 hold on;
%                 wid=TimeLyrics(i+1)-TimeLyrics(i);
%                 rectangle('Position',[TimeLyrics(i) TimeLyrics(i+1) wid wid],'EdgeColor','r','LineWidth',2);
%             end
%         end
%     end
% end
k=1;
for i=1:length(Lyrics)
    if isempty(Lyrics{i})~=1    
        for j=1:length(Lyrics)
            if strcmp(Lyrics{i},Lyrics{j})==1
%                 hold on;
%                 chk=ones(1,length(TimeStamps));
%                 chk=chk.*TimeLyrics(j);
%                 [~,loc]=min(abs(TimeStamps-chk));
%                 for k=1:length(TimeStamps)
                    
                wid=TimeLyrics(i+1)-TimeLyrics(i);
                rectangle('Position',[TimeLyrics(j) TimeLyrics(i) wid wid],'EdgeColor','r','LineWidth',2);
%                 end
            end
        end
    end
end

% for i=1:length(TimeLyrics)-1
%     hold on;
%     if isempty(Lyrics{i})~=1
%         wid=TimeLyrics(i+1)-TimeLyrics(i);
%         rectangle('Position',[TimeLyrics(i) TimeLyrics(i+1) wid wid],'EdgeColor','r','LineWidth',1);
%     
% end


end