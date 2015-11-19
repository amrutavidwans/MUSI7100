% fclose all;
% close all;
clear all;
clc;

addpath(genpath('/Users/Amruta/Documents/MS GTCMT/Sem1/Research Project 7100/MATLAB_SM-Toolbox_1.0/'));

oldpath=('/Users/Amruta/Documents/MS GTCMT/Sem1/Research Project 7100/Audios_16kHz/');

filename = 'Christina Perri - A Thousand Years.wav';
% [f_audio,sideinfo] = wav_to_audio('', 'data_music/', filename);
[f_audio,sideinfo] = wav_to_audio('',oldpath, filename);
paramPitch.winLenSTMSP = 4410;
% paramPitch.visualize = 1;      % visualise pitch based feature inbuit in SM toolbox
[f_pitch] = audio_to_pitch_via_FB(f_audio,paramPitch);
paramCENS.winLenSmooth = 11;
paramCENS.downsampSmooth = 5;
% paramCENS.visualize = 1;      % visualize the chroma feature iin SMtoolbox
% paramCENS.featureRate = 2;
[f_CENS] = pitch_to_CENS(f_pitch,paramCENS);

S = features_to_SM(f_CENS,f_CENS);
paramVis.colormapPreset = 2;
paramVis.featureRate=2;
% visualizeSM(S,paramVis);

paramVis.colormapPreset = 2;
paramSM.smoothLenSM = 20;
paramSM.tempoRelMin = 0.5;
paramSM.tempoRelMax = 2;
paramSM.tempoNum = 7;
paramSM.forwardBackward = 1;
paramSM.circShift = [0:11];
[S,I,parameter_ftosm] = features_to_SM(f_CENS,f_CENS,paramSM);

paramThres.threshTechnique = 1;
paramThres.threshValue = 0.75;
paramThres.applyBinarize = 1;
[S_thres,paramThres] = threshSM(S,paramThres);
% visualizeSM(S_thres,paramVis);

% paramThres.threshTechnique = 2;
% paramThres.threshValue = 0.15;
% paramThres.applyBinarize = 0;
% paramThres.applyScale = 1;
% paramThres.penalty = -2;
% [S_final] = threshSM(S,paramThres);  
% paramVis.imagerange = [-2,1];
% paramVis.colormapPreset = 3;
% % paramVis.print=1;
% paramVis.figureName='SM_final';
% % paramVis.featureRate=2;
% handleFigure = visualizeSM(S_final,paramVis);
% 
% paramThres.threshTechnique = 1;
% paramThres.threshValue = 0.5;
% [S_thres_final2,paramThres2] = threshSM(S_final,paramThres);
% visualizeSM(S_thres_final2,paramVis);

TimeStamps= (0:(length(S_thres))-1)/paramVis.featureRate;

[B,L] = bwboundaries(S_thres,'noholes');
figure;imagesc(TimeStamps,TimeStamps,label2rgb(L, @jet, [.5 .5 .5])); axis xy;
figure;imagesc(label2rgb(L, @jet, [.5 .5 .5])); axis xy;
hold on
for k = 1:length(B)
boundary = B{k};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%% filter out the edges with v low pixels inside it

for i=1:length(B)
chk_p(i)=length(B{i});
end
maxchkB=max(chk_p);

hardThreshBnd=round(0.10*maxchkB);
NewBnds={};
for i=1:length(B)
   if length(B{i})>hardThreshBnd && length(B{i})<maxchkB && B{i}(1,1)<B{i}(1,2)
      NewBnds=[NewBnds B{i}];
   end
end

NewBnds=NewBnds';

%% plot new vals
figure;imagesc(label2rgb(L, @jet, [.5 .5 .5])); axis xy;
hold on
for k = 1:length(NewBnds)
boundary = NewBnds{k};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%% take pixels start points

% NewBndsOneSide={};
% for i=1:length(NewBnds)
%    for j=1:length(NewBnds{i})-1
%        if NewBnds{i}(j,1) <= NewBnds{i}(j+1,1)
%            NewBndsOneSide{i}(j,:)=NewBnds{i}(j,:);
%        end
%    
%    end
% end
% NewBndsOneSide=NewBndsOneSide';

NewBndsOneSide=zeros(length(NewBnds),2);
for i=1:length(NewBnds)
    NewBndsOneSide(i,:)=NewBnds{i}(1,:);
end
% NewBndsOneSide=NewBndsOneSide';

%% pick the y co-ordinates of the start points with x co-ordinate same

% ThresRep=2;
% TempRep=[];
% RepChk={};
% chk1=-ThresRep:ThresRep;
% cnt1=0;
% cnt2=1;

% for i=1:length(NewBndsOneSide)
%     cnt1=cnt1+1;
%     cnt2=1;
%     tempstr=NewBndsOneSide(i,1);
%     RepChk{cnt1}(cnt2,:)=NewBndsOneSide(i,:);
%     
%     for j=1:length(NewBndsOneSide)-1
%         
%         if j~=i %&& abs(NewBndsOneSide(j,1)-2) < ThresRep 
%             TempRep=[];
%             chk2=repmat(NewBndsOneSide(j,1),1,length(chk1));
%             chk=abs(chk1-chk2);
%             if sum(chk<=ThresRep)>=1
%                 TempRep=[TempRep ; NewBndsOneSide(j,:)];
%             end
%         end
%     end
%         [rwTR,clTR]=size(TempRep);
%         if isempty(TempRep)~=1
%             cnt2=cnt2+1;
%             RepChk{cnt}(cnt2,:)=NewBndsOneSide(j,:);
%         
%         
%     RepChk{cnt}=NewBnds{i}(1,:);
%         end
% end


ThresRep=2;
TempRep=[];
RepChk={};
CopyNewBndsOS=sortrows(NewBndsOneSide,1);
chk1=-1*ThresRep:ThresRep;
cnt1=0;
cnt2=1;
lenpst=length(CopyNewBndsOS);
lenprsnt=length(CopyNewBndsOS)-1;

while isempty(CopyNewBndsOS)~=1 && lenprsnt<lenpst && lenprsnt~=1
    
%     cnt2=0;
    
%     for i=1:length(CopyNewBndsOS)
        cnt3=1;
        tempstr=CopyNewBndsOS(1,1);
%         RepChk{cnt1}(cnt2,:)=NewBndsOneSide(i,:);
        TempRep=[];
        TempRep=[TempRep ; CopyNewBndsOS(1,:)];
%         CopyNewBndsOS(cnt1,:)=[];
        indx=[];
        indx(cnt3)=1;
        
        for j=2:length(CopyNewBndsOS)
            
%             if j~=i %&& abs(NewBndsOneSide(j,1)-2) < ThresRep 
                
%                 chk2=repmat(tempstr,1,length(chk1));
                chk=abs(tempstr-CopyNewBndsOS(j,1));
                if sum(chk<=ThresRep)>=1
                    cnt3=cnt3+1;
                    TempRep=[TempRep ; CopyNewBndsOS(j,:)];
                    indx(cnt3)=j;
                end
%             end
        end
        [rwTR,clTR]=size(TempRep);
        if isempty(TempRep)~=1 && rwTR>1
            cnt1=cnt1+1;
            RepChk{cnt1}=TempRep;
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenpst=rwcp;
            CopyNewBndsOS(indx,:)=[];
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenprsnt=rwcp;
        else 
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenpst=rwcp;
            CopyNewBndsOS(indx,:)=[];
            [rwcp,clcp]=size(CopyNewBndsOS);
            lenprsnt=rwcp;
        end
end

% hold on
% for k = 1:length(RepChk)
% boundary = RepChk{k};
% plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
% end

        
    
    
