function [p,C]=SimpleDtw(D)

% p: path through matrix (Np x 2)
% C: cost matrix (same dimension as D)
% D: distance matrix (Nb x Na)

% [0,1,1,2,0;
%  1,2,0,1,1;
%  2,3,1,0,2;
%  1,0,2,3,1]

% initialization
[rwD,clD]=size(D);

CostMatrix=zeros(rwD,clD);
PathStore=cell(rwD,clD);

CostMatrix(1:2,1:2)=D(1:2,1:2);
CostMatrix(1,:)=cumsum(D(1,:));
CostMatrix(:,1)=cumsum(D(:,1));
PathStore{1,1}=[1,1];
PathStore{1,2}=[1,1];
PathStore{2,1}=[1,1];
for i=2:clD
    PathStore{1,i}=[1,i-1];
end
for i=2:rwD
    PathStore{i,1}=[i-1,1];
end 


% cost matrix computation
for i = 2:rwD
    for j=2:clD
        Val=zeros(1,3);
        Val(1,1)=D(i,j)+CostMatrix(i-1,j-1);
        Val(1,2)=D(i,j)+CostMatrix(i,j-1);
        Val(1,3)=D(i,j)+CostMatrix(i-1,j);
        [chk,loc]=min(Val);
        
        CostMatrix(i,j)=chk;
        switch loc
            case 1
                PathStore{i,j}=[i-1,j-1];
            case 2
                PathStore{i,j}=[i,j-1];
            case 3
                PathStore{i,j}=[i-1,j];
        end
        
    end
end

C=CostMatrix;

% back tracking for extracting the path indices
pstart=PathStore{rwD,clD};
p{1}=pstart;
for i=2:max(rwD,clD)
    chk=p{i-1};
    p{i}=PathStore{chk(1),chk(2)};
end

end