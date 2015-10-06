% Function for cluster purity
function [r_e,acp,r_a,asp,K]=clust_purity(idx1,gt_lbl)

sz=length(idx1);
% gt_lbl=zeros(1,sz(1));
% 
% gt_lbl=label_segment_forTaan(gt,sz);

conf_mat=confusionmat(gt_lbl,idx1);

% conf_mat=conf_mat(1:2,:);

Na=max(gt_lbl);

Ne=max(idx1);

N=length(idx1);


% Code to find the ACP

r_e=zeros(1,Ne);
r_a=zeros(1,Na);


% Matrix for re
% j varying from 1 to Ne Estimated
% i varying from 1 to Na Annotated
for j=1:Ne
  for i=1:Na
        r_e(j)=r_e(j)+((conf_mat(i,j)^2)/(sum(conf_mat(:,j))^2));
  end
end

acp=0;
for j=1:Ne
    acp=acp+r_e(j)*sum(conf_mat(:,j));
end
acp=acp/N;

%%%%%%%%%%%%%%%%%%%%


for i=1:Na
  for j=1:Ne
        r_a(i)=r_a(i)+((conf_mat(i,j)^2)/(sum(conf_mat(i,:))^2));
  end
end

asp=0;
for i=1:Na
    asp=asp+r_a(i)*sum(conf_mat(i,:));
end
asp=asp/N;

  
K=sqrt(asp*acp);

end


% Code to find ASP
    