function [prec, rec]= PrecRec(orig, estim, thresh)

tp=0;
fp=0;
fn=0;

for lp=1:length(orig)
    [~,loc]=min(abs(estim-orig(lp)));
    if abs(estim(loc)-orig(lp))<=thresh
        orig(lp)=0;
        estim(loc)=0;
%         orig=orig(orig~=0);
        estim=estim(estim~=0);
        tp=tp+1;
    elseif abs(estim(loc)-orig(lp))>thresh
        orig(lp)=Inf;
%         orig=orig(orig~=0);
        fn=fn+1;
    end
end
fp=length(estim);
prec=tp/(tp+fp);
rec=tp/(tp+fn);

end