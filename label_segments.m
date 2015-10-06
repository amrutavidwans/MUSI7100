%%%%% Function label segments

% Given input : Peaks, timestamps
% Output : idx like output

function idx_label=label_segments(bndry,Tymlyrc)
    % Number of labels
    
    sz=size(bndry);
    nol=max(sz);
    
    lbl_ct=1;
%     flag=0;
    for i=1:length(Tymlyrc)
        
        if(lbl_ct<=nol)
            th=bndry(lbl_ct);
        else
            th=length(Tymlyrc);
        end
        
        if(Tymlyrc(i)<=th)
            idx_label(i)=lbl_ct;
        elseif Tymlyrc(i)~=nol
            lbl_ct=lbl_ct+1;
            idx_label(i)=lbl_ct;
        else
            idx_lbl(i)=lbl_ct;
        end
    end
    
end
            
        
    