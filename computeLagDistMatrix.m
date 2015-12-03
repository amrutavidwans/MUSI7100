%% compute lag distance matrix 
% input:
%   SDM: numSamples by numSamples float matrix, self-distance matrix
% output:
%   R: numSamples by numSamples float matrix, lag-distance matrix
% Note: R should be a triangular matrix, xaxis = time, yaxis = lag
%       for more details, please refer to Figure 2 in the reference
%       "Paulus et al., Audio-based Music Structure Analysis, 2010"

function R = computeLagDistMatrix(SDM)

[rwSDM,clSDM]=size(SDM);
R=zeros(rwSDM,clSDM);
for i= 1:rwSDM
    for k = 1:clSDM
        if (i-k > 0)
            R(i, i-k) = SDM( i, k);  %% x-axis,column == lag, y-axis, row == time
        end
    end
end
R =  R';
end

