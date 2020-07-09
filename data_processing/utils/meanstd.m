function [meanv, stdv] = meanstd(data,dim,true_val)

    if nargin < 2
        dim = 2;
    end

    if nargin < 3
        true_val = zeros(size(data));
    end
    
    % Mean values
    meanv = squeeze(nanmean(data, dim));
    size(meanv)
    
    % Std deviations
    stdv = squeeze(nanstd(data-true_val, 1, dim));

end