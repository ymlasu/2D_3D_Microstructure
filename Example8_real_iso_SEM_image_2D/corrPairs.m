function [corrVec,covVec] = corrPairs(samples,mu)
% corrPairs - non-Fourier implementation of sample autocorrelation
%
% Inputs:
% -------
%   samples       - nx by nSamples matrix of data
%   mu (optional) - mean of random field
%
% Outputs:
% --------
%   acf - vector of correlations

% Copyright Jack Pierce-Brown 2018

    [nLags,nSamples] = size(samples);

    if nargin < 2
%         mu = mean(mean(samples));
        mu = mean(samples,2);
    end
        
    samplesAdj = samples - mu;
% samplesAdj = samples;
    covVec = zeros(nLags,1);

    divisor = nLags;

%     for iLag = 0:(nLags-1)
%         sum1 = 0;
%         for iSample = 1:nSamples
%             for ix = 1:(nLags - iLag)
%                 sum1 = sum1 + samplesAdj(ix,iSample)*samplesAdj(ix+iLag,iSample);
%             end
%         end
% %         sum1 =  sum(sum(samplesAdj.^2));
% %         covVec(iLag+1) = sum1/(nSamples*divisor);
%  covVec(iLag+1) = sum1/(nSamples*(divisor-iLag));
%     end
%     corrVec = covVec/covVec(1);
    
    
    
    for iLag = 0:(nLags-1)
        
        covVec(iLag+1) = sum(sum(samplesAdj(1:end-iLag,:).*samplesAdj(iLag+1:end,:)))/(nSamples*(divisor-iLag));
    end
    corrVec = covVec/covVec(1);
    
end

