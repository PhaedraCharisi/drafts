function [feat, names] = wbic_stats(wbic)
% MAYBE SOME FEATURES ARE NOT NEEDED. BUT THESE STATS COULD BE USED FOR STH
% ELSE? MAYBE THE RECONSTRUCTED SIGNAL

%WBIC_STATS Compact statistical summary from a WBIC vector.
%   [feat, names] = wbic_stats(wbic_vec) returns a small set of
%   descriptive statistics from the WBIC vector.
%
%   Features:
%     1) mean
%     2) standard deviation
%     3) skewness
%     4) kurtosis
%     5) median
%     6) interquartile range (IQR)
%     7) min
%     8) max
%
%   Input:
%     wbic_vec - 1xN or Nx1 numeric vector
%
%   Output:
%     feat  - 1x8 vector of feature values
%     names - 1x8 string array of feature names

    if ~isvector(wbic)
        error('wbic_vec must be a vector.');
    end

    x = wbic(:); % ensure column vector (so the code works the same whether you give it a row or column vector)

    mu   = mean(x);
    sig  = std(x);
    sk   = skewness(x);
    ku   = kurtosis(x);
    med  = median(x);
    q    = iqr(x);
    mn   = min(x);
    mx   = max(x);

    feat  = [mu, sig, sk, ku, med, q, mn, mx];
    names = ["mean","std","skewness","kurtosis",...
             "median","iqr","min","max"];
end