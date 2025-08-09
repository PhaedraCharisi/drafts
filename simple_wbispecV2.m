function [wbis, wbic, sumwbic, freqs] = simple_wbispecV2(x, fs, wbicthres, doplot)
% HYBRID_WBISPEC Computes wavelet bispectrum with auxiliary scale alignment
%   x         : 1D signal (column vector)
%   fs        : sampling frequency
%   wbicthres : true to threshold bicoherence
%   doplot    : true to visualize

% Default arguments
if nargin < 3, wbicthres = false; end
if nargin < 4, doplot = false; end

% Ensure column vector
x = x(:);

% Normalize signal
x = x - mean(x);
x = x / max(abs(x));

% Frequency setup for PADS dataset
freqs = 2.5:0.5:18;
wname = 'cmor1-1';
fc = centfrq(wname);
scales = fc * fs ./ freqs;
N = length(scales);

% Compute wavelet coefficients
coefs_all = cwt(x, scales, wname);
for i = 1:N
    coefs_all(i, :) = coefs_all(i, :) / sqrt(scales(i));
end
coefs = coefs_all;

% Compute auxiliary scales for f1 + f2
scales_aux = zeros(N, N);
indx_aux = zeros(N, N);

for i = 1:N
    for j = 1:N
        s_aux = (scales(i) * scales(j)) / (scales(i) + scales(j));  % interaction scale
        scales_aux(i,j) = s_aux;
        [~, idx] = min(abs(scales - s_aux));  % find closest available scale
        indx_aux(i,j) = idx;
    end
end

% Build coefs_aux matrix (N x N x T)
T = size(coefs, 2);
coefs_aux = zeros(N, N, T);

for i = 1:N
    for j = 1:N
        coefs_aux(i,j,:) = coefs(indx_aux(i,j), :);
    end
end

% Prepare for bispectrum computation
wbis = zeros(N, N);
AsumSq = (abs(coefs).^2) * (abs(coefs)').^2;
coefs_conj_aux = conj(coefs_aux);

% Compute wavelet bispectrum
for t = 1:T
    z = coefs(:, t);
    for i = 1:N
        for j = 1:N
            wbis(i,j) = wbis(i,j) + z(i) * z(j) * coefs_conj_aux(i,j,t);
        end
    end
end

% Compute bicoherence
wbic = abs(wbis).^2 ./ (AsumSq + eps);
if wbicthres
    wbic(wbic < 0.01) = 0;
end

% Summed bicoherence
sumwbic = sum(wbic, 2);

% Optional plotting
if doplot
    t = (1:T) / fs;
    figure
    subplot(311); plot(t, x(1:T)); title('Input Signal'); xlabel('Time (s)');
    subplot(312); imagesc(t, freqs, abs(coefs)); axis xy;
    title('Scalogram'); xlabel('Time (s)'); ylabel('Freq (Hz)');
    subplot(313); plot(freqs, sumwbic); title('Summed Bicoherence'); xlabel('Frequency (Hz)');

    figure
    subplot(211); imagesc(freqs, freqs, abs(wbis)); axis xy; title('Wavelet Bispectrum'); colorbar;
    subplot(212); imagesc(freqs, freqs, abs(wbic)); axis xy; title('Wavelet Bicoherence'); colorbar;

    figure
    [F1, F2] = meshgrid(freqs, freqs);
    surf(F1, F2, abs(wbis), 'EdgeColor', 'none');
    xlabel('f₁ (Hz)'); ylabel('f₂ (Hz)'); zlabel('|WBIS|');
    title('3D Wavelet Bispectrum'); colorbar; view(45, 30);
end
end
