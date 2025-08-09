function [wbis, wbic, sumwbic, freqs] = simple_wbispecV1(x, fs, wbicthres, doplot)
% SIMPLE_WBISPEC Computes wavelet bispectrum & bicoherence using nearest frequency matching
%   x         : 1D signal (e.g., acceleration on one axis)
%   fs        : sampling rate (e.g., 100 Hz)
%   wbicthres : true to apply thresholding
%   doplot    : true to visualize

% Default parameters
if nargin < 3, wbicthres = false; end
if nargin < 4, doplot = false; end

% Normalize input
x = x(:);  % ensure column vector
x = x - mean(x);
x = x / max(abs(x));

% Frequency range (for PADS: tremor etc.)
freqs = 2.5:0.5:18;
wname = 'cmor1-1';
scales = centfrq(wname) * fs ./ freqs;
N = length(freqs);

% Continuous Wavelet Transform
coefs = cwt(x, scales, wname);
for i = 1:N
    coefs(i,:) = coefs(i,:) / sqrt(scales(i));
end

% Prepare bispectrum matrices
wbis = zeros(N, N);
AsumSq = (abs(coefs).^2) * (abs(coefs)').^2;
coefs_conj = conj(coefs);

% Wavelet bispectrum using nearest frequency match
for t = 1:length(x)
    z = coefs(:, t);
    for i = 1:N
        for j = 1:N
            f_sum = freqs(i) + freqs(j);
            [~, k] = min(abs(freqs - f_sum));  % nearest matching frequency index
            wbis(i,j) = wbis(i,j) + z(i) * z(j) * coefs_conj(k, t);
        end
    end
end

% Wavelet bicoherence
wbic = abs(wbis).^2 ./ (AsumSq + eps);
if wbicthres
    wbic(wbic < 0.01) = 0;
end

% Summed bicoherence (useful feature)
sumwbic = sum(wbic, 2);

% Plot
if doplot
    t = (1:length(x)) / fs;
    figure
    subplot(311); plot(t, x); title('Signal'); xlabel('Time (s)');
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
