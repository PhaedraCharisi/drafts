function [wbis, wbic, sumwbic, freqs] = simple_wbispecV0(x, fs, wbicthres, doplot)
% SIMPLE_WBISPEC Computes wavelet bispectrum & bicoherence for PADS data
%   x         : 1D signal (e.g., acceleration on one axis)
%   fs        : sampling rate (100 Hz for PADS)
%   wbicthres : true to apply thresholding
%   doplot    : true to visualize

% Default parameters
if nargin < 3, wbicthres = false; end
if nargin < 4, doplot = false; end

% Normalize input
x = x - mean(x);
x = x / max(abs(x));

% Frequency range for PADS (useful for tremor 3–12 Hz)
freqs = 2.5:0.5:18;     % THAT'S BETTER THAN 3-12Hz I think 3:0.5:12;
wname = 'cmor1-1';           % Morlet wavelet
scales = centfrq(wname) * fs ./ freqs;

% Continuous Wavelet Transform
coefs = cwt(x, scales, wname);

% Normalize wavelet energy
for i = 1:length(scales)
    coefs(i,:) = coefs(i,:) / sqrt(scales(i));
end

% Prepare bispectrum matrices
N = length(freqs);
wbis = zeros(N, N);
AsumSq = (abs(coefs).^2) * (abs(coefs)').^2;
coefs_conj = conj(coefs);

% Wavelet bispectrum accumulation
for t = 1:length(x)
    z = coefs(:, t);
    for i = 1:N
        for j = 1:N
            k = i + j;  % index for f_i + f_j
            if k <= N
                wbis(i,j) = wbis(i,j) + z(i) * z(j) * coefs_conj(k, t);
            end
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

% THAT'S BETTER, NEEDS ADJUSTING THOUGH
% for index1 = 1:numOfRows
%     sumwbic(index1) = sum(wbic(index1, 1:index1));
% end

% Plot
if doplot
    figure
    subplot(311); plot((1:length(x))/fs, x); title('Signal');
    subplot(312); imagesc((1:length(x))/fs, freqs, abs(coefs)); axis xy; title('Scalogram');
    subplot(313); plot(freqs, sumwbic); title('Summed Bicoherence');

    figure
    subplot(211); imagesc(freqs, freqs, abs(wbis)); axis xy; title('Wavelet Bispectrum'); colorbar;
    subplot(212); imagesc(freqs, freqs, abs(wbic)); axis xy; title('Wavelet Bicoherence'); colorbar;

    % 3D Wavelet Bispectrum
    figure
    [F1, F2] = meshgrid(freqs, freqs);
    surf(F1, F2, abs(wbis), 'EdgeColor', 'none');
    xlabel('f₁ (Hz)');
    ylabel('f₂ (Hz)');
    zlabel('|WBIS(f₁, f₂)|');
    title('3D Wavelet Bispectrum');
    colorbar;
    view(45, 30); % adjust viewing angle
end

end