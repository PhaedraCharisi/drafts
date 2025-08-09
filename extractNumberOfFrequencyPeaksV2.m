function no_freq_peaks = extractNumberOfFrequencyPeaksV2(signal,sampling_rate,pltspr)
% MAYBE THIS FUNCTION IS NOT NEEDED. (It's a modified version of a function
% from a 2023 adsp project)
    %This function calculates approximately the number of frequency peaks 
    %of THE SIGNAL (V2, Phaedra)?(a detailed coefficinet sample). To do that it calculates the 
    %envelope graph and smooths it using a moving average filter with 
    %window size equal to 10. Finally calculates the number of peaks in the
    %frequency domain taking into account the parameters minPeakDistance 
    %and minPeakPorminance. (V2: REMOVED ARGUMENT Period: Time Period Of The Detailed Coeff.)

    if nargin==2
        pltspr = false;
    end

    signal = signal.*hamming(length(signal));
    
    signal_fft = abs(fftshift(fft(signal)));
    signal_fft = signal_fft(floor(length(signal_fft)/2):end);

    freq_axis = linspace(0, sampling_rate/2, length(signal_fft));
    [~, max_idx] = max(signal_fft);
    dominant_freq = freq_axis(max_idx);
    if dominant_freq > 0
        period = 1/dominant_freq;
    else
        no_freq_peaks = 0;
        return
    end

    f0freq = 1/period;
    NoFreqSamples = f0freq*length(signal)/sampling_rate;

    localpeaks = {};
    localpeaksind = {};

    sig = signal_fft;

    windowsize = 10;
    sig = filter(ones(1,windowsize)/windowsize,1,sig);

    if length(sig)<3
        no_freq_peaks = 0;
        return
    end
        
    [localpeaks{end+1},localpeaksind{end+1}] = findpeaks(sig);

    if isempty(localpeaks{end})
        no_freq_peaks = 0;
        return
    end

    x = localpeaksind{end};
    sig = localpeaks{end};

    minpeakdistance = round(NoFreqSamples/5)*length(localpeaks{end})/length(signal_fft);
    minpeakprominence = 1.5*mean(signal_fft);

    if isnan(minpeakdistance) || length(sig)<3
        no_freq_peaks = 0;
        return
    end
    
    [localpeaks{end+1},localpeaksind{end+1}] = findpeaks(localpeaks{end},MinPeakDistance=minpeakdistance,MinPeakProminence=minpeakprominence);

    no_freq_peaks = numel(localpeaks{end});

    switch pltspr
        case 'true'
            figure
            nexttile
            plot(signal_fft)
            nexttile
            plot(x,localpeaks{end-1})
    end

end
