% DEMO
data = readmatrix('C:\Users\xaris\Documents\pads-parkinsons-disease-smartwatch-dataset-1.0.0\movement\timeseries\317_TouchNose_LeftWrist.txt', 'Delimiter', ',');
%data = readmatrix('C:\Users\xaris\Documents\pads-parkinsons-disease-smartwatch-dataset-1.0.0\movement\timeseries\060_RelaxedTask_RightWrist.txt', 'Delimiter', ',');
%data = readmatrix('C:\Users\xaris\Documents\pads-parkinsons-disease-smartwatch-dataset-1.0.0\movement\timeseries\060_Relaxed_LeftWrist.txt', 'Delimiter', ',');

signal = data(50:end,2); % e.g. acc_x
[wbis,wbic, sumwbic, freqs] = wbispecV2(signal,100,1,1); % (1->true, 0->false)
size(wbis)
size(wbic)
size(sumwbic)
size(freqs)

[a,b,c,d] = extract_prominent_wb_peaks(wbis);
h = wbispec_entropy(wbis);

s = total_sumwbic(sumwbic);

%%

data = readmatrix('C:\Users\xaris\Documents\pads-parkinsons-disease-smartwatch-dataset-1.0.0\movement\timeseries\060_HoldWeight_RightWrist.txt', 'Delimiter', ',');
data = data(50:end,:);

signal = data(:,2);
[wbis,wbic, sumwbic, freqs] = wbispecV2(signal,100,0,1); % (1->true, 0->false)

[a,b,c,d] = extract_prominent_wb_peaks(wbis);
h = wbispec_entropy(wbis);
h

[e,f,g,l] = extract_prominent_wb_peaks(wbic);
m = wbispec_entropy(wbic);
m
s = total_wbic(sumwbic);
s

[a,b,c,d] = extract_prominent_wb_peaksV2(wbis,1);

metrics = wb_pair_metrics(wbis) % [pairCount, pairFraction, pairEnergyShare, meanDiagDistance]