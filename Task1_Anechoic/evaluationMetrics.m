% specify audio files to evaluate

processed_audio_file = 'result_WIENER_GSC_highnoise.wav'; 
reference_audio_file = '1272-141231-0012.flac'; 

% output text file for logging
dt = datetime("now", 'Format', 'yyyy-MM-dd HH:mm:ss');
log_file = strcat('evaluation_metrics_', processed_audio_file, '_', string(dt), '.txt');

% load and pre-process audio
fprintf('loading audio files\n');

if ~isfile(processed_audio_file) || ~isfile(reference_audio_file)
    error('audio files not found');
end

[sig_proc, fs_p] = audioread(processed_audio_file);
[sig_ref, fs_r] = audioread(reference_audio_file);

% ensure sampling rates match (16 kHz required for PESQ compliance)
target_fs = 16000;

if fs_p ~= target_fs
    sig_proc = resample(sig_proc, target_fs, fs_p);
end
if fs_r ~= target_fs
    sig_ref = resample(sig_ref, target_fs, fs_r);
end
fs = target_fs;

% ensure mono
if size(sig_proc, 2) > 1, sig_proc = mean(sig_proc, 2); end
if size(sig_ref, 2) > 1, sig_ref = mean(sig_ref, 2); end

% temporal alignment
% beamformers introduce delay and metrics assume signals are perfectly aligned.
% use cross-correlation to find the lag.

fprintf('aligning audio\n');
[c, lags] = xcorr(sig_ref, sig_proc);
[~, I] = max(abs(c));
lag = lags(I);

if lag > 0
    % reference is ahead of processed so pad processed
    sig_proc_aligned = [zeros(lag, 1); sig_proc];
    sig_ref_aligned = sig_ref;
elseif lag < 0
    % processed is ahead so pad reference (or cut processed)
    sig_ref_aligned = [zeros(-lag, 1); sig_ref];
    sig_proc_aligned = sig_proc;
else
    sig_proc_aligned = sig_proc;
    sig_ref_aligned = sig_ref;
end

% truncate to same length
len = min(length(sig_proc_aligned), length(sig_ref_aligned));
ref = sig_ref_aligned(1:len);
deg = sig_proc_aligned(1:len);

% scale invariance (normalize volume of degraded to match reference)
% project 'deg' onto 'ref' to find the best scaling factor alpha
alpha = (ref' * deg) / (deg' * deg);
deg = deg * alpha;

% calculate metrics
results = struct();

fprintf('calculating STOI...\n');
if exist('stoi', 'file')
    % stoi requires MATLAB audio toolbox
    results.STOI = stoi(ref, deg, fs);
else
    results.STOI = NaN;
    warning('STOI function not found or failed. Install MATLAB Audio Toolbox');
end

fprintf('attempting to calculate PESQ\n');
used_visqol = false;
if exist('pesq_mex', 'file')
    results.PESQ = pesq_mex(ref, deg, fs, 'both');
    results.Metric_Type = 'PESQ (MOS-LQO)';
else
    fprintf('PESQ not installed. Compile from https://github.com/ludlows/pesq-mex. Attempting ViSQOL...\n');
end
% ViSQOL (virtual speech quality objective listener)
results.VISQOL = visqol(ref, deg,fs);    

fprintf('Calculating OSINR...\n');
% calculate error signal
noise_residual = deg - ref;
p_signal = sum(ref.^2);
p_noise = sum(noise_residual.^2);
if p_noise == 0
    results.OSINR = Inf; 
else
    results.OSINR = 10 * log10(p_signal / p_noise);
end

fprintf('Calculating Log-Spectral Distortion (LSD)...\n');

% calculate lsd
results.LSD = calculate_lsd(ref, deg, fs);

% display results
clc;
fprintf('=================================================\n');
fprintf('EVALUATION RESULTS\n');
fprintf('Processed File: %s\n', processed_audio_file);
fprintf('=================================================\n');
fprintf('1. STOI (Intelligibility):  %.4f  (Range: 0-1)\n', results.STOI);
fprintf('2. PESQ narrowband (Quality):      %.4f  (Range: 1-4.5)\n', results.PESQ(1));
fprintf('2. PESQ wideband (Quality):      %.4f  (Range: 1-4.5)\n', results.PESQ(2));
fprintf('3. ViSQOL (Quality):      %.4f  (Range: 1-5)\n', results.VISQOL);
fprintf('4. OSINR (Accuracy):   %.2f dB (Higher is better)\n', results.OSINR);
fprintf('5. Timbre Distortion (LSD): %.2f dB (Lower is better)\n', results.LSD);
fprintf('=================================================\n');

% save
fid = fopen(log_file, 'w');
fprintf(fid, 'Evaluation Date: %s\n', datestr(now));
fprintf(fid, 'File Processed: %s\n\n', processed_audio_file);
fprintf(fid, 'Metrics:\n');
fprintf(fid, 'STOI:  %.4f (0-1)\n', results.STOI);
fprintf(fid, 'PESQ narrowband:  %.4f (1-4.5)\n', results.PESQ(1));
fprintf(fid, 'PESQ wideband:  %.4f (1-4.5)\n', results.PESQ(2));
fprintf(fid, 'ViSQOL:  %.4f (1-5)\n', results.VISQOL);
fprintf(fid, 'OSINR: %.2f dB (Higher is better)\n', results.OSINR);
fprintf(fid, 'LSD:   %.2f dB (Lower is better)\n ', results.LSD);
fclose(fid);
fprintf('Results saved to "%s".\n', log_file);

% log spectral distortion helper function
function lsd = calculate_lsd(ref, deg, fs)
    % STFT settings
    nfft = 2048;
    win = hamming(nfft);
    hop = nfft/2;
    
    [S_ref, ~] = spectrogram(ref, win, hop, nfft, fs);
    [S_deg, ~] = spectrogram(deg, win, hop, nfft, fs);
    
    % power spectrograms
    P_ref = abs(S_ref).^2 + 1e-10;
    P_deg = abs(S_deg).^2 + 1e-10;
    
    % log power difference
    log_diff = (10*log10(P_ref) - 10*log10(P_deg)).^2;
    
    % average over frequency (rows), then over time (cols)
    lsd = mean(sqrt(mean(log_diff, 1)));
end