% decision-directed wiener filter to remove residual white gaussian noise 
% load beamformed audio
input_file = 'result_GSC_zoom_0.wav'; 

if ~isfile(input_file)
    error('file %s not found. please run simulation2.m then gsc.m first .', input_file);
end

[y, fs] = audioread(input_file);

% ensure mono input (wiener filters are usually single-channel post-processors)
if size(y, 2) > 1
    y = y(:, 1); 
end

%  stft setup
win_len = 512;
overlap = 256; % increase to increase quality and computation time
fft_len = 512;

% convert to frequency domain
[Y, f, t_stft] = stft(y, fs, 'window', hamming(win_len), 'overlaplength', overlap, 'fftlength', fft_len);
[num_bins, num_frames] = size(Y);

% output buffer
s_est = zeros(num_bins, num_frames);

% noise estimation (static wgn assumption)
% since the noise is wgn, we can estimate the noise power
% from the first few frames (assuming speech hasn't started immediately).
init_frames = 6; % use first 6 frames (112 ms) for noise fingerprint
noise_pow = mean(abs(Y(:, 1:init_frames)).^2, 2); 

% flooring: ensure noise power isn't zero to avoid division errors
noise_pow = max(noise_pow, 1e-10);

alpha = 0.98;   % smoothing factor (0.98 is standard)
gain_min = 0.1; % gain floor - can try adjusting this to remove more noise

% initialize "previous clean estimate" for recursive step
S_prev_abs_sq = max(abs(Y(:,1)).^2 - noise_pow, 0); 

for n = 1:num_frames
    %  current noisy magnitude squared
    Y_curr = Y(:, n);
    Y_abs_sq = abs(Y_curr).^2;
    
    %  a-posteriori snr (gamma)
    % how much signal we see right now compared to the noise floor
    SNR_post = Y_abs_sq ./ noise_pow;
    
    %  a-priori snr (xi) - the decision directed step
    % formula: xi(n) = alpha * (g(n-1)^2 * gamma(n-1)) + (1-alpha) * max(gamma(n)-1, 0)

    SNR_prio = alpha * (S_prev_abs_sq ./ noise_pow) + (1 - alpha) * max(SNR_post - 1, 0);
           
    % compute wiener gain
    % G = xi / (xi + 1)
    G = SNR_prio ./ (SNR_prio + 1);
    
    % apply gain floor (prevent dead silence/artifacts)
    G = max(G, gain_min);
    
    % apply gain to signal
    % keep the phase of the noisy signal (y_curr), only modify magnitude.
    S_curr = G .* Y_curr;
    
    %  store for output and next iteration
    S_est(:, n) = S_curr;
    S_prev_abs_sq = abs(S_curr).^2; 
end

% inverse stft & output
y_clean = istft(S_est, fs, 'Window', hamming(win_len), 'OverlapLength', overlap, 'FFTLength', fft_len);

% length matching
len = min(length(y), length(y_clean));
y_clean = y_clean(1:len);
y_raw = y(1:len);

% normalize
y_clean = y_clean / max(abs(y_clean)) * 0.9;

% plot time domain
t = (0:len-1)/fs;
figure('Name', 'Wiener Filter Results');
subplot(2,1,1);
plot(t, y_raw); 
title('Beamformer Output with WGN'); grid on; ylim([-1 1]);

subplot(2,1,2);
plot(t, y_clean); 
title('Wiener Filter - Denoised'); grid on; ylim([-1 1]);

% plot spectrogram
figure('Name', 'Spectrogram Comparison');
subplot(2,1,1);
spectrogram(y_raw, 256, 128, 256, fs, 'yaxis');
title('Before Denoising');

subplot(2,1,2);
spectrogram(y_clean, 256, 128, 256, fs, 'yaxis');
title('After Denoising');

% metrics calculation
noise_residue = mean(y_clean(1:1600).^2); % measure noise in first 112 ms
noise_orig = mean(y_raw(1:1600).^2);
noise_reduction_db = 10*log10(noise_orig / noise_residue);

fprintf('Estimated Noise Reduction: %.2f dB\n', noise_reduction_db);

audiowrite('result_WIENER.wav', y_clean, fs);
fprintf('Audio saved to "result_WIENER.wav"\n');