% task 4: generalized sidelobe canceller (gsc)
% a computationally efficient implementation of the mvdr beamformer.
% structure: 
%   -fixed beamformer (upper rail): preserves target.
%   -blocking matrix (lower rail): blocks target, collects interference.
%   -adaptive canceller: subtracts lower rail from upper rail.

% load data
filename = 'output_mixture.wav';
if ~isfile(filename)
    error('file output_mixture.wav not found. run simulation2.m first.');
end
[x, fs] = audioread(filename);

% array geometry & alignment
d = 0.08;           % mic spacing (m)
c = 340;            % speed of sound (m/s)
mic_locs = [-d/2, d/2];
theta_target = 90;  % target direction

% stft setup
win_len = 1024;
overlap = round(0.75 * win_len);
fft_len = 1024;

[s, f, t_stft] = stft(x, fs, 'window', hamming(win_len), 'overlaplength', overlap, 'fftlength', fft_len);
[num_bins, num_frames, num_mics] = size(s);

% output buffer
y_gsc_stft = zeros(num_bins, num_frames);


% frequency domain GSC processing
for k = 1:num_bins
    freq = f(k);
    
    % pre-steering (alignment)
    % we phase-shift the inputs so the target appears to come from 0 deg (broadside).
    % this simplifies the blocking matrix significantly.
    
    k_wave = 2*pi*freq/c;
    steering_vec = exp(1j * k_wave * mic_locs' * cosd(theta_target));
    
    % apply phase shift to align signals
    X_k = squeeze(s(k, :, :)).'; % (Mics x Time)
    X_aligned = X_k ./ steering_vec; % Element-wise division removes delay
    
    % the upper rail (fixed beamformer)
    % since signals are aligned, a simple average preserves the target.
    d_upper = mean(X_aligned, 1); % (1 x Time) - "desired signal + noise"
    
    % the lower rail (blocking matrix)
    % we need to block the target. since target signals are identical in x_aligned,
    % subtracting them (mic 1 - mic 2) cancels the target perfectly.
    % in gsc terms, blocking matrix b = [1; -1].
    
    % this contains only interference + noise (Target is gone)
    u_lower = X_aligned(1,:) - X_aligned(2,:); 
    
    % adaptive cancellation (unconstrained)
    % we want to find a weight 'w' such that: output = d_upper - w * u_lower
    % is minimized.
    % w = e[u * d'] / e[u * u']  (cross-variance / auto-variance)
    
    % calculate statistics (averaging over time frames)
    R_du = (u_lower * d_upper') / num_frames; % cross-correlation
    R_uu = (u_lower * u_lower') / num_frames; % auto-correlation
    
    % robustness: add tiny noise to denominator to avoid divide-by-zero
    R_uu = R_uu + 1e-6 * mean(abs(u_lower).^2); 
    
    % calculate scalar weight
    w_adapt = R_du / R_uu; 
    
    % step e: final subtraction
    Y_gsc_stft(k, :) = d_upper - (w_adapt' * u_lower);
end

% inverse STFT
y_gsc = istft(Y_gsc_stft, fs, 'Window', hamming(win_len), 'OverlapLength', overlap, 'FFTLength', fft_len);

% length matching
len = min(length(x), length(y_gsc));
y_gsc = y_gsc(1:len);
x_ref = x(1:len, 1);

% normalize
y_gsc = y_gsc / max(abs(y_gsc)) * 0.9;

% visualization
t_vec = (0:len-1)/fs;

figure('Name', 'GSC Beamformer Results');
subplot(3,1,1);
plot(t_vec, x_ref); title('Original 2 Mic input'); grid on; ylim([-1 1]);
subplot(3,1,2);
plot(t_vec, y_gsc); title('GSC Output'); grid on; ylim([-1 1]);

% spectrogram
subplot(3,1,3);
spectrogram(y_gsc, 256, 128, 256, fs, 'yaxis');
title('GSC Spectrogram');

audiowrite('result_GSC.wav', y_gsc, fs);
fprintf('Audio saved to "result_GSC.wav"\n');