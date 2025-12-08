% beamforming test

% load Data
filename = 'output_mixture.wav';
if ~isfile(filename)
    error('output_mixture.wav not present. Run simulation2.m first.');
end

[x, fs] = audioread(filename); 
% x is n x 2 matrix (time x mics)

% array geometry
d = 0.08;       % Mic spacing (m)
c = 340;        % Speed of sound (m/s)
theta_target = 90; % Target Angle (degrees)
theta_interf = 40;

% method 1: delay-and-sum beamformer
% Since the target is at 90 deg (broadside) and the mics are symmetric
% around the center, the TDOA (Time Difference of Arrival) is 0.
% Therefore, we simply average the channels.

y_das = mean(x, 2); % Simple average preserves target, reduces uncorrelated noise

% method 2: frequency-domain MVDR Beamformer
% calculate weights for every frequency bin that minimize total
% output power while ensuring the Target Direction has a gain of 1.

% STFT Parameters
win_len = 512;
overlap = 256;
fft_len = 512;

% perform STFT on input
[S, f, t] = stft(x, fs, 'Window', hamming(win_len), 'OverlapLength', overlap, 'FFTLength', fft_len);
% s dimensions: (Frequency x Time x Channels)

[num_bins, num_frames, num_mics] = size(S);
Y_mvdr_stft = zeros(num_bins, num_frames);

% loop through each frequency bin
for k = 1:num_bins
    freq = f(k);
    
    % create Steering Vector (v) for Target Direction (90 deg)
    % For a linear array on x-axis, delay tau = (d * cos(theta)) / c
    % At 90 degrees, cos(90) = 0, so delay is 0.
    % In general: v = exp(-1j * 2*pi * freq * tau_vector)
    % Here, v is simply [1; 1] because the target arrives simultaneously.
    v = ones(num_mics, 1);
    
    % extract snapshots for this frequency across all time frames
    % X_k: (Num_Mics x Num_Frames)
    X_k = squeeze(S(k, :, :)).'; 
    
    % estimate spatial covariance matrix (R)
    % average over the whole signal to find stationary interference statistics
    R = (X_k * X_k') / num_frames; 
    
    % add slight diagonal loading for numerical stability
    R = R + 0.001 * trace(R) * eye(num_mics);
    
    % calculate MVDR Weights
    % w = (R^-1 * v) / (v' * R^-1 * v)
    R_inv = inv(R);
    numerator = R\v;
    denominator = v' / R * v;
    w = numerator / denominator;
    
    % apply Weights to the signal at this frequency
    % output = w_hermitian * Input
    Y_mvdr_stft(k, :) = w' * X_k;
end

% inverse STFT to get time-domain signal
y_mvdr = istft(Y_mvdr_stft, fs, 'Window', hamming(win_len), 'OverlapLength', overlap, 'FFTLength', fft_len);

% normalize output volume to avoid clipping
y_mvdr = y_mvdr / max(abs(y_mvdr)) * 0.9;

% save audio
audiowrite('result_DAS.wav', y_das, fs);
audiowrite('result_MVDR.wav', y_mvdr, fs);

% plot time domain
t_vec = (0:length(x)-1)/fs;
% adjust length of MVDR (ISTFT padding)
min_len = min([length(x), length(y_das), length(y_mvdr)]);
y_mvdr = y_mvdr(1:min_len);
y_das = y_das(1:min_len);
x_ref = x(1:min_len, 1);

figure('Name', 'Beamforming Comparison');
subplot(3,1,1);
plot(t_vec(1:min_len), x_ref);
title('Raw Mixture'); grid on; ylabel('Amplitude');

subplot(3,1,2);
plot(t_vec(1:min_len), y_das);
title('DAS Output'); grid on; ylabel('Amplitude');

subplot(3,1,3);
plot(t_vec(1:min_len), y_mvdr);
title('MVDR Output'); grid on; ylabel('Amplitude');
xlabel('Time (s)');

% Spectrogram comparison
figure('Name', 'Spectrogram Analysis');
subplot(3,1,1);
spectrogram(x_ref, 256, 128, 256, fs, 'yaxis');
title('Mixture');

subplot(3,1,2);
spectrogram(y_das, 256, 128, 256, fs, 'yaxis');
title('Delay-and-Sum');

subplot(3,1,3);
spectrogram(y_mvdr, 256, 128, 256, fs, 'yaxis');
title('MVDR');

fprintf('1. "result_DAS.wav" created (Baseline).\n');
fprintf('2. "result_MVDR.wav" created (Advanced).\n');
