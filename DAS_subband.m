% delay and sum with dyadic filter banks

% load data
filename = 'output_mixture.wav';
if ~isfile(filename)
    error('file output_mixture.wav not found. run simulation2.m first.');
end

[x, fs] = audioread(filename); 

% stft parameters (the "filter bank")
% use a smaller window for better time resolution, or larger for freq resolution.
window_len = 512;
overlap = 256;
fft_len = 512;

% perform stft (decomposition into frequency bands)
% s dimensions: (frequency x timeframes x channels)
[s, f, t] = stft(x, fs, 'window', hamming(window_len), 'overlaplength', overlap, 'fftlength', fft_len);

[num_bins, num_frames, num_mics] = size(s);

% design the frequency weights (the "sub-band" part)
% we want to apply different gains to different frequencies.
% male speech fundamental: ~85hz - 180hz
% key intelligibility: 300hz - 3400hz

weights = ones(num_bins, 1);

for k = 1:num_bins
    freq = f(k);
    
    % define bands
    if abs(freq) < 100
        % low rumble / dc -> cut significantly
        weights(k) = 0.1;
    elseif abs(freq) >= 100 && abs(freq) <= 3500
        % core speech range -> keep unity gain
        weights(k) = 1.0;
    else
        % high frequency (> 3.5khz) -> attenuate noise/hiss
        % since both male/female speech is here, but energy is lower,
        % we can dampen the high-frequency sensor noise.
        weights(k) = 0.2; 
    end
end

% plot the weight mask for visualization
figure;
plot(f, weights);
title('sub-band weighting mask');
xlabel('frequency (hz)'); ylabel('gain');
grid on; xlim([0 fs/2]);

% apply beamforming and weighting
% output container in frequency domain
y_stft = zeros(num_bins, num_frames);

for frame = 1:num_frames
    % extract mic signals for this time step
    % x_f is (numbins x nummics)
    x_f = squeeze(s(:, frame, :));
    
    % delay-and-sum beamforming
    % since we are broadside (90 deg), the delay is 0. 
    % average the complex coefficients which combines the signals constructively
    beamformed_bin = mean(x_f, 2); 
    
    % apply sub-band weights
    % apply the mask from above
    enhanced_bin = beamformed_bin .* weights;
    
    % store
    y_stft(:, frame) = enhanced_bin;
end

% reconstruction (inverse stft)
y_out = istft(y_stft, fs, 'window', hamming(window_len), 'overlaplength', overlap, 'fftlength', fft_len);

% normalize volume
y_out = y_out / max(abs(y_out)) * 0.9;

% ensure length matches original
n = min(length(x), length(y_out));
y_out = y_out(1:n);
x_ref = x(1:n, 1);

% visualization and save
figure('name', 'robust sub-band results');
subplot(2,1,1);
spectrogram(x_ref, 256, 128, 256, fs, 'yaxis');
title('original');

subplot(2,1,2);
spectrogram(y_out, 256, 128, 256, fs, 'yaxis');
title('sub-band filtered beamformer');

audiowrite('DAS_subband.wav', y_out, fs);
fprintf(' Saved "result_DAS_subband.wav".\n');