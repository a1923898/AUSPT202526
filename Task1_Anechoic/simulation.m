% anechoic chamber simulation
% simulation Parameters & geometry setup
c = 340;                % speed of sound (m/s)
fs = 16000;             % sampling rate (Hz)
room_dim = [4.9, 4.9, 4.9]; % room dimensions (meters)

% microphone positions (linear array)
% mic 1: (2.41, 2.45, 1.5)
% mic 2: (2.49, 2.45, 1.5)
mic_pos = [2.41, 2.45, 1.5; ...
    2.49, 2.45, 1.5]'; % 3x2 Matrix (x,y,z per column)

% source 1: target (male Speech) - azimuth 90 deg (front)
% position: (2.45, 3.45, 1.5)
src_target_pos = [2.45; 3.45; 1.5];

% source 2: interference (female Speech) - azimuth 40 deg (right-front)
% position: (3.22, 3.06, 1.5)
src_interf_pos = [3.22; 3.06; 1.5];

% load and pre-process audio sources
file_target = 'target_signal.flac';
file_interf = 'interference_signal1.flac';

% check if files exist
if ~isfile(file_target) || ~isfile(file_interf)
    error('Audio files not found. Please ensure the .flac files are in the current directory.');
end

% read audio
[sig_target, fs_t] = audioread(file_target);
[sig_interf, fs_i] = audioread(file_interf);

% convert to mono if necessary
if size(sig_target, 2) > 1, sig_target = mean(sig_target, 2); end
if size(sig_interf, 2) > 1, sig_interf = mean(sig_interf, 2); end

% resample to 16 kHz
sig_target = resample(sig_target, fs, fs_t);
sig_interf = resample(sig_interf, fs, fs_i);

% truncate to same length (use the shorter length)
N = min(length(sig_target), length(sig_interf));
sig_target = sig_target(1:N);
sig_interf = sig_interf(1:N);

% normalize inputs initially to unit power for consistency before propagation
sig_target = sig_target / rms(sig_target);
sig_interf = sig_interf / rms(sig_interf);

% simulate anechoic propagation (ISM Order = 0)
% in an ideal anechoic chamber (RT60 = 0), only the direct path exists.
% simulate using fractional delay filters and 1/distance attenuation.

fprintf('Simulating Target Propagation...\n');
target_at_mics = simulate_propagation(sig_target, src_target_pos, mic_pos, fs, c);

fprintf('Simulating Interference Propagation...\n');
interf_at_mics = simulate_propagation(sig_interf, src_interf_pos, mic_pos, fs, c);

% mixture conditions
% condition a: SIR = 0 dB (signal-to-interference ratio) at the microphones.
% calculate the power of the target and interference as received by the array
% and scale the interference to match the target power.

P_target = mean(rms(target_at_mics).^2); % Average power across mics
P_interf = mean(rms(interf_at_mics).^2);

% calculate scaling factor for SIR = 0 dB
desired_SIR_dB = 0;
gain_interf = sqrt(P_target / P_interf * 10^(-desired_SIR_dB/10));

% apply scaling
interf_at_mics_scaled = interf_at_mics * gain_interf;

% create the clean mix (signal + interference)
clean_mix = target_at_mics + interf_at_mics_scaled;

% condition B: SNR = 5 dB (signal-to-noise ratio)
% noise type: uncorrelated white gaussian noise (sensor noise)

% add different levels of WGN
% add white gaussian noise with snr 5
snr = 5;
raw_mixture = awgn(clean_mix, snr, 'measured');
snr10 = awgn(clean_mix, 10, 'measured');
snr20 = awgn(clean_mix, 20, 'measured');
snr40 = awgn(clean_mix, 40, 'measured');
snr60 = awgn(clean_mix, 60, 'measured');
noise = raw_mixture - clean_mix;

% normalise to prevent clipping
max_val = max(abs(raw_mixture(:)));
if max_val > 0.95
    scaling_factor = 0.95 / max_val;
    mixture_signal = raw_mixture * scaling_factor;
else
    mixture_signal = raw_mixture;
end

max_val = max(abs(clean_mix(:)));
if max_val > 0.95
    scaling_factor = 0.95 / max_val;
    clean_mix_final = clean_mix * scaling_factor;
else
    clean_mix_final = clean_mix;
end

max_val = max(abs(snr10(:)));
if max_val > 0.95
    scaling_factor = 0.95 / max_val;
    snr10_final = snr10 * scaling_factor;
else
    snr10_final = snr10;
end

max_val = max(abs(snr20(:)));
if max_val > 0.95
    scaling_factor = 0.95 / max_val;
    snr20_final = snr20 * scaling_factor;
else
    snr20_final = snr20;
end

max_val = max(abs(snr40(:)));
if max_val > 0.95
    scaling_factor = 0.95 / max_val;
    snr40_final = snr40 * scaling_factor;
else
    snr40_final = snr40;
end

max_val = max(abs(snr60(:)));
if max_val > 0.95
    scaling_factor = 0.95 / max_val;
    snr60_final = snr60 * scaling_factor;
else
    snr60_final = snr60;
end

% play mixture
fprintf('Playing mixture')
player1 = audioplayer(mixture_signal, fs);
playblocking(player1);
% output and visualization
t = (0:N-1)/fs;

figure('Name', 'Simulation Results (Mic 1)');
subplot(3,1,1);
plot(t, mixture_signal(:,1));
title('Mic 1: Final Mixture (Target + Interference + Noise)');
xlabel('Time (s)'); grid on;

subplot(3,1,2);
plot(t, target_at_mics(:,1));
title('Mic 1: Target Component Only');
xlabel('Time (s)'); grid on;

subplot(3,1,3);
plot(t, interf_at_mics_scaled(:,1));
title('Mic 1: Interference Component Only (Zoomed In)');
xlabel('Time (s)'); grid on;
ylim([-6 6]);
xlim([0.1 0.2]);
figure('Name', 'Simulation Results (Mic 2)');
subplot(3,1,1);
plot(t, mixture_signal(:,2));
title('Mic 2: Final Mixture (Target + Interference + Noise)');
xlabel('Time (s)'); grid on;

subplot(3,1,2);
plot(t, target_at_mics(:,2));
title('Mic 2: Target Component Only');
xlabel('Time (s)'); grid on;

subplot(3,1,3);
plot(t, interf_at_mics_scaled(:,2));
title('Mic 2: Interference Component Only (Zoomed In)');
xlabel('Time (s)'); grid on;
ylim([-6 6]);
xlim([0.1 0.2]);

figure('Name', 'Simulation Results (Both Mics)');
colororder("reef");
plot(t, mixture_signal(:,1));
hold on;
plot(t, mixture_signal(:,2));
title('Final Mixture (Target + Interference + Noise)');
legend('Mic 1', 'Mic 2');
xlabel('Time (s)'); grid on;

% play/save
audiowrite('mixture_signal.wav', mixture_signal, fs);
audiowrite('mixture_signal_10_snr.wav', snr10_final, fs);
audiowrite('mixture_signal_20_snr.wav', snr20_final, fs);
audiowrite('mixture_signal_40_snr.wav', snr40_final, fs);
audiowrite('mixture_signal_60_snr.wav', snr60_final, fs);
% save result without wgn for comparison
audiowrite('clean_mix.wav', clean_mix_final, fs);

fprintf('saved audio as mixture_signal.wav, mixture_signal_10_snr.wav, mixture_signal_20_snr.wav, mixture_signal_40_snr.wav, mixture_signal_60_snr.wav, clean_mix.wav\n');
% display calculated metrics
fprintf('\nSimulation Verification\n');
actual_SIR = 10*log10(mean(rms(target_at_mics).^2) / mean(rms(interf_at_mics_scaled).^2));
actual_SNR = 10*log10(mean(rms(clean_mix).^2) / mean(rms(noise).^2));
fprintf('Target SIR (should be ~0 dB): %.2f dB\n', actual_SIR);
fprintf('Target SNR (should be ~5 dB): %.2f dB\n', actual_SNR);

% Beamform
filtered_signal_GSC = GSC(mixture_signal, fs);
audiowrite('GSC_filtered.wav', filtered_signal_GSC, fs);




% play beamformed signal
fprintf('\nSignal after GSC beamforming')
player2 = audioplayer(filtered_signal_GSC, fs);
playblocking(player2);

% Reduce WG noise with Wiener filter
processed_signal = wiener(filtered_signal_GSC, fs);




% play denoised signal
fprintf('\nPlaying final processed signal after Wiener decision-directed filter')
player3 = audioplayer(processed_signal, fs);
playblocking(player3);

% Save final result
audiowrite('processed_signal.wav', processed_signal, fs);
fprintf('final audio saved to processed_signal.wav\n');

% Run evaluation metrics
evaluationMetrics(processed_signal, fs);

% play 40db SNR
fprintf('\nPlaying mixture with 40dB SNR')
player1 = audioplayer(snr40_final, fs);
playblocking(player1);
% output and visualization
t = (0:N-1)/fs;

% Beamform the one with 40dB SNR
filtered_signal_40dB = GSC(snr40_final, fs);

% play beamformed signal with 40 SNR
fprintf('\n40 dB Signal after GSC beamforming')
player2 = audioplayer(filtered_signal_40dB, fs);
playblocking(player2);

% Reduce WG noise of 40db SNR with Wiener filter
processed_signal_40db = wiener(filtered_signal_40dB, fs);

% play denoised signal from 40dB mix
fprintf('\nPlaying final processed signal after Wiener decision-directed filter on 40dB SNR mix\n')
player3 = audioplayer(processed_signal_40db, fs);
playblocking(player3);

fprintf('\nEvaluation metrics for the signal with 40db SNR\n')

% Run evaluation metrics
evaluationMetrics(processed_signal_40db, fs);

% helper Function: anechoic propagation
function mic_signals = simulate_propagation(src_sig, src_pos, mic_pos, fs, c)

num_mics = size(mic_pos, 2);
N = length(src_sig);
mic_signals = zeros(N, num_mics);

for m = 1:num_mics
    % calculate Distance
    dist = norm(src_pos - mic_pos(:,m));

    % calculate delay (samples)
    delay_sec = dist / c;
    delay_samples = delay_sec * fs;

    % calculate 1/r attenuation
    attenuation = 1 / dist;

    % apply fractional delay

    if exist('delayseq', 'file')
        delayed_sig = delayseq(src_sig, delay_samples);
    else
        % fallback: integer delay
        d_int = round(delay_samples);
        delayed_sig = [zeros(d_int,1); src_sig(1:end-d_int)];
    end

    % apply attenuation
    mic_signals(:, m) = delayed_sig * attenuation;
end
end