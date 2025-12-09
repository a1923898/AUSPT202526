% task 2: reverberant room 
% simulates reverberation by calculating reflections explicitly
% without using the audio toolbox 'acousticRoomResponse' function for
% paedagogical purposes
clc; clear all; 
% simulation parameters
c = 340;                    % speed of sound (m/s)
fs = 16000;                 % sampling rate (hz)
room_dim = [4.9, 4.9, 4.9]; % room dimensions [lx, ly, lz]
target_rt60 = 0.5;          % desired reverberation time

% calculate wall absorption (alpha) using sabiBeamformer (Broadband LCMV)ne's formula
v = prod(room_dim);
s = 2 * (room_dim(1)*room_dim(2) + room_dim(1)*room_dim(3) + room_dim(2)*room_dim(3));
alpha = (0.161 * v) / (s * target_rt60);
fprintf('calculated wall absorption (alpha): %.3f\n', alpha);

% calculate wall reflection coefficient (beta)
% beta is the fraction of amplitude remaining after one bounce.
% beta = sqrt(1 - alpha)
if alpha > 1, alpha = 0.99; end % safety cap
beta = sqrt(1 - alpha);

fprintf('calculated wall reflection coeff (beta): %.3f\n', beta);

% geometry setup
% mics
mic_pos = [2.41, 2.45, 1.5; ...
           2.49, 2.45, 1.5]'; % 3x2 matrix (x,y,z per column)

% sources
src_target_pos = [2.45; 3.45; 1.5];
src_interf_pos = [3.22; 3.06; 1.5];

% generate impulse responses
% define the order of reflections.
% order = 0 : direct only, anechoic
% order = 1 : 6 reflections off walls
% order = 10+ : thousands of reflections
% set the order high enough to ensure it captures all reflections that 
% arrive before the sound dies down up to the time determined by the rt60.
% max distance = c * rt60 ≈ 340*0.5= 170 
% => approx order≈ max distance / room dimension => order≈170/4.9≈34
refl_order = 34;

fprintf('generating manual rir for target (order %d)...\n', refl_order);
rir_target = roomImpulseResponse(src_target_pos, mic_pos, room_dim, beta, c, fs, refl_order);

fprintf('generating manual rir for interferer (order %d)...\n', refl_order);
rir_interf = roomImpulseResponse(src_interf_pos, mic_pos, room_dim, beta, c, fs, refl_order);

% normalize rirs to prevent clipping during convolution later
rir_target = rir_target ./ max(abs(rir_target(:)));
rir_interf = rir_interf ./ max(abs(rir_interf(:)));

% plot the generated impulse responses
figure;

plot((0:length(rir_target)-1)/fs, rir_target(:,1));
title('Reverberant Room Impulse Response (Mic 1)');
xlabel('time (s)'); ylabel('amplitude'); grid on;
figure
plot((0:length(rir_target)-1)/fs, rir_target(:,2));
title('Reverberant Room Impulse Response (mic 2)');
xlabel('time (s)'); ylabel('amplitude'); grid on;
% load audio and convolve
target_signal = 'target_signal.flac';
interference_signal = 'interference_signal1.flac';

% load and preprocess (same as before)
[sig_target, fs_t] = audioread(target_signal);
[sig_interf, fs_i] = audioread(interference_signal);

if size(sig_target, 2) > 1, sig_target = mean(sig_target, 2); end
if size(sig_interf, 2) > 1, sig_interf = mean(sig_interf, 2); end

sig_target = resample(sig_target, fs, fs_t);
sig_interf = resample(sig_interf, fs, fs_i);

n = min(length(sig_target), length(sig_interf));
sig_target = sig_target(1:n) / rms(sig_target(1:n));
sig_interf = sig_interf(1:n) / rms(sig_interf(1:n));

% apply the reverb by convolving the signals with the room impulse response
% fftfilt does the entire process and is executed in parallel on gpu but on
% smartphone would be better to use fft and ifft
fprintf('convolving audio with room impulse responses using fft\n');
target_at_mics = [fftfilt(rir_target(:,1), sig_target), fftfilt(rir_target(:,2), sig_target)];
interf_at_mics = [fftfilt(rir_interf(:,1), sig_interf), fftfilt(rir_interf(:,2), sig_interf)];

% trim to original length
target_at_mics = target_at_mics(1:n, :);
interf_at_mics = interf_at_mics(1:n, :);

% mix with conditions (same as task 1)
% SIR = 0db
p_t = mean(rms(target_at_mics).^2);
p_i = mean(rms(interf_at_mics).^2);
gain = sqrt(p_t / p_i);
interf_at_mics = interf_at_mics * gain;

clean_mix = target_at_mics + interf_at_mics;

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



diff_target = max(abs(rir_target(:,1) - rir_target(:,2)));
diff_interf = max(abs(rir_interf(:,1) - rir_interf(:,2)));
figure
plot((0:length(rir_target)-1)/fs, rir_interf(:,1));
title('interference rir mic 1 ');
xlabel('time (s)'); ylabel('amplitude'); grid on;
figure
plot((0:length(rir_target)-1)/fs, rir_interf(:,2));
title('interference rir mic 2 ');
xlabel('time (s)'); ylabel('amplitude'); grid on;

fprintf('difference between target mics:       %.4e (expected ~0 due to symmetry)\n', diff_target);
fprintf('difference between interference mics: %.4e (expected > 0)\n', diff_interf);

% play/save
audiowrite('mixture_signal.wav', mixture_signal, fs);
audiowrite('mixture_signal_10_snr.wav', snr10_final, fs);
audiowrite('mixture_signal_20_snr.wav', snr20_final, fs);
audiowrite('mixture_signal_40_snr.wav', snr40_final, fs);
audiowrite('mixture_signal_60_snr.wav', snr60_final, fs);
% save result without wgn for comparison
audiowrite('clean_mix.wav', clean_mix_final, fs);
% play mixture
fprintf('Playing mixture')
player1 = audioplayer(mixture_signal, fs);
playblocking(player1);
fprintf('saved audio as mixture_signal.wav, mixture_signal_10_snr.wav, mixture_signal_20_snr.wav, mixture_signal_40_snr.wav, mixture_signal_60_snr.wav, clean_mix.wav\n');

% Beamform
filtered_signal_GSC = GSC(mixture_signal, fs);
audiowrite('GSC_filtered.wav', filtered_signal_GSC, fs);
% play beamformed signal
fprintf('Signal after GSC beamforming')
player2 = audioplayer(filtered_signal_GSC, fs);
playblocking(player2);
% Frost Beamformer (Broadband LCMV)
mic_pos_toolbox = mic_pos; 

target_angle = 90; % Azimuth

[filtered_signal_FROST, w] = apply_frost_beamformer(mixture_signal, fs, mic_pos_toolbox, target_angle);

figure('Name', 'FROST Beamforming Results');
subplot(2,1,1);
plot(mixture_signal(:,1));
title('Noisy Reverberant Input (Mic 1)');
grid on;

subplot(2,1,2);
plot(filtered_signal_FROST);
title('Beamformed Output (Frost with Diagonal Loading)');
grid on;

% Save FROST result
audiowrite('result_FROST.wav', filtered_signal_FROST, fs);
fprintf('audio saved to result_FROST.wav\n');

% Reduce WG noise with Wiener filter
processed_signal = wiener(filtered_signal_GSC, fs);
% play denoised signal
fprintf('Playing final processed signal after Wiener decision-directed filter')
player3 = audioplayer(processed_signal, fs);
playblocking(player3);
% Save final result
audiowrite('processed_signal.wav', processed_signal, fs);
fprintf('final audio saved to processed_signal.wav\n');

% Run evaluation metrics
evaluationMetrics(processed_signal, fs);
fprintf('Try reducing the variable snr for the White Gaussian Noise SNR dB to see if there are improvements to the metrics')

% function to calculate echoic room impulse response
function rir = roomImpulseResponse(src_pos, mic_pos, room_dim, beta, c, fs, order)
    % inputs:
    %   src_pos: [x;y;z] source coordinates
    %   mic_pos: [x,y,z] mic coords (3xnummics)
    %   room_dim: [lx, ly, lz]
    %   beta: wall reflection coefficient (scalar)
    %   order: max reflection order (integer, e.g. 36)
    
    num_mics = size(mic_pos, 2);
    
    % estimate max length of rir based on order and room size
    % max distance approx = order * roomdiagonal
    max_dist = order * norm(room_dim) + norm(room_dim); % add room dim to prevent matrix index errors
    len_samples = ceil((max_dist/c) * fs);
    rir = zeros(len_samples, num_mics);
    
    lx = room_dim(1); ly = room_dim(2); lz = room_dim(3);
    
    % loop through virtual rooms
    % iterate through integer indices (nx, ny, nz)
    % nx=0 is the real room. nx=1 is the neighbor, etc.
    
    fprintf('calculating image sources up to order %d... ', order);
    
    for nx = -order:order
        for ny = -order:order
            for nz = -order:order
                
                % check reflection order for this specific image
                % order = |nx| + |ny| + |nz| (euclidian distance in room grid)
                current_order = abs(nx) + abs(ny) + abs(nz);
                if current_order > order
                    continue; 
                end
                
                % determine image source coordinates
                           
                % use ism coordinate formulas:
                rx = src_pos(1); ry = src_pos(2); rz = src_pos(3);
                
                % calculate image position              
                % 1. calculate base offset: [2*nx*lx, 2*ny*ly, 2*nz*lz]
                % 2. add source pos (if even) or mirrored source pos (if odd)
                
                vec_n = [nx, ny, nz];
                dim_vec = [lx, ly, lz];
                pos_vec = [rx, ry, rz];
                
                image_pos = zeros(3,1);
                
                for d = 1:3
                    n = vec_n(d);% Grid index (e.g., -1, 0, 1...)
                    l = dim_vec(d);% Dimension (Length/Width/Height)
                    p = pos_vec(d);% Original Source Coordinate                    
                    
                    % unified formula for image source coordinate:
                    % 1. n*L: moves us to the correct virtual room block
                    % 2. (-1)^n  and the rest handles the mirroring vs shifting logic
                    image_pos(d) = n*l + 0.5*l*(1 - (-1)^n) + p*(-1)^n;                    
                end
                
                
                % calculate contribution to each mic
                for m = 1:num_mics
                    dist = norm(image_pos - mic_pos(:,m));                    
                    % delay (samples)
                    d_samples = round((dist/c) * fs);
                    % check bounds (d_samples > 0)
                    if d_samples > 0 && d_samples <= len_samples
                        % amplitude attenuation
                        % 1. geometric spreading (1/distance)
                        % 2. wall absorption (beta^order)                        
                        gain = (1/dist) * (beta^current_order);                        
                        % add to accumulator
                        rir(d_samples, m) = rir(d_samples, m) + gain;
                    end
                end
            end
        end
    end
    fprintf('done.\n');
end