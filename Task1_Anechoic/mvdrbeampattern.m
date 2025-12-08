% MVDR Beam patterns

d = 0.08;           % Mic spacing
c = 340;            % Speed of sound
fs = 16000;
target_angle = 90-90;  % Target
interf_angle = 40-90;  % Interference

% Frequencies to plot
freq_list = round(linspace(200, 8000, 10)); 

% Load Audio
[x, fs] = audioread('output_mixture.wav');

% STFT Setup
fft_len = 1024;
[S, f_stft, t] = stft(x, fs, 'Window', hamming(1024), 'OverlapLength', 512, 'FFTLength', fft_len);

% Prepare Figure 
figure('Name', 'Multi-Frequency MVDR', 'Color', 'w', 'Position', [100, 100, 800, 700]);

% Create the Polar Axes object specifically
ax = polaraxes; 
hold(ax, 'on'); % Hold these specific polar axes

% Define colors for the lines (Blue -> Red)
colors = parula(length(freq_list));

% Mic positions relative to center
pos_mics = [-d/2; d/2];
angles = 0:1:360; 


% Loop through frequencies and plot
for i = 1:length(freq_list)
    freq = freq_list(i);
    
    % Calculate MVDR Weights for this freq
    [~, bin_idx] = min(abs(f_stft - freq));
    X_k = squeeze(S(bin_idx, :, :)).';
    
    % Covariance with regularization
    R = (X_k * X_k') / size(X_k, 2);
    R = R + 0.001 * trace(R) * eye(2); 
    
    % Steering Vector (Target 90 deg -> 0 delay)
    v_steer = [1; 1];
    
    % Weights
    R_inv = inv(R);
    w = (R_inv * v_steer) / (v_steer' * R_inv * v_steer);
    
    % Calculate Response Pattern
    response_pow = zeros(size(angles));
    for j = 1:length(angles)
        theta = angles(j);
        tau = (pos_mics * cosd(theta)) / c;
        v_test = exp(-1j * 2 * pi * freq * tau);
        response_pow(j) = abs(w' * v_test)^2;
    end
    
    % Normalization
    response_db = 10*log10(response_pow / max(response_pow)); % Norm to 0dB
    min_db = -40; % Floor
    response_db = max(response_db, min_db); % Clip
    
    rho_plot = response_db - min_db; % Shift positive for plotting
    
    % Plot Line
    polarplot(ax, deg2rad(angles), rho_plot, 'Color', colors(i,:), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('%d Hz', freq));
end

% Formatting
title(ax, 'MVDR Beam Patterns (200Hz - 8000Hz)');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'counterclockwise';

% Fix Axis Labels (Shift back to dB)
rticks = 0:10:40; 
ax.RTick = rticks;
ax.RTickLabel = string(rticks + min_db) + " dB";

% Add Source Markers
% plot them at radius 40 (which corresponds to 0 dB)
polarplot(ax, deg2rad(target_angle), 40, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9, 'DisplayName', 'Target');
polarplot(ax, deg2rad(interf_angle), 40, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9, 'DisplayName', 'Interference');

legend(ax, 'show', 'Location', 'eastoutside');
