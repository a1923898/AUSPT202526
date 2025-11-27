%%parameter definition 
fs = 16e3;                                 % sampling rate = 16 kHz
room_dim = [4.9 4.9 4.9];                  % the room's width, length and height in meters
t_length_play=10;
c=340; % speed of sound in air m/s


num_play_samples = min(t_length_play*fs, numel(target_signal));


%% audio sources
% Male Voice (Target)
[target_signal, fs_file] = audioread('1272-141231-0012.flac'); % random male speech from the resources provided, 
% feel free to replace with your own however code currently assumes a stereo signal then converts to mono!!!!!!


if size(target_signal,2) > 1
    target_signal = mean(target_signal,2); 
end

if fs_file ~= fs
    target_signal = resample(target_signal, fs, fs_file); 
end
target_signal = target_signal / max(abs(target_signal));        % normalise
num_play_samples = min(numel(target_signal), t_length_play*fs); % trim samples to the minimum between length of the signal or 10 seconds
target_signal = target_signal(1:num_play_samples);               % trim to length

target_positions = [2.45 3.45 1.5];
target_angles = [90 0];                     % azimuth = 90 degree, elevation = 0 degree
soundsc(target_signal(1:num_play_samples),fs) %plays originak speech file once normalised


%% room impulse response first stab:
mic_positions = [2.41 2.45 1.5; ...
                 2.49 2.45 1.5];

%ir = acousticRoomResponse(room_dim,target_positions,mic_positions(k,:));
%if wanting to analyse behaviour of only the k-th mic ^^^^^^^^



ir = acousticRoomResponse(room_dim, target_positions, mic_positions, ...
    SampleRate=fs, SoundSpeed=c, ImageSourceOrder=0, ...
    NumStochasticRays=0, MaxNumRayReflections=0, ...
    MaterialAbsorption=1, MaterialScattering=0);
[num_mics, num_ir_samples] = size(ir); %room impulse response which I THINK is anechoic but need further investigation
t = (0:num_ir_samples-1) / fs; 



%% plot room impulse response
subplot(1,2,1);
plot(t, ir(1,:), 'b', 'LineWidth', 1.2);

xlabel("Time (s)");
ylabel("Amplitude");
title("Impulse Response – Mic 1 [2.41 2.45 1.5]");
grid on;

subplot(1,2,2);
plot(t, ir(2,:), 'r', 'LineWidth', 1.2);
xlabel("Time (s)");
ylabel("Amplitude");
title("Impulse Response – Mic 2 [2.49 2.45 1.5]");
grid on;

%% plotting room i have no idea why this wont work????
%h = figure;
%plotRoom(room_dim,mic_positions,target_positions,h);


mic_signals = zeros(num_mics, num_play_samples + num_ir_samples - 1);

for m = 1:num_mics
    mic_signals(m,:) = conv(target_signal, ir(m,:), 'full');   % use microphone to "record" signal
end

pause(10);
mic_signals = mic_signals(:, 1:num_play_samples);  % cut to play length
soundsc(mic_signals(1,:), fs); % play mic1 recording


figure;
subplot(2,1,1);
plot(mic_signals(1,:)); title('Mic 1'); %received audio signal for mic 1

subplot(2,1,2);
plot(mic_signals(2,:)); title('Mic 2'); %received audio signal for mic 2


figure;

subplot(2,1,1);
spectrogram(mic_signals(1,:), 256, 200, 512, fs, 'yaxis');
title('Mic 1 Spectrogram');

subplot(2,1,2);
spectrogram(mic_signals(2,:), 256, 200, 512, fs, 'yaxis');
title('Mic 2 Spectrogram');



%%data file creation, be careful if modifying
simulation.fs = fs;
simulation.mic_signals = mic_signals;       
simulation.target_signal = target_signal;    
simulation.ir = ir;                          
simulation.mic_positions = mic_positions;    
simulation.source_positions = target_positions;
fname = 'simData.mat';
save(fname, 'simulation', '-v7.3');
