


%%parameter definition 
fs = 16e3;                                 % sampling rate = 16 kHz
room_dim = [4.9 4.9 4.9];                  % the room's width, length and height in meters
t_length_play=10;
c=340; % speed of sound in air m/s




%% audio sources
% Male Voice (Target)
[target_signal, fs_target] = audioread('1272-141231-0012.flac'); % random male speech from the resources provided, 
% feel free to replace with your own audio file however code currently assumes a stereo signal then converts to mono!!!!!!


if size(target_signal,2) > 1
    target_signal = mean(target_signal,2); %stereo to mono
end

if fs_target ~= fs
    target_signal = resample(target_signal, fs, fs_target); %resample at 16khHz
end
target_signal = target_signal / max(abs(target_signal));        % normalise
num_play_samples = min(numel(target_signal), t_length_play*fs); % trim samples to the minimum between length of the signal or 10 seconds
target_signal = target_signal(1:num_play_samples);               % trim to length

target_positions = [2.45, 3.45, 1.5];
target_angles = [90 0];                     % azimuth = 90 degree, elevation = 0 degree
% soundsc(target_signal(1:num_play_samples),fs) %plays original speech file once normalised

%Female Voice (Interference)
[interf_signal, fs_interf] = audioread('Prime_Minister_Gillard_of_Australia_at_News_Conference_with_President_Obama.flac');   % a female speech I downloaded from wikimedia, 
% might be too specific feel free to replace

if size(interf_signal,2) > 1
    interf_signal = mean(interf_signal,2); 
end

if fs_interf ~= fs
    interf_signal = resample(interf_signal, fs, fs_interf); %resample at 16khHz
end
interf_signal = interf_signal(1:num_play_samples); %cut length to same as target (IMPLIES interference longer than target!)
interf_signal = interf_signal / max(abs(interf_signal));

interf_position = [3.22 3.06 1.5];
interf_angles = [40 0];                     % azimuth = 40 degree, elevation = 0 degree
% pause(10);
% soundsc(interf_signal(1:num_play_samples),fs);
disp(length(target_signal));
disp(length(interf_signal));


%% room impulse response first stab:
mic_positions = [2.41 2.45 1.5; ...
                 2.49 2.45 1.5];

%ir = acousticRoomResponse(room_dim,target_positions,mic_positions(k,:),.....);
%if wanting to analyse behaviour of only the k-th mic ^^^^^^^^



ir_target = acousticRoomResponse(room_dim, target_positions, mic_positions, ...
    SampleRate=fs, SoundSpeed=c, ImageSourceOrder=0, ...
    NumStochasticRays=0, MaxNumRayReflections=0, ...
    MaterialAbsorption=1, MaterialScattering=0);

ir_interf = acousticRoomResponse(room_dim, interf_position, mic_positions, ...
    SampleRate=fs, SoundSpeed=c, ImageSourceOrder=0, ...
    NumStochasticRays=0, MaxNumRayReflections=0, ...
    MaterialAbsorption=1, MaterialScattering=0);

[num_mics, num_ir_samples] = size(ir_target); %room impulse response which I THINK is anechoic but need further comfirmation
t = (0:num_ir_samples-1) / fs; 



%% plot room impulse response
% subplot(1,2,1);
% plot(t, ir_target(1,:), 'b', 'LineWidth', 1.2);
% 
% xlabel("Time (s)");
% ylabel("Amplitude");
% title("Target Impulse Response – Mic 1 [2.41 2.45 1.5]");
% grid on;
% 
% subplot(1,2,2);
% plot(t, ir_target(2,:), 'r', 'LineWidth', 1.2);
% xlabel("Time (s)");
% ylabel("Amplitude");
% title("Target Impulse Response – Mic 2 [2.49 2.45 1.5]");
% grid on;
% 
% figure;
% subplot(1,2,1);
% plot(t, ir_interf(1,:), 'b', 'LineWidth', 1.2);
% 
% xlabel("Time (s)");
% ylabel("Amplitude");
% title("Interference Impulse Response – Mic 1 [2.41 2.45 1.5]");
% grid on;
% 
% subplot(1,2,2);
% plot(t, ir_interf(2,:), 'r', 'LineWidth', 1.2);
% xlabel("Time (s)");
% ylabel("Amplitude");
% title("Interference Impulse Response – Mic 2 [2.49 2.45 1.5]");
% grid on;


% %% plotting room i have no idea why this wont work???? 
% h = figure;
% plotRoom(room_dim,mic_positions,[target_positions;interf_position],h);
% 

mic_signals = zeros(num_mics, num_play_samples + num_ir_samples - 1);

for m = 1:num_mics
    target_contribution = conv(target_signal, ir_target(m,:), 'full');   % use microphone to "record" target signal
    
    interf_contribution = conv(interf_signal, ir_interf(m,:), 'full'); 

    mic_signals(m,:) = target_contribution + interf_contribution;
end

pause(10);
mic_signals = mic_signals(:, 1:num_play_samples);  % cut to play length
% soundsc(mic_signals(1,:), fs); % play mic1 recording

% 
% figure;
% subplot(2,1,1);
% plot(mic_signals(1,:)); title('Mic 1'); %received audio signal for mic 1
% 
% subplot(2,1,2);
% plot(mic_signals(2,:)); title('Mic 2'); %received audio signal for mic 2
% 
% 
% figure;
% 
% subplot(2,1,1);
% spectrogram(mic_signals(1,:), 256, 200, 512, fs, 'yaxis');
% title('Mic 1 Spectrogram (Target & Interference)');
% 
% subplot(2,1,2);
% spectrogram(mic_signals(2,:), 256, 200, 512, fs, 'yaxis');
% title('Mic 2 Spectrogram (Target & Interference)');



%%data file creation, be careful if modifying
simulation.fs = fs;
simulation.mic_signals = mic_signals;       
simulation.target_signal = target_signal;
simulation.interf_signal = interf_signal;
simulation.ir_target = ir_target;
simulation.ir_interf = ir_interf;
simulation.mic_positions = mic_positions;    
simulation.source_positions = [target_positions; interf_position];
fname = 'simData.mat';
save(fname, 'simulation', '-v7.3');




blms = dsp.BlockLMSFilter;


%[y, e, w] = blms(mic_signals(1,:)', target_signal)


% 
% [y,err] = ftf(mic_signals(1,:)', target_signal);

filt = dsp.FIRFilter;


filt.Numerator = designLowpassFIR(FilterOrder=10,CutoffFrequency=0.25)

ftf = dsp.FastTransversalFilter(10);


[y,e] = ftf(mic_signals(1,:)', target_signal);

soundsc(y,fs)












