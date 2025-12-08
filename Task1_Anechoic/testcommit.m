
%%FIR BP
bp = [300 3500] / (fs/2);
order = 160;
b = fir1(order, bp, hann(order+1));
out_FIR = filtfilt(b,1,aligned_sum); % zero-phase for offline
soundsc(out_FIR,fs)

%%IIR butter
[b_iir,a_iir] = butter(6, [300 3500]/(fs/2), 'bandpass');
out_IIR = filter(b_iir, a_iir, aligned_sum);
soundsc(out_IIR, fs)
freqz(b_iir,a_iir,2048,fs)


%%Wiener Post-Filter
%askdakjdaskdaskdbsadbask

