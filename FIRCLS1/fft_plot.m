function Out =  fft_plot(Vector,N,FrequencySampling,description)
FFT = fft(Vector);
P2 = abs(FFT/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = FrequencySampling*(0:(N/2))/N;
Out = figure;
plot(f,P1,'LineWidth',2);
title('One-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend(description)
end