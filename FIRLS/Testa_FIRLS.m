close all
clear
clc

%signal fixed point = sfi(value,word length, fraction length)

%Frequency Sampling
FrequencySampling = 1.6e3;

%Time sampling
TimeSampling = 1/FrequencySampling;

%number of samples
N = 1024;   

%Time vector
t = 0:TimeSampling:((N-1)*TimeSampling);

FrequencySignal = 100;

%Signal = cos(2*pi*FrequencySignal*t);
%Signal = zeros(1,length(t));


%Signal frequency
FrequencySignal1 = 2;
FrequencySignal2 = 80;
FrequencySignal3 = 110;
G1 = 0.3;
G2 = 0.4;
G3 = 0.3;
Signal1 = cos(2*pi*FrequencySignal1*t);
Signal2 = cos(2*pi*FrequencySignal2*t);
Signal3 = cos(2*pi*FrequencySignal3*t);
Signal = G1*Signal1 + G2*Signal2 + G3*Signal3;


%}
%Inserting some noise and offset to the signal and changing the variable name
NoiseGain = 0.01;
Offset = 0.1;
Vector = Signal + Offset +  NoiseGain*randn(size(Signal));   
fft_plot(Vector,N,FrequencySampling,'Original Signal')

%Vector = Signal;
%Generating the theoretical output to perform the comparison between the signals 
IntegratedSignal = 0*G1*sin(2*pi*FrequencySignal1*t)/(2*pi*FrequencySignal1)+G2*sin(2*pi*FrequencySignal2*t)/(2*pi*FrequencySignal2)+G3*sin(2*pi*FrequencySignal3*t)/(2*pi*FrequencySignal3);

figure
plot(t,Vector,'k','LineWidth',3)
hold on
plot(t,Signal,'r','LineWidth',1)
legend('Noised Signal','Original Signal')

%filters (hide this section)

%Apply the highpass filter
CutFrequency1 = 5;
NormCOF =  (CutFrequency1*2)/FrequencySampling;
order = 90;
b = firls(order,[0 NormCOF NormCOF 1],[0 0 1 1]);
a = 1;
figure
freqz(b,a,512);                 % Frequency response of filter
a = 1;

figure
subplot(1,2,1)
plot(t,Vector,'k','LineWidth',3)
Vector =  filtfilt(b,a,Vector);

hold on
plot(t,Vector,'r','LineWidth',1)
legend('Noised Signal','Filtered Signal')

subplot(1,2,2)
plot(t,Vector,'k','LineWidth',3)
hold on
plot(t,Signal,'r','LineWidth',1)
legend('Filtered Signal','Original Signal')

fft_plot(Vector,N,FrequencySampling,'Filtered Signal')


VectorOut = integral_trap(Vector,TimeSampling);


fft_plot(VectorOut,N,FrequencySampling,'Integrated')

%Calculates the rms of the signal
rms2 = rms(Vector)
rmms = rms(VectorOut)

%plot the first signal
figure
subplot(1,2,1)
plot(t, Vector)
legend('m/s^2')
title('before integration')
%plot the second signal
subplot(1,2,2)
plot(t, VectorOut)
legend('m/s')
title('after integration')

%plot the theoretical signal and the integrated signal
figure



subplot(1,2,1)
plot(t, VectorOut,'LineWidth',2)
hold on
plot(t, IntegratedSignal,'LineWidth',1)
legend('Integrated','theoretical')
title('Theoretical Signal Vs Integrated Signal')

CutFrequency2 = 25;
NormCOF =  (CutFrequency2*2)/FrequencySampling;
b = firls(order,[0 NormCOF NormCOF 1],[0 0 1 1]);
a = 1;

FilteredVector =  filtfilt(b,a,VectorOut);

subplot(1,2,2)
plot(t, VectorOut,'LineWidth',2)
hold on
plot(t, FilteredVector,'LineWidth',1)
legend('Only integrated','Filtered integrated')
title('Integrated Signal Vs Filtered Integrated Signal')

fft_plot(FilteredVector,N,FrequencySampling,'Integrated and Filtered');
t_trans = N/16;
figure
plot(t, IntegratedSignal,'LineWidth',2)
hold on
plot(t, FilteredVector,'LineWidth',1)
legend('Theoretical','pratical')
title('Theoretical Signal Vs Pratical Signal')

max_error = max(abs(IntegratedSignal(t_trans:(end-t_trans)) - FilteredVector(t_trans:(end-t_trans))))
percentual_error = (100*max_error)/max(IntegratedSignal)

%Print the results on a file
%{
f=fopen('Testa_FIRLS.txt','w');
fprintf(f,'Resultados do Teste do Filtro FIRLS\n');
fprintf(f,'(OBS. Os parâmetros do filtro são os mínimos para o funcionamento desejado)\n\n');
fprintf(f,'Ordem do filtro: %.0f\n',order);
fprintf(f,'Técnica de filtragem: filtfilt\n');
fprintf(f,'Frequencia de corte pré integração: %.0f Hz\n',CutFrequency1);
fprintf(f,'Frequencia de corte pós integração: %.0f Hz\n',CutFrequency2);
fprintf(f,'Amplitude máxima do sinal teórico: %.15f\n',max(IntegratedSignal));
fprintf(f,'Erro máximo após integração (desconsiderando transitório): %.15f\n',max_error);
fprintf(f,'Erro máximo percentual após integração (desconsiderando transitório): %.2f %%\n',percentual_error);
fclose(f);

%}


%{
figure(1)
print -dpng fft_original_signal

figure(2)
print -dpng original_and_noised_signals

figure(3)
print -dpng filter_parameters_pre

figure(4)
print -dpng comparison_filtered_and_original_signals

figure(5)
print -dpng fft_after_first_filter

figure(6)
print -dpng fft_after_integration

figure(7)
print -dpng waveform_before_and_after_integration

figure(8)
print -dpng waveform_before_and_after_second_filter

figure(9)
print -dpng fft_integrated_and_filtered_signal

figure(10)
print -dpng comparison_between_theoretical_and_filtered_signal

%}




