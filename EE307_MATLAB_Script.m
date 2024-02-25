% Course Code : EE 307
% Course Name : Communication Systems
% Group Number : 5
% Project Title : Pulse Width Modulation and Demodulation (PWM-DM)

% Clear Screen, Terminal and Workspace
close all;
clear all;
clc;

% PART 1 : PULSE WIDTH MODULATION (PWM)
time_limit = 7;
sampling_frequency = 2750;
time_step = 10^-5;

time = 0 : time_step : time_limit;
sawtooth_signal = sawtooth(2*pi*sampling_frequency*time);
message_signal = 0.75 * sin(2*pi*1*time);
array_length = length(sawtooth_signal);

pwm_signal = zeros(1,array_length);
for j = 1:array_length
    if (message_signal(j) >= sawtooth_signal(j))
        pwm_signal(j) = 1;
    end
end

figure();
p = plot(time, pwm_signal, 'r', time, message_signal, 'g' , time, sawtooth_signal, 'k');
p(1).LineWidth = 1.0;
axis([0 time_limit -1.5 1.5]);
ylabel('Amplitude');
xlabel('Time');
title('PWM Wave, Sampling Signal and Message Signal');
legend('Pulse Width Modulated Wave', 'Message Signal', 'Sampling (Sawtooth) Signal')
grid on;



% PART 2 - PULSE WIDTH DEMODULATION
pulse_train_signal = zeros(1,array_length);
ramp_and_hold_signal = zeros(1,array_length);
positive_edge = zeros();
threshold_value = [];

for l = 2:array_length
    if (pwm_signal(l) == 1) && (pwm_signal(l-1) == 0) && (l > 1)
        positive_edge = [positive_edge l];
    end
end
positive_edge(end) = [];

for m = positive_edge
    ramp_count_value = 1;
    while (pwm_signal(m + ramp_count_value) == 1)
        ramp_and_hold_signal(m + ramp_count_value) = (time_step)*ramp_count_value;
        ramp_count_value = ramp_count_value + 1;
    end
    val = (time_step)*ramp_count_value;
    threshold_value = [threshold_value val];
    while (pwm_signal(m + ramp_count_value) == 0)
        ramp_and_hold_signal(m + ramp_count_value) = val;
        ramp_count_value = ramp_count_value + 1;
    end
end

for k = positive_edge
    if (k > 10)
        pulse_train_signal(k-10:k-1) = max(threshold_value);
    end
end

pam_signal = pulse_train_signal + ramp_and_hold_signal - max(threshold_value);
pam_signal = max(pam_signal,0);

[numerator, denominator] = butter(5,0.02);
demod_signal = filter(numerator, denominator, pam_signal);
figure();
plot(time, demod_signal*4e3, 'r', time, message_signal, 'k');
axis([0 time_limit -1 1]);
ylabel('Amplitude');
xlabel('Time');
title('Demodulated Signal and Original Message Signal');
legend('Demodulated Signal', 'Original Message Signal');
grid on;