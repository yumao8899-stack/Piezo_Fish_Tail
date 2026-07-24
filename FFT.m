function [f,P1]=FFT(time_vector, displacement_signal)
displacement_signal=displacement_signal*1e3;
dt = time_vector(2) - time_vector(1); % 计算时间步长 (采样周期)
Fs = 1 / dt;                          % 采样频率 (Hz)
% 信号长度 L
L = length(displacement_signal);
% 对位移信号进行FFT
Y = fft(displacement_signal);
% FFT的结果是复数，包含幅值和相位。我们先计算其幅值。
% P2 是双边频谱，包含了正频率和负频率的信息。
P2 = abs(Y/L);

% 我们通常只关心正频率部分，所以取前半部分并乘以2（除了直流分量）, P1 是单边频谱。
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
% f 是与单边频谱 P1 对应的频率向量。
f = Fs * (0:(L/2)) / L;
end