close all; 
clear all;
clc;


[nodes, elements] = processmeshfile('Situations/squareplate/squareplate40x40.msh');
nodes = supportedbcsquareplate(nodes);


F = zeros(size(nodes, 1) * 3, 1);

for i=1:size(nodes, 1)
   F(3 * (i - 1) + 1) = 300000.8190239317/size(nodes, 1);
end



tic;
m = material(215*10^9 , 0.01, 0.3, 5/6, 8050);


[bcvector] = calculatebcvector(nodes);
[K] = assembleglobalstiffnessmatrix(m, bcvector, nodes, elements);
K = (K + K')/2;

[M] = assembleglobalmassmatrix(m, bcvector, nodes, elements);

M = (M + M')/2;

F = removerowcolumn(F, bcvector, 0);

nyquistrate = [195.377454586989] * 2;

samplingrate = ceil(nyquistrate * 50);
dt = 1/samplingrate;

duration = 0.3;
duration  = ceil(duration * samplingrate) / samplingrate;
timesteps = samplingrate * duration;



[psi, b] = eigs(K, M, 7, 'smallestabs');
b = sqrt(b) / (2 * pi);


newu = zeros(size(F));
newu = (psi(:, 1) + psi(:, 2) + psi(:, 3) + psi(:, 4) + psi(:, 5) + psi(:, 6)) * 100;

newudot = newu;
newuddot = newu;

ucollated = zeros(size(nodes, 1) * 3, timesteps);
ucollated(:, 1) = insertzero(newu, bcvector);

time = 0;

Machnumber = 5;
speedofsound = 303;

for step=1:timesteps
    time = [time;step * dt];
    currentF = F .* 0;
%     if(time(step + 1) < 0.25)
% %         currentF = F .* sin(195.377454586989 * 2 * pi * time(step+1));
%         currentF = F;
%     else
%         currentF = zeros(size(F));
%     end

    [newu, newudot, newuddot] = newmarkstep(newu, newudot, newuddot, K, M, currentF, dt);
    
    ucollated(:, step+1) = insertzero(newu, bcvector);
    disp(['step: ', num2str(step), ' max: ', num2str(max(abs(newu)))]);
end

ws = ucollated(1:3:end, :);
mag = 0.25 / max(max(abs(ws)));

% figure;
% view(45, 45);
% axis([0, 1, 0, 1, -1, 1]);
% for step=2400:timesteps
%     p = plotnodes(nodes, ws(:, step), mag);
%     pause(0.1);
%     disp(['time: ', num2str(time(step))]);
%     delete(p);
% end

%wsfreqs = ws(457, :).';
wsfreqs = ws.';

% samples = 1:dt:timesteps+1;
% samples = samples.';

figure;
plot(time, wsfreqs);

% figure;
% melSpectrogram(wsfreqs, samplingrate);


figure;
freqs = abs(fft(wsfreqs));
freqs = sum(freqs.').';
% freqs = freqs(:, 487);

freqs = freqs(1:timesteps/2);
f = samplingrate * (0:timesteps/2-1)/timesteps;
freqs = freqs / (timesteps/2);

plot(f, freqs);