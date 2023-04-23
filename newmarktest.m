close all; 
clear all;
clc;


[nodes, elements] = processmeshfile('Situations/squareplate/squareplate21x21.msh');
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

Ka = assembleglobalaerodynamicmatrix(bcvector, nodes, elements);

[psi, b] = eigs(K, M, 7, 'smallestabs');
b = sqrt(b) / (2 * pi);

nyquistrate = 101.1856 * 2;

samplingrate = ceil(nyquistrate * 10);
dt = 1/samplingrate;

duration = 2;
duration  = ceil(duration * samplingrate) / samplingrate;
timesteps = samplingrate * duration;


lambda_dimensionless = 870;
lambda = lambda_dimensionless / ((1000/1000)^3 / (m.E * m.thickness^3/(12 * (1 - m.nu^2))));

[a, b] = eigs((K + lambda * Ka), M, 50, 'smallestabs');

% newu = zeros(size(nodes, 1) * 5);
newu = (real(a(:, 1))) * 100;

newudot = zeros(size(newu));
newuddot = zeros(size(newu));

ucollated = zeros(size(nodes, 1) * 3, timesteps);
ucollated(:, 1) = insertzero(newu, bcvector);

time = 0;
forcehistory = [];
for step=1:timesteps
    time = [time;step * dt];
    currentF = lambda * Ka * newu;
    forcehistory = [forcehistory, currentF];
    [newu, newudot, newuddot] = newmarkstep(newu, newudot, newuddot, K, M, currentF, dt);
    
    ucollated(:, step+1) = insertzero(newu, bcvector);
    disp(['step: ', num2str(step), '/', num2str(timesteps), ' max: ', num2str(max(abs(newu)))]);
end

ws = ucollated(1:3:end, :);
mag = 0.25 / max(max(abs(ws)));

figure;
view(45, 45);
axis([0, 1, 0, 1, -1, 1]);


figure('units','normalized','outerposition',[0 0 1 1]);
view(45, 45);
axis([0, 1, 0, 1, -1, 1]);

xlabel('x');
ylabel('y');
zlabel('z');

index = 1;
for step=3500:1:timesteps
%     p = plotnodes(nodes, ws(:, step), mag);
    p = plotmeshdeflection(nodes, elements, ucollated(:, step), mag);
    pause(0.01);
    disp(['time: ', num2str(time(step))]);
    frames(index) = getframe(gcf);
    delete(p);
    index = index + 1;
end

video = VideoWriter('animation.avi', 'Uncompressed AVI');
video.FrameRate = 60;
open(video);
writeVideo(video, frames);
close(video);

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