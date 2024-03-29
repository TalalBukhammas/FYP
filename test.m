close all; 
clear all;
clc;


[nodes, elements] = processmeshfile('Situations/squareplate/squareplate21x21.msh');
% nodes = supportedbccircularplate(nodes);
nodes = supportedbcsquareplate(nodes);
nodes(:, 1) = nodes(:, 1) * 2;
figure();
hold on;
set(gcf, 'WindowState', 'maximized');
view(45, 45);
plotmesh(nodes, elements);


F = zeros(size(nodes, 1) * 3, 1);

for i=1:size(nodes, 1)
   F(3 * (i - 1) + 1) = 10 * 1000 * 2/size(nodes, 1);
end

% P = 1;
% elementnodevectors = zeros(8, 2);
% for i=1:size(elements, 1)
%     %Loop through nodes in each element
%     for n=1:8
%        elementnodevectors(n, :) = nodes(elements(i, n), 1:2); 
%     end
%     coords = [elementnodevectors(1, :);elementnodevectors(2, :);elementnodevectors(3, :);elementnodevectors(4, :);elementnodevectors(5, :);elementnodevectors(6, :);elementnodevectors(7, :);elementnodevectors(8, :)];
% 
%     localforce = FSDTForceVector(coords, P);
%     
%     for n=1:8
%        F(3 * (elements(i, n) - 1) + 1) = F(3 * (elements(i, n) - 1) + 1) + localforce(n); 
%     end
% end

%plotloads(nodes, F, 0.5, 20);

axis([0, 1, 0, 1, 0, 1]);
axis([-1, 1, -1, 1, 0, 1]);

daspect([1, 1, 1]);
xlabel('x');
ylabel('y');
zlabel('z');

hold off;

figure();
hold on;
set(gcf, 'WindowState', 'maximized');
view(45, 45);

tic;
%m = material(10920, 0.01, 0.3, 5/6, 1);
m = material(215*10^9 , 0.1, 0.3, 5/6, 8050);

[bcvector] = calculatebcvector(nodes);
[K] = assembleglobalstiffnessmatrixdof(m, bcvector, nodes, elements, 3);
K = (K + K')/2;

[M] = assembleglobalmassmatrix(m, bcvector, nodes, elements);
M = (M + M')/2;

F = removerowcolumn(F, bcvector, 0);
u = K \ F;
u = insertzero(u, bcvector);
F = insertzero(F, bcvector);
toc
ws = [];
for i=1:size(u, 1)/3
   ws = [ws;u(3 * (i - 1) + 1)];
end

mag = 1/max(ws) * 0.25;

plotmeshdeflection(nodes, elements, u, mag);

axis([0, 1, 0, 1, 0, 1]);
axis([-1, 1, -1, 1, 0, 1]);
daspect([1, 1, 1]);
xlabel('x');
ylabel('y');
zlabel('z');

figure();
hold on;
set(gcf, 'WindowState', 'maximized');
view(45, 45);
% 
% opts.maxit = 5000;
% opts.tol = 5000;
% [a, b] = eig(full(K), full(M));

[a, b] = eigs(K, M, 5, 'smallestabs');

mvalue = max(b);
eigendiagonal = zeros(size(b, 1), 2);
for i=1:size(eigendiagonal, 1)
    eigendiagonal(i, 1) = i;
    eigendiagonal(i, 2) = b(i, i);
end

eigendiagonal = sortrows(eigendiagonal, 2);

modeindex = eigendiagonal(5, 1);

omega = diag(sqrt(b));

nondimomega = omega(1) * sqrt(m.density/m.G);

modeshape1 = insertzero(a(:, modeindex), bcvector);

modeshape1ws = [];
for i=1:size(modeshape1, 1)/3
   modeshape1ws = [modeshape1ws;modeshape1(3 * (i - 1) + 1)];
end
mag = 1/max(abs(modeshape1ws)) * 0.25;

plotmeshdeflection(nodes, elements, modeshape1, mag)

axis([0, 1, 0, 1, 0, 1]);
daspect([1, 1, 1]);
xlabel('x');
ylabel('y');
zlabel('z');

hold off;
