clear;
general = importdata('./output/scalars.txt');
N=size(general,1);

figure(1);
plot(general(1:N,1), general(1:N,4));
title ('Utot');
xlabel ('t');
ylabel ('E');
grid;