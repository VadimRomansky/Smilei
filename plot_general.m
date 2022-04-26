clear;
general = importdata('./output_gamma0.3_sigma0.0002_theta80/scalars.txt');
N=size(general,1);

v = 0.9;
theta = 35*pi/180;
gamma = 1.0/sqrt(1.0 - v*v);
sigma = 0.004;
E1 = general(1,9)/general(1,6);

Eb = general(1,9)/(1 + v*v*sin(theta)*sin(theta));
Etot = general(1,6)*gamma/(gamma - 1);
measuredSigma = 2*Eb/Etot;

figure(1);
plot(general(1:N,1), general(1:N,4), 'red', general(1:N,1), general(1:N,6), 'green', general(1:N,1), general(1:N,9), 'blue');
title ('Utot');
legend('Utot', 'Ukin', 'Uelm');
xlabel ('t');
ylabel ('E');
grid;

figure(5);
plot(general(1:N,1), general(1:N,34), 'color', 'red');
title ('Exmax');
xlabel ('t');
ylabel ('E');
grid;

figure(6);
plot(general(1:N,1), general(1:N,46));
title ('Bxmax');
xlabel ('t');
ylabel ('B');
grid;