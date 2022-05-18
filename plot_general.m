clear;
general = importdata('./output/scalars.txt');
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


%note - numbers are shifted due number of species!!!!
startExField = 54;
figure(5);
hold on;
plot(general(1:N,1), general(1:N,startExField),'red');
plot(general(1:N,1), general(1:N,startExField+4),'green');
plot(general(1:N,1), general(1:N,startExField+8),'blue');
plot(general(1:N,1), general(1:N,startExField+12),'black');
plot(general(1:N,1), general(1:N,startExField+16),'magenta');
plot(general(1:N,1), general(1:N,startExField+20),'cyan');
title ('Max field');
legend('maxEx', 'maxEy', 'maxEz', 'maxBx', 'maxBy', 'maxBz');
xlabel ('t');
ylabel ('B');
grid;

startEnergy = 19;
figure(6);
hold on;
plot(general(1:N,1), general(1:N,startEnergy),'red');
plot(general(1:N,1), general(1:N,startEnergy+5),'green');
plot(general(1:N,1), general(1:N,startEnergy+10),'blue');
plot(general(1:N,1), general(1:N,startEnergy+15),'black');
plot(general(1:N,1), general(1:N,startEnergy+20),'magenta');
plot(general(1:N,1), general(1:N,startEnergy+25),'cyan');
title ('Species kin energy');
legend('protons', 'positrons cold', 'positrons hot', 'electrons right', 'electrons hot', 'electrons cold');
xlabel ('t');
ylabel ('E_k');
grid;