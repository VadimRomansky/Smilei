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

%note - different in different smilei versions
Nspecies = 2;


%note - numbers are shifted due number of species!!!!
startExField = 16 + 5*Nspecies + 6 + 3;
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

startEnergy = 20;
figure(6);
hold on;
plot(general(1:N,1), general(1:N,startEnergy),'red');
plot(general(1:N,1), general(1:N,startEnergy+5),'green');
if(Nspecies > 2)
    plot(general(1:N,1), general(1:N,startEnergy+10),'blue');
end;
if(Nspecies > 3)
    plot(general(1:N,1), general(1:N,startEnergy+15),'black');
end;
if(Nspecies > 4)
    plot(general(1:N,1), general(1:N,startEnergy+20),'magenta');
end;
if(Nspecies > 5)
    plot(general(1:N,1), general(1:N,startEnergy+25),'cyan');
end;
title ('Species kin energy');
legend('protons', 'positrons cold', 'positrons hot', 'electrons right', 'electrons hot', 'electrons cold');
xlabel ('t');
ylabel ('E_k');
grid;


startElmEnergy = 17 + 5*Nspecies;
figure(7);
hold on;
plot(general(1:N,1), general(1:N,startElmEnergy),'red');
plot(general(1:N,1), general(1:N,startElmEnergy+1),'green');
plot(general(1:N,1), general(1:N,startElmEnergy+2),'blue');
plot(general(1:N,1), general(1:N,startElmEnergy+3),'black');
plot(general(1:N,1), general(1:N,startElmEnergy+4),'magenta');
plot(general(1:N,1), general(1:N,startElmEnergy+5),'cyan');
title ('Uelm');
legend('UEx', 'UEy', 'UEz','UBx','UBy', 'UBz');
xlabel ('t');
ylabel ('E');
grid;