clear;
general = importdata('./output/scalars.txt');
N=size(general,1);

gamma = 1.5;
v = 0.75;
sigma = 0.004;
E1 = general(1,9)/general(1,6);
Eb = general(1,9)/(1 + v*v);
Etot = general(1,6)*gamma/(gamma - 1);
E2 = Eb/Etot;

figure(1);
plot(general(1:N,1), general(1:N,4));
title ('Utot');
xlabel ('t');
ylabel ('E');
grid;

figure(2);
plot(general(1:N,1), general(1:N,5));
title ('Uexp');
xlabel ('t');
ylabel ('E');
grid;

figure(3);
plot(general(1:N,1), general(1:N,6));
title ('Ukin');
xlabel ('t');
ylabel ('E');
grid;

figure(4);
plot(general(1:N,1), general(1:N,9));
title ('Uelm');
xlabel ('t');
ylabel ('E');
grid;

figure(5);
plot(general(1:N,1), general(1:N,34));
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