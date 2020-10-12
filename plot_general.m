clear;
general = importdata('./output/scalars.txt');
N=size(general,1);

figure(1);
plot(1:N, general(1:N,1));
title ('t');
xlabel ('N');
ylabel ('t');
grid;

figure(2);
plot(general(1:N,1), general(1:N,4));
title ('Utot');
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