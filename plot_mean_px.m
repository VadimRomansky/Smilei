clear;
directory_name = './output/';
file_name = 'ParticleBinning10';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 20;
name = info.Datasets(1).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minP = -100.0;
maxP = 100.0;
dp = (maxP - minP)/Np;

dt = 4000*0.5*0.1;

meanP(1:Ndata) = 0;
for k = 1:Ndata,
    name = info.Datasets(k).Name;
    fp= hdf5read(full_name, name);
    particles = 0;
    for i = 1:Nx,
        for j = 1:Np,
            vx = minP + (j-0.5)*dp;
            meanP(k) = meanP(k) + vx*fp(j,i);
            particles = particles+fp(j,i);
        end;
    end;
    meanP(k) = meanP(k)/particles;
end;

time(1:Ndata) = 0;
for i = 1:Ndata,
    time(i) = (i-1)*dt;
end;

figure(1);
plot(time(1:Ndata),meanP(1:Ndata),'red');
title('Px');
xlabel('t');
ylabel('Px');
grid;
