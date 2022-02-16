clear;
directory_name = './output/';
file_name = 'ParticleBinning2';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 20;
name = info.Datasets(1).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minE = 0.001;
maxE = 5000.0;
factor = (maxE/minE)^(1.0/(Np-1));
energy(1) = minE;
for i = 2:Np,
    energy(i) = energy(i-1)*factor;
end;

dt = 4000*0.5*0.1;

meanE(1:Ndata) = 0;
for k = 1:Ndata,
    name = info.Datasets(k).Name;
    fp= hdf5read(full_name, name);
    particles = 0;
    for i = 1:Nx,
        for j = 1:Np,
            meanE(k) = meanE(k) + energy(j)*fp(j,i);
            particles = particles+fp(j,i);
        end;
    end;
    meanE(k) = meanE(k)/particles;
end;

time(1:Ndata) = 0;
for i = 1:Ndata,
    time(i) = (i-1)*dt;
end;

figure(1);
plot(time(1:Ndata),meanE(1:Ndata),'red');
title('E');
xlabel('t');
ylabel('E');
grid;
