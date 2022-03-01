clear;
directory_name = './output/';
file_name = 'ParticleBinning9';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 20;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);
samplingfactor = 20;

minV = -1.0;
maxV = 1.0;
dv = (maxV - minV)/Np;

dt = 4000*0.5*0.1;

meanV(1:Ndata) = 0;
for k = 1:Ndata,
    name = info.Datasets(k).Name;
    fp= hdf5read(full_name, name);
    particles = 0;
    %for i = 40000/samplingfactor:50000/samplingfactor,
    for i = 1:Nx,
        for j = 1:Np,
            vx = minV + (j-0.5)*dv;
            meanV(k) = meanV(k) + vx*fp(j,i);
            particles = particles+fp(j,i);
        end;
    end;
    meanV(k) = meanV(k)/particles;
end;

time(1:Ndata) = 0;
for i = 1:Ndata,
    time(i) = (i-1)*dt;
end;

figure(1);
plot(time(1:Ndata),meanV(1:Ndata),'red');
title('V');
xlabel('t');
ylabel('V');
grid;
