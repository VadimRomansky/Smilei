clear;
directory_name = './output/';
file_name = 'ParticleBinning10';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 20;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata/2)).Name;
name3 = info.Datasets(fix(Ndata)).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);
fp3 = hdf5read(full_name, name3);

Np=size(fp1,1);
Nx=size(fp1,2);

minV = -1.0;
maxV = 1.0;
dv = (maxV - minV)/Np;

meanV(1:Nx,1:3) = 0;
concentration(1:Nx,1:3) = 0;
samplingfactor = 20;
%for i = 1:Nx,
 for i = 400000/samplingfactor:420000/samplingfactor,
    for j = 1:Np,
        vx = minV + (j-0.5)*dv;
        meanV(i,1) = meanV(i,1) + vx*fp1(j,i);
        meanV(i,2) = meanV(i,2) + vx*fp2(j,i);
        meanV(i,3) = meanV(i,3) + vx*fp3(j,i);
        concentration(i,1) = concentration(i,1) + fp1(j,i);
        concentration(i,2) = concentration(i,2) + fp2(j,i);
        concentration(i,3) = concentration(i,3) + fp3(j,i);
    end;
end;
for i = 1:Nx,
    for j = 1:3,
        meanV(i,j) = meanV(i,j)/concentration(i,j);
    end;
end;

figure(1);
plot(1:Nx,concentration(1:Nx,1),'red',1:Nx,concentration(1:Nx,2),'green',1:Nx,concentration(1:Nx,3),'blue');
title('n');
xlabel('x');
ylabel('n');
grid;

figure(2);
plot(1:Nx,meanV(1:Nx,1),'red',1:Nx,meanV(1:Nx,2),'green',1:Nx,meanV(1:Nx,3),'blue');
title('Vx');
xlabel('x');
ylabel('Vx');
grid;