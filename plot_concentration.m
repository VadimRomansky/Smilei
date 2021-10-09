clear;
directory_name = './output/';
file_name = 'ParticleBinning0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 30;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(2).Name;
name3 = info.Datasets(Ndata).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);
fp3 = hdf5read(full_name, name3);
N=size(fp1,1);

dx = 0.2;
factor = 1.0/(dx*dx);
Ny = 200;

figure(1);
plot(1:N,fp2(1:N)*factor/200,1:N,fp3(1:N)*factor/200);
grid;