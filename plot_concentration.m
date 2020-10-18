clear;
directory_name = './output1/';
file_name = 'ParticleBinning01';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name1 = info.Datasets(1).Name;
name2 = info.Datasets(Ndata).Name;
name3 = info.Datasets(25).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);
fp3 = hdf5read(full_name, name3);
N=size(fp1,1);

figure(1);
plot(1:N,fp1(1:N),1:N,fp3(1:N),1:N,fp2(1:N));
grid;