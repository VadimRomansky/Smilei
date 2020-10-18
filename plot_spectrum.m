clear;
directory_name = './output1/';
file_name = 'ParticleBinning63';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name1 = info.Datasets(1).Name;
name2 = info.Datasets(Ndata).Name;
name3 = info.Datasets(fix(Ndata/4)).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);
fp3 = hdf5read(full_name, name3);

N=size(fp1,1);

minE = 0.1;
maxE = 1000;
factor = (maxE/minE)^(1.0/(N-1));

energy(1:N) = 0;
energy(1) = minE;
for i = 2:N,
    energy(i) = energy(i-1)*factor;
end;

figure(1);
plot(energy(1:N),fp1(1:N),'red',energy(1:N), fp3(1:N),'green',energy(1:N),fp2(1:N),'blue');
grid;