clear;
directory_name = './output1/';
file_name = 'ParticleBinning63';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name1 = info.Datasets(1).Name;
name2 = info.Datasets(Ndata).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);

Np=size(fp1,1);
Nx=size(fp1,2);

minE = 0.1;
maxE = 1000;
factor = (maxE/minE)^(1.0/(Np-1));

energy(1:Np) = 0;
energy(1) = minE;
for i = 2:Np,
    energy(i) = energy(i-1)*factor;
end;

Fp1(1:Np)=0;
Fp2(1:Np)=0;

samplingFactor = 20;

startx = fix(20000/samplingFactor)+1;
endx = fix(30000/samplingFactor);

for i=1:Np,
    for j=startx:endx,
        Fp1(i)=Fp1(i)+fp1(i,j)*(energy(i) + 1)*(energy(i) + 1);
        Fp2(i)=Fp2(i)+fp2(i,j)*(energy(i) + 1)*(energy(i) + 1);
    end;
end;

figure(1);
plot(energy(1:Np),Fp1(1:Np),'red',energy(1:Np),Fp2(1:Np),'green');
title('F(E) E^2');
xlabel('E/me c^2');
ylabel('F(E) E^2');
grid;