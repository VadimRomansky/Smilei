clear;
directory_name = './output/';
file_name = 'ParticleBinning0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Color = {'red','blue','green','black','magenta', [1.0,0.6,0]};

N=size(fp,1);
Ns = 5;

Fp(1:N,1:Ns) = 0;

for k=1:Ns,
    Nt = fix(Ndata*k/Ns);
    name = info.Datasets(Nt).Name;
    fp= hdf5read(full_name, name);
    for i=1:N,
        Fp(i,k)=fp(i);
    end;
end;

figure(1);
hold on;
for k = 1:Ns,
    plot(1:N,Fp(1:N,k),'color',Color{k});
end;
title('N(x)');
xlabel('x');
ylabel('N');
grid;