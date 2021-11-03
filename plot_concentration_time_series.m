clear;
directory_name = './output_theta0-90_gamma0.3sigma0.0002/';
file_name = 'ParticleBinning01';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 11;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Color = {'red','blue','green','black','magenta', [1.0,0.6,0]};

N=size(fp,1);
Ns = 5;
Ny = 200;
dx = 0.2;
x(1:N) = (1:N)*dx;

Fp(1:N,1:Ns) = 0;

for k=1:Ns,
    Nt = fix(1 + (Ndata-1)*k/Ns);
    name = info.Datasets(Nt).Name;
    fp= hdf5read(full_name, name)/(Ny*dx*dx);
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

output(1:N/5,1:(Ns+1)) = 0;
for i = 1:N/5,
    output(i,1) = x(i);
    for j = 1:Ns,
        output(i,j+1) = Fp(i,j);
    end;
end;

dlmwrite('concentrations.dat', output, 'delimiter',' ');