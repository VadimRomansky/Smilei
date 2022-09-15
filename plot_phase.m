clear;
directory_name = './output/';
file_name = 'ParticleBinning2';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 1;
name3 = info.Datasets(Ndata).Name;

fp = hdf5read(full_name, name3);

Np=size(fp,1);
Nx=size(fp,2);

minEe = 0.001;
maxEe = 5000;
minEp = 0.1;
maxEp = 5000;
minE = minEe;
maxE = maxEe;
factor = (maxE/minE)^(1.0/(Np-1));

me = 1;
mp = 100;
m = me;

energy(1:Np) = 0;
de(1:Np) = 0;
energy(1) = minE;
for i = 2:Np,
    energy(i) = energy(i-1)*factor;
end;
de(1) = energy(2) - energy(1);
for i = 2:Np,
    de(i) = energy(i) - energy(i-1);
end;

Fp(1:Np)=0;

samplingFactor = 20;

startx = fix(10000/samplingFactor)+1;
endx = fix(60000/samplingFactor);

for i=1:Np,
    for j=startx:endx,
        Fp(i)=Fp(i)+fp(i,j)/de(i);
    end;
end;


set(0, 'DefaultLineLineWidth', 2);
figure(1);
[X, Y] = meshgrid(1:Nx, energy);
h = pcolor(X, Y, log(fp));
set(gca, 'YScale', 'log');
%set(gca, 'ZScale', 'log');
set(h, 'EdgeColor', 'none')
shading interp;
colormap('jet');
colorbar;
title('F(E,x)');
xlabel('Ekin/me c^2');
ylabel('F(E)');
%legend('Fe', name,'Location','southeast');
%legend('t=1000','t=500','t=0')
grid;