clear;
directory_name = './output/';
file_name = 'ParticleBinning7';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Color = {'red','blue','green','black','magenta'};

Np=size(fp,1);
Nx=size(fp,2);

Ns=5;

Fp(1:Np,1:Ns)=0;

samplingFactor = 20;

startx(1:Ns) = 0; 
startx(1)= 1;
endx(1) = fix(10000/samplingFactor);
for i =2:Ns,
    startx(i) = endx(i-1);
    endx(i) = startx(i)+fix(10000/samplingFactor);
end;

for k=1:Ns,
    for i=1:Np,
        for j=startx(k):endx(k),
            Fp(i,k)=Fp(i,k)+fp(i,j)*i*i;
        end;
    end;
end;

figure(1);
hold on;
for k = 1:Ns,
    plot(1:Np,Fp(1:Np,k),'color',Color{k});
end;
title('F(E) E^2');
xlabel('E/me c^2');
ylabel('F(E) E^2');
grid;