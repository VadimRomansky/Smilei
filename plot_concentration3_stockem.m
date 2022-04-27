clear;
directory_name = './output/';
file_name1 = 'ParticleBinning0';
file_name2 = 'ParticleBinning1';
file_name3 = 'ParticleBinning3';
file_number = '.h5';
full_name1 = strcat(directory_name, file_name1, file_number);
full_name2 = strcat(directory_name, file_name2, file_number);
full_name3 = strcat(directory_name, file_name3, file_number);

info = h5info(full_name1);
Ndata = size(info.Datasets,1);
%Ndata = 15;
info = h5info(full_name1);
name1 = info.Datasets(Ndata).Name;
info = h5info(full_name2);
name2 = info.Datasets(Ndata).Name;
info = h5info(full_name3);
name3 = info.Datasets(Ndata).Name;

fp1= hdf5read(full_name1, name1);
fp2 = hdf5read(full_name2, name2);
fp3 = hdf5read(full_name3, name3);
N=size(fp1,1);

dx = 0.1;
factor = 1.0/(dx*dx);
Ny = 20;

xsw(1:Ndata) = 0;
t(1:Ndata) = 0;

for i = 1:Ndata,
    t(i) = (i-1)*0.09*80000;
    name = info.Datasets(i).Name;
    fp = hdf5read(full_name1, name)*factor/Ny;
    for j = N:-1:1,
        if (fp(j) > 2)
            xsw(i) = j*dx;
            break;
        end;
    end;
end;

vsw(1:Ndata) = 0;
vsw1(1:Ndata) = 0;
for i = 2:Ndata,
    vsw(i) = (xsw(i) - xsw(i-1))/(t(i) - t(i-1));
    vsw1(i) = xsw(i)/t(i);
end;

set(0, 'DefaultLineLineWidth', 1);
figure(1);
hold on;
xlim([250000 550000]);
ylim([0 15]);
%plot((1:N)*dx,fp1(1:N)*factor/Ny,'red',(1:N)*dx,fp2(1:N)*factor/Ny,'green',(1:N)*dx,fp3(1:N)*factor/Ny,'blue');
plot((1:N),fp2(1:N)*factor/Ny,'red');
plot((1:N),fp3(1:N)*factor/Ny,'magenta');
plot((1:N),fp1(1:N)*factor/Ny,'blue');
legend('protons','positrons','electrons');
grid;

% figure(2);
% plot(t(1:Ndata), xsw(1:Ndata));
% grid;
% 
% figure(3);
% plot(t(1:Ndata), vsw(1:Ndata));
% grid;
% 
% figure(4);
% plot(t(1:Ndata), vsw1(1:Ndata));
% grid;

output(1:Ndata-1,4) = 0;
for i = 1:Ndata-1,
    output(i,1) = t(i+1);
    output(i,2) = xsw(i+1);
    output(i,3) = vsw(i+1);
    output(i,4) = vsw1(i+1);
end;
dlmwrite('shockwave.dat',output,'delimiter',' ');