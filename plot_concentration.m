clear;
%directory_name = './output_gamma0.3_sigma0.0002_theta30/';
directory_name = './output/';
file_name = 'ParticleBinning1';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 3;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata/2)+1).Name;
name3 = info.Datasets(Ndata).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);
fp3 = hdf5read(full_name, name3);
N=size(fp1,1);

dx = 0.1;
factor = 1.0/(dx*dx);
Ny = 200;

xsw(1:Ndata) = 0;
t(1:Ndata) = 0;

for i = 1:Ndata,
    t(i) = (i-1)*dx*0.45*80000;
    name = info.Datasets(i).Name;
    fp = hdf5read(full_name, name)*factor/Ny;
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

figure(1);
hold on;
%xlim([150000 300000]);
%ylim([0 4]);
%plot((1:N)*dx,fp1(1:N)*factor/Ny,'red',(1:N)*dx,fp2(1:N)*factor/Ny,'green',(1:N)*dx,fp3(1:N)*factor/Ny,'blue');
plot((1:N),fp1(1:N)*factor/Ny,'red');
plot((1:N),fp2(1:N)*factor/Ny,'green');
plot((1:N),fp3(1:N)*factor/Ny,'blue');
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