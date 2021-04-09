clear;
directory_name = './output_theta0-90/';
file_name = 'Fields0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);

Color = {'red','blue','green','black','magenta', [1.0,0.6,0]};
%h5disp(full_name);
Ndata = size(info.Groups.Groups,1);
%datasets = info.Groups.Groups(1).Datasets;
%name1x = strcat(info.Groups.Groups(1).Name, '/Bx');
%name1y = strcat(info.Groups.Groups(1).Name, '/By');
%name1z = strcat(info.Groups.Groups(1).Name, '/Bz');
name2x = strcat(info.Groups.Groups(Ndata).Name, '/Bx');
name2y = strcat(info.Groups.Groups(Ndata).Name, '/By');
name2z = strcat(info.Groups.Groups(Ndata).Name, '/Bz');

%Bx1= hdf5read(full_name, name1x);
%By1= hdf5read(full_name, name1y);
%Bz1= hdf5read(full_name, name1z);
Bx= hdf5read(full_name, name2x);
By= hdf5read(full_name, name2y);
Bz= hdf5read(full_name, name2z);

Ny=size(Bx,1);
Nx=size(Bx,2);


Bxa(1:Nx,1:Ndata) = 0;
Bya(1:Nx,1:Ndata) = 0;
Bza(1:Nx,1:Ndata) = 0;
Bnorma(1:Nx,1:Ndata) = 0;

sampling = 4;
beta = 0.75;
diag_every = 20000;
xstart(1:Ndata) = 0;
xend(1:Ndata) = 0;
xstart(1) = 50000;
xend(1) = xstart(1)+100;
for i = 2:Ndata,
    xstart(i) = xstart(1) - fix((i-1)*beta*0.45*diag_every/sampling);
    xend(i) = xstart(i) + 100;
end;

Bnormaa(1:Ndata) = 0;

for k=1:Ndata,
    Nt = k;
    name2x = strcat(info.Groups.Groups(Nt).Name, '/Bx');
    name2y = strcat(info.Groups.Groups(Nt).Name, '/By');
    name2z = strcat(info.Groups.Groups(Nt).Name, '/Bz');
    Bx= hdf5read(full_name, name2x);
    By= hdf5read(full_name, name2y);
    Bz= hdf5read(full_name, name2z);
    for i = 1:Nx,
        for j = 1:Ny,
            Bxa(i,k) = Bxa(i,k) + Bx(j,i)/Ny;
            Bya(i,k) = Bya(i,k) + By(j,i)/Ny;
            Bza(i,k) = Bza(i,k) + Bz(j,i)/Ny;
            Bnorma(i,k) = Bnorma(i,k) + By(j,i)*By(j,i) + Bz(j,i)*Bz(j,i);
        end;
        Bnorma(i,k) = sqrt(Bnorma(i,k)/Ny);
    end;
end;

for k = 1:Ndata,
    for i = xstart(k):xend(k),
        for j = 1:Ny,
            Bnormaa(k) = Bnormaa(k) + By(j,i)*By(j,i) + Bz(j,i)*Bz(j,i);
        end;
    end;
    Bnormaa(k) = sqrt(Bnormaa(k)/(100*Ny));
end;

set(0,'DefaultFigureColormap',feval('jet'));

figure(1);
hold on;
plot(1:Ndata, Bnormaa(1:Ndata));
title ('B_{\perp}');
xlabel ('t');
ylabel ('B_{\perp}');
grid ;
