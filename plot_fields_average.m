clear;
directory_name = './output/';
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

Ns = 5;

Bxa(1:Nx,1:Ns) = 0;
Bya(1:Nx,1:Ns) = 0;
Bza(1:Nx,1:Ns) = 0;
Bnorma(1:Nx,1:Ns) = 0;

for k=1:Ns,
    Nt = fix(Ndata*k/Ns);
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

set(0,'DefaultFigureColormap',feval('jet'));

figure(1);
hold on;
for k = 1:Ns,
    plot((1:Nx), Bxa(1:Nx,k),'color',Color{k});
end;
title ('Bx');
xlabel ('x {\omega}_{pi} /c');
ylabel ('Bx');
grid ;

figure(2);
hold on;
for k = 1:Ns,
    plot((1:Nx), Bya(1:Nx,k),'color',Color{k});
end;
title ('By');
xlabel ('x {\omega}_{pi} /c');
ylabel ('By');
grid ;

figure(3);
hold on;
for k = 1:Ns,
    plot((1:Nx), Bza(1:Nx,k),'color',Color{k});
end;
title ('Bz');
xlabel ('x {\omega}_{pi} /c');
ylabel ('Bz');
grid ;

figure(4);
hold on;
for k = 1:Ns,
    plot((1:Nx), smooth(Bnorma(1:Nx,k),100),'color',Color{k});
end;
legend('700 {\omega}_{pi}^{-1}', '1400 {\omega}_{pi}^{-1}', '2100 {\omega}_{pi}^{-1}', '2800 {\omega}_{pi}^{-1}', '3500 {\omega}_{pi}^{-1}','Location','northwest');
%legend('600 wpi^{-1}', '1200 wpi^{-1}', '1800 wpi^{-1}', '2400 wpi^{-1}', '3000 wpi^{-1}', '3600 wpi^{-1}','Location','northwest');
title ('B_{\perp}');
xlabel ('x {\omega}_{pi} /c');
ylabel ('B_{\perp}');
grid ;
