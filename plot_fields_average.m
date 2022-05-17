clear;
directory_name = './output/';
file_name = 'Fields0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);

Color = {'red','blue','green','black','magenta', [1.0,0.6,0]};
%h5disp(full_name);
Ndata = size(info.Groups.Groups,1);
%Ndata = 1;
%datasets = info.Groups.Groups(1).Datasets;
%name1x = strcat(info.Groups.Groups(1).Name, '/Bx');
%name1y = strcat(info.Groups.Groups(1).Name, '/By');
%name1z = strcat(info.Groups.Groups(1).Name, '/Bz');
name2x = strcat(info.Groups.Groups(Ndata).Name, '/Bx');
name2y = strcat(info.Groups.Groups(Ndata).Name, '/By');
name2z = strcat(info.Groups.Groups(Ndata).Name, '/Bz');
name3x = strcat(info.Groups.Groups(Ndata).Name, '/Ex');
name3y = strcat(info.Groups.Groups(Ndata).Name, '/Ey');
name3z = strcat(info.Groups.Groups(Ndata).Name, '/Ez');


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
Exa(1:Nx,1:Ns) = 0;
Eya(1:Nx,1:Ns) = 0;
Eza(1:Nx,1:Ns) = 0;
Bnorma(1:Nx,1:Ns) = 0;

for k=1:Ns,
    Nt = fix((Ndata-1)*k/Ns)+1;
    name2x = strcat(info.Groups.Groups(Nt).Name, '/Bx');
    name2y = strcat(info.Groups.Groups(Nt).Name, '/By');
    name2z = strcat(info.Groups.Groups(Nt).Name, '/Bz');
    name3x = strcat(info.Groups.Groups(Nt).Name, '/Ex');
    name3y = strcat(info.Groups.Groups(Nt).Name, '/Ey');
    name3z = strcat(info.Groups.Groups(Nt).Name, '/Ez');
    Bx= hdf5read(full_name, name2x);
    By= hdf5read(full_name, name2y);
    Bz= hdf5read(full_name, name2z);
    Ex= hdf5read(full_name, name3x);
    Ey= hdf5read(full_name, name3y);
    Ez= hdf5read(full_name, name3z);
    for i = 1:Nx,
        for j = 1:Ny,
            Bxa(i,k) = Bxa(i,k) + Bx(j,i)/Ny;
            Bya(i,k) = Bya(i,k) + By(j,i)/Ny;
            Bza(i,k) = Bza(i,k) + Bz(j,i)/Ny;
            Exa(i,k) = Exa(i,k) + Ex(j,i)/Ny;
            Eya(i,k) = Eya(i,k) + Ey(j,i)/Ny;
            Eza(i,k) = Eza(i,k) + Ez(j,i)/Ny;

%             Bxa(i,k) = Bxa(i,k) + Bx(i,j)/Ny;
%             Bya(i,k) = Bya(i,k) + By(i,j)/Ny;
%             Bza(i,k) = Bza(i,k) + Bz(i,j)/Ny;
%             Exa(i,k) = Exa(i,k) + Ex(i,j)/Ny;
%             Eya(i,k) = Eya(i,k) + Ey(i,j)/Ny;
%             Eza(i,k) = Eza(i,k) + Ez(i,j)/Ny;
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
    plot((1:Nx), Bnorma(1:Nx,k),'color',Color{k});
end;
legend('700 {\omega}_{pi}^{-1}', '1400 {\omega}_{pi}^{-1}', '2100 {\omega}_{pi}^{-1}', '2800 {\omega}_{pi}^{-1}', '3500 {\omega}_{pi}^{-1}','Location','northwest');
%legend('600 wpi^{-1}', '1200 wpi^{-1}', '1800 wpi^{-1}', '2400 wpi^{-1}', '3000 wpi^{-1}', '3600 wpi^{-1}','Location','northwest');
title ('B_{\perp}');
xlabel ('x {\omega}_{pi} /c');
ylabel ('B_{\perp}');
grid ;

figure(5);
hold on;
for k = 1:Ns,
    plot((1:Nx), Exa(1:Nx,k),'color',Color{k});
end;
title ('Ex');
xlabel ('x {\omega}_{pi} /c');
ylabel ('Ex');
grid ;

figure(6);
hold on;
for k = 1:Ns,
    plot((1:Nx), Eya(1:Nx,k),'color',Color{k});
end;
title ('Ey');
xlabel ('x {\omega}_{pi} /c');
ylabel ('Ey');
grid ;

figure(7);
hold on;
for k = 1:Ns,
    plot((1:Nx), Eza(1:Nx,k),'color',Color{k});
end;
title ('Ez');
xlabel ('x {\omega}_{pi} /c');
ylabel ('Ez');
grid ;
