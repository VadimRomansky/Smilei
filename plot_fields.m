clear;
directory_name = './output/';
file_name = 'Fields0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
%h5disp(full_name);
Ndata = size(info.Groups.Groups,1);
%Ndata = 1;
%datasets = info.Groups.Groups(1).Datasets;
name1x = strcat(info.Groups.Groups(1).Name, '/Bx');
name1y = strcat(info.Groups.Groups(1).Name, '/By');
name1z = strcat(info.Groups.Groups(1).Name, '/Bz');
name2x = strcat(info.Groups.Groups(Ndata).Name, '/Bx');
name2y = strcat(info.Groups.Groups(Ndata).Name, '/By');
name2z = strcat(info.Groups.Groups(Ndata).Name, '/Bz');
name3x = strcat(info.Groups.Groups(Ndata).Name, '/Ex');
name3y = strcat(info.Groups.Groups(Ndata).Name, '/Ey');
name3z = strcat(info.Groups.Groups(Ndata).Name, '/Ez');

Bx1= hdf5read(full_name, name1x);
By1= hdf5read(full_name, name1y);
Bz1= hdf5read(full_name, name1z);

%B0 = sqrt(Bx1(10,10)*Bx1(10,10) + By1(10,10)*By1(10,10) + Bz1(10,10)*Bz1(10,10));
Bx= hdf5read(full_name, name2x);
By= hdf5read(full_name, name2y);
Bz= hdf5read(full_name, name2z);
Ex= hdf5read(full_name, name3x);
Ey= hdf5read(full_name, name3y);
Ez= hdf5read(full_name, name3z);

Ny=size(Bx,1);
Nx=size(Bx,2);

set(0,'DefaultFigureColormap',feval('jet'));

% figure(1);
% colormap Jet;
% [X, Y] = meshgrid(1:Nx, 1:Ny);
% surf(X, Y, Bx);
% shading interp;
% title ('Bx');
% xlabel ('x');
% ylabel ('y');
% zlabel ('Bx');
% grid ;
% 
% sampling = 20;
% startx = fix(Nx/10);
% endx = fix(Nx/5);
% tempB(1:Ny, 1:(endx-startx + 1)) = 0;
% for i = startx:endx,
%     for j = 1:Ny,
%         tempB(j,i-startx+1) = sqrt(Bx(j,i)*Bx(j,i) + By(j,i)*By(j,i) + Bz(j,i)*Bz(j,i))/B0;
%     end;
% end;
% figure(2);
% colormap Jet;
% [X, Y] = meshgrid((startx:endx)*0.2*sampling, (1:Ny)*0.4*sampling);
% surf(X, Y, tempB);
% shading interp;
% colorbar;
% title ('B');
% xlabel ('x \omega_e /c');
% ylabel ('y \omega_e /c');
% zlabel ('B/B_0');
% grid ;

figure(3);
colormap Jet;
[X, Y] = meshgrid(1:Nx, 1:Ny);
surf(X, Y, Bx);
shading interp;
title ('Bx');
xlabel ('x');
ylabel ('y');
zlabel ('Bx');
grid ;

figure(4);
colormap Jet;
[X, Y] = meshgrid(1:Nx, 1:Ny);
surf(X, Y, By);
shading interp;
title ('By');
xlabel ('x');
ylabel ('y');
zlabel ('By');
grid ;

figure(5);
colormap Jet;
[X, Y] = meshgrid(1:Nx, 1:Ny);
surf(X, Y, Bz);
shading interp;
title ('Bz');
xlabel ('x');
ylabel ('y');
zlabel ('Bz');
grid ;

figure(6);
colormap Jet;
[X, Y] = meshgrid(1:Nx, 1:Ny);
surf(X, Y, Ex);
shading interp;
title ('Ex');
xlabel ('x');
ylabel ('y');
zlabel ('Ex');
grid ;

figure(7);
colormap Jet;
[X, Y] = meshgrid(1:Nx, 1:Ny);
surf(X, Y, Ey);
shading interp;
title ('Ey');
xlabel ('x');
ylabel ('y');
zlabel ('Ey');
grid ;

figure(8);
colormap Jet;
[X, Y] = meshgrid(1:Nx, 1:Ny);
surf(X, Y, Ez);
shading interp;
title ('Ez');
xlabel ('x');
ylabel ('y');
zlabel ('Ez');
grid ;
