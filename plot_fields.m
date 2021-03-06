clear;
directory_name = './output/';
file_name = 'Fields0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
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
figure(2);
colormap Jet;
[X, Y] = meshgrid(1:Nx, 1:Ny);
surf(X, Y, By);
shading interp;
title ('By');
xlabel ('x');
ylabel ('y');
zlabel ('By');
grid ;

% figure(3);
% colormap Jet;
% [X, Y] = meshgrid(1:Nx, 1:Ny);
% surf(X, Y, Bz);
% shading interp;
% title ('Bz');
% xlabel ('x');
% ylabel ('y');
% zlabel ('Bz');
% grid ;
