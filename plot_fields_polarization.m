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

startx = 4000;
size = Ny-1;
smallBx(1:size,1:Ny-1) = 0;
smallBy(1:size,1:Ny-1) = 0;
smallBz(1:size,1:Ny-1) = 0;
for i = 1:size,
    for j = 1:Ny-1;
        smallBx(i,j) = Bx(j,startx+i);
        smallBy(i,j) = By(j,startx+i);
        smallBz(i,j) = Bz(j,startx+i);
    end;
end;

fourierBx = fft2(smallBx);
fourierBy = fft2(smallBy);
fourierBz = fft2(smallBz);
I(1:size,1:Ny-1)=0;
V(1:size,1:Ny-1)=0;
Hi(1:size,1:Ny-1)=0;
temp(1:size,1:Ny-1) = 0;
scalar(1:size,1:Ny-1) = 0;
divergence(1:size-1,1:Ny-2) = 0;
for i = 1:size-1,
    for j = 1:Ny-2;
        divergence(i,j) = (smallBx(i+1,j) + smallBx(i+1,j+1) - smallBx(i,j+1) - smallBx(i+1,j+1))/0.2 + (smallBy(i,j+1) + smallBy(i+1,j+1) - smallBy(i,j) - smallBy(i+1,j))/0.2;
    end;
end;
for i = 1:size,
    for j = 1:Ny-1;
        scalar(i,j) = ((i-1)/size)*fourierBx(i,j)+((j-1)/(Ny-1))*fourierBy(i,j);
        k2 = 1.0*(i-1)*(i-1)/(size*size) + (j-1)*(j-1)/((Ny-1)*(Ny-1));
        normk = sqrt(k2);
        cosphi = (i-1)/(size*normk);
        sinphi = (j-1)/((Ny-1)*normk);
        if(normk == 0)
            cosphi = 1.0;
            sinphi = 0;
        end;
        fourierBx(i,j) = fourierBx(i,j)*(i-1)*scalar(i,j)/k2;
        fourierBy(i,j) = fourierBx(i,j)*(j-1)*scalar(i,j)/k2;
        fourierBnorm = -fourierBy(i,j)*cosphi + fourierBx(i,j)*sinphi;
        I(i,j) = abs(fourierBx(i,j)*fourierBx(i,j)) + abs(fourierBy(i,j)*fourierBy(i,j)) + abs(fourierBz(i,j)*fourierBz(i,j));
        phi1 = angle(fourierBnorm);
        phi3 = angle(-fourierBy(i,j));
        phi2 = angle(fourierBz(i,j));
        V(i,j) = 2*abs(fourierBz(i,j))*sqrt(abs(fourierBx(i,j)*fourierBx(i,j)) + abs(fourierBy(i,j)*fourierBy(i,j)))*sin(phi1 - phi2);
        temp(i,j) = V(i,j)/I(i,j);
        Hi(i,j) = asin(V(i,j)/I(i,j))/2;
    end;
end;

figure(1);
title ('Hi');
xlabel ('kx');
ylabel ('ky');
zlabel ('Hi');
grid ;
colormap Jet;
[X, Y] = meshgrid(1:Ny-1, 1:size);
surf(X, Y, Hi);
shading interp;
 
