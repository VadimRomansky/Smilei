clear;
directory_name = './output/';
file_electron_name = 'ParticleBinning6';
file_proton_name = 'ParticleBinning7';
file_number = '';
file_extension = '.h5';
field_file_name = 'Fields0';
full_electron_name = strcat(directory_name, file_electron_name, file_number, file_extension);
full_proton_name = strcat(directory_name, file_proton_name, file_number, file_extension);
field_name = strcat(directory_name, field_file_name, file_number, file_extension);
info = h5info(full_electron_name);
Ndata = size(info.Datasets,1);
Ndata = 7;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata)).Name;
fp1= hdf5read(full_proton_name, name1);
fp2 = hdf5read(full_proton_name, name2);
fe1= hdf5read(full_electron_name, name1);
fe2 = hdf5read(full_electron_name, name2);

me = 1.0;
massFactor = 100;
mp = me*massFactor;

v = 0.1;
gamma = 1.0/sqrt(1.0 - v*v);

dx = 0.2;

samplingFactor = 20;
fieldsSamplingFactor = 4;
start = 5000;
fin = 10000;
startFieldx = fix(start/fieldsSamplingFactor)+1;
endFieldx = fix(fin/fieldsSamplingFactor);
startx = fix(start/samplingFactor)+1;
endx = fix(fin/samplingFactor);

Np=size(fp1,1);
Nx=size(fp1,2);
Ny1 = 200;

minElectronE = 0.001;
maxElectronE = 1000;
factorElectron = (maxElectronE/minElectronE)^(1.0/(Np-1));

energyElectron(1:Np) = 0;
deElectron(1:Np) = 0;
energyElectron(1) = minElectronE;
for i = 2:Np,
    energyElectron(i) = energyElectron(i-1)*factorElectron;
end;
deElectron(1) = energyElectron(2) - energyElectron(1);
for i = 2:Np,
    deElectron(i) = energyElectron(i) - energyElectron(i-1);
end;

minProtonE = 0.1;
maxProtonE = 5000;
factorProton = (maxProtonE/minProtonE)^(1.0/(Np-1));

energyProton(1:Np) = 0;
deProton(1:Np) = 0;
energyProton(1) = minProtonE;
for i = 2:Np,
    energyProton(i) = energyProton(i-1)*factorProton;
end;
deProton(1) = energyProton(2) - energyProton(1);
for i = 2:Np,
    deProton(i) = energyProton(i) - energyProton(i-1);
end;

Fp1(1:Np)=0;
Fp2(1:Np)=0;

norm = 0;
for i = 1:Np,
    norm = norm + (fe1(i,10)/(dx*dx*Ny1))/samplingFactor;
end;

fluxKineticEnergy = 0;
fluxTotalEnergy = 0;

for i=1:Np,
    for j=startx:endx,
        fluxKineticEnergy = fluxKineticEnergy + fp1(i,j)*me*energyProton(i);
        fluxKineticEnergy = fluxKineticEnergy + fe1(i,j)*me*energyElectron(i);
        fluxTotalEnergy = fluxTotalEnergy + fp1(i,j)*(me*energyProton(i) + mp);
        fluxTotalEnergy = fluxTotalEnergy + fe1(i,j)*(me*energyElectron(i) + me);
    end;
end;

info = h5info(field_name);
%Ndata = size(info.Groups.Groups,1);
name1x = strcat(info.Groups.Groups(1).Name, '/Bx');
name1y = strcat(info.Groups.Groups(1).Name, '/By');
name1z = strcat(info.Groups.Groups(1).Name, '/Bz');
name2x = strcat(info.Groups.Groups(Ndata).Name, '/Bx');
name2y = strcat(info.Groups.Groups(Ndata).Name, '/By');
name2z = strcat(info.Groups.Groups(Ndata).Name, '/Bz');
name1ex = strcat(info.Groups.Groups(1).Name, '/Ex');
name1ey = strcat(info.Groups.Groups(1).Name, '/Ey');
name1ez = strcat(info.Groups.Groups(1).Name, '/Ez');
name2ex = strcat(info.Groups.Groups(Ndata).Name, '/Ex');
name2ey = strcat(info.Groups.Groups(Ndata).Name, '/Ey');
name2ez = strcat(info.Groups.Groups(Ndata).Name, '/Ez');

Bx1= hdf5read(field_name, name1x);
By1= hdf5read(field_name, name1y);
Bz1= hdf5read(field_name, name1z);
Bx= hdf5read(field_name, name2x);
By= hdf5read(field_name, name2y);
Bz= hdf5read(field_name, name2z);

Ex1= hdf5read(field_name, name1ex);
Ey1= hdf5read(field_name, name1ey);
Ez1= hdf5read(field_name, name1ez);
Ex= hdf5read(field_name, name2ex);
Ey= hdf5read(field_name, name2ey);
Ez= hdf5read(field_name, name2ez);

Ny = size(Bx1,1);

magneticInitialEnergy = 0;
magneticEnergy = 0;
for i = startFieldx:endFieldx,
    for j = 1:Ny,
        magneticInitialEnergy = magneticInitialEnergy + (Bx1(j,i)*Bx1(j,i) + By1(j,i)*By1(j,i) + Bz1(j,i)*Bz1(j,i))*(0.5*dx*dx);
        magneticEnergy = magneticEnergy + (Bx(j,i)*Bx(j,i) + By(j,i)*By(j,i) + Bz(j,i)*Bz(j,i))*(0.5*dx*dx);
    end;
end;
magneticInitialEnergy = magneticInitialEnergy*fieldsSamplingFactor*fieldsSamplingFactor;
magneticEnergy = magneticEnergy*fieldsSamplingFactor*fieldsSamplingFactor;
electricInitialEnergy = 0;
electricEnergy = 0;
for i = startFieldx:endFieldx,
    for j = 1:Ny,
        electricInitialEnergy = electricInitialEnergy + (Ex1(j,i)*Ex1(j,i) + Ey1(j,i)*Ey1(j,i) + Ez1(j,i)*Ez1(j,i))*(0.5*dx*dx);
        electricEnergy = electricEnergy + (Ex(j,i)*Ex(j,i) + Ey(j,i)*Ey(j,i) + Ez(j,i)*Ez(j,i))*(0.5*dx*dx);
    end;
end;
electricInitialEnergy = electricInitialEnergy*fieldsSamplingFactor*fieldsSamplingFactor;
electricEnergy = electricEnergy*fieldsSamplingFactor*fieldsSamplingFactor;

initialSigma = 2*magneticInitialEnergy/fluxTotalEnergy;

electronTotalEnergy = 0;
electronKineticEnergy = 0;
protonTotalEnergy = 0;
protonKineticEnergy = 0;
electronAcceleratedTotalEnergy= 0;
protonAcceleratedTotalEnergy= 0;
electronAcceleratedKineticEnergy= 0;
protonAcceleratedKineticEnergy= 0;

electronAcceleratedLevel = 2*gamma;
protonAcceleratedLevel = 2*gamma;

for i=1:Np,
    for j=startx:endx,
        protonKineticEnergy = protonKineticEnergy + fp2(i,j)*me*energyProton(i);
        electronKineticEnergy = electronKineticEnergy + fe2(i,j)*me*energyElectron(i);
        protonTotalEnergy = protonTotalEnergy + fp2(i,j)*(me*energyProton(i) + mp);
        electronTotalEnergy = electronTotalEnergy + fe2(i,j)*(me*energyElectron(i) + me);
        electronGamma = 1.0 + energyElectron(i);
        protonGamma = 1.0 + energyProton(i)*me/mp;
        if(electronGamma > electronAcceleratedLevel)
            electronAcceleratedKineticEnergy = electronAcceleratedKineticEnergy + fe2(i,j)*me*energyElectron(i);
            electronAcceleratedTotalEnergy = electronAcceleratedTotalEnergy + fe2(i,j)*(me*energyElectron(i) + me);
        end;
        if(protonGamma > protonAcceleratedLevel)
            protonAcceleratedKineticEnergy = protonAcceleratedKineticEnergy + fp2(i,j)*me*energyElectron(i);
            protonAcceleratedTotalEnergy = protonAcceleratedTotalEnergy + fp2(i,j)*(me*energyElectron(i) + mp);
        end;
    end;
end;

epsilonB = magneticEnergy/fluxTotalEnergy;
epsilonE = electricEnergy/fluxTotalEnergy;
epsilon_e = electronTotalEnergy/fluxTotalEnergy;
epsilon_e_accelerated = electronAcceleratedTotalEnergy/fluxTotalEnergy;
epsilon_p = protonTotalEnergy/fluxTotalEnergy;
epsilon_p_accelerated = protonAcceleratedTotalEnergy/fluxTotalEnergy;

epsilonBnonrel = magneticEnergy/fluxKineticEnergy;
epsilonEnonrel = electricEnergy/fluxKineticEnergy;
epsilon_enonrel = electronTotalEnergy/fluxKineticEnergy;
epsilon_e_acceleratednonrel = electronAcceleratedTotalEnergy/fluxKineticEnergy;
epsilon_pnonrel = protonTotalEnergy/fluxKineticEnergy;
epsilon_p_acceleratednonrel = protonAcceleratedTotalEnergy/fluxKineticEnergy;
