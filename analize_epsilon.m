clear;
directory_name = './B30mass100/';
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
%Ndata = 7;
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
kB = 1.3806488*10^-16;
c = 2.99792458*10^10;

samplingFactor = 20;
fieldsSamplingFactor = 4;
start = 10000;
fin = 20000;
startFieldx = fix(start/fieldsSamplingFactor)+1;
endFieldx = fix(fin/fieldsSamplingFactor);
startx = fix(start/samplingFactor)+1;
endx = fix(fin/samplingFactor);

Np=size(fp1,1);
Nx=size(fp1,2);
Ny1 = 200;

minElectronE = 0.1;
maxElectronE = 5000;
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

mprotonreal = 1.67*10^-24;
melectronreal = 0.91*10^-27;
Fproton(1:Np)=0;
Felectron(1:Np)=0;
FprotonMaxwell(1:Np) = 0;
FelectronMaxwell(1:Np) = 0;
Te = 4.236*10^10;
Tp = 1.667*10^11;
%norm = 0;
%for i = 1:Np,
%    norm = norm + (fe1(i,10)/(dx*dx*Ny1))/samplingFactor;
%end;

normproton = 0.0;
normelectron = 0.0;
for i = 1:Np,
    for j=startx:endx,
        normproton = normproton + fp2(i,j)*me/mp;
        normelectron = normelectron + fe2(i,j);
        Fproton(i)=Fproton(i)+fp2(i,j)/deProton(i);
        Felectron(i)=Felectron(i)+fe2(i,j)/deElectron(i);
    end;
    
    theta = kB*Tp/(mprotonreal*c*c);
    bes = besselk(2, 1/theta);
    gam = energyProton(i)*me/mp + 1;
    beta = sqrt(1.0 - 1.0/(gam*gam));
    exp1 = exp(-gam/theta);   
    FprotonMaxwell(i) = (1.0/(theta*bes))*exp1*gam*gam*beta;
    
    theta = kB*Te/(melectronreal*c*c);
    bes = besselk(2, 1/theta);
    gam = energyElectron(i) + 1;
    beta = sqrt(1.0 - 1.0/(gam*gam));
    exp1 = exp(-gam/theta);   
    FelectronMaxwell(i) = (1.0/(theta*bes))*exp1*gam*gam*beta;
end;

for i = 1:Np,
    FprotonMaxwell(i) = FprotonMaxwell(i)*normproton;
    FelectronMaxwell(i) = FelectronMaxwell(i)*normelectron;
end;

%figure(1);
%hold on;
%plot(energyProton(1:Np),Fproton(1:Np), energyProton(1:Np), FprotonMaxwell(1:Np));

fluxKineticEnergy = 0;
fluxTotalEnergy = 0;
fluxRestEnergy = 0;
for i=1:Np,
    for j=startx:endx,
        fluxKineticEnergy = fluxKineticEnergy + fp1(i,j)*me*energyProton(i);
        fluxKineticEnergy = fluxKineticEnergy + fe1(i,j)*me*energyElectron(i);
        fluxTotalEnergy = fluxTotalEnergy + fp1(i,j)*(me*energyProton(i) + mp);
        fluxTotalEnergy = fluxTotalEnergy + fe1(i,j)*(me*energyElectron(i) + me);
        fluxRestEnergy = fluxRestEnergy + fp1(i,j)*mp;
        fluxRestEnergy = fluxRestEnergy + fe1(i,j)*me;
    end;
end;

%fluxRestEnergy = fluxTotalEnergy - fluxKineticEnergy;
fluxKineticEnergy2 = fluxRestEnergy*0.375*0.375;

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
protonRestEnergy = 0;
protonTotalEnergy = 0;
protonKineticEnergy = 0;
electronAcceleratedTotalEnergy= 0;
protonAcceleratedTotalEnergy= 0;
electronAcceleratedKineticEnergy= 0;
protonAcceleratedKineticEnergy= 0;

electronAcceleratedLevel = 8;
protonAcceleratedLevel = 2*gamma;

maxwellElectronKineticEnergy = 0;
maxwellElectronTotalEnergy = 0;
maxwellProtonKineticEnergy = 0;
maxwellProtonTotalEnergy = 0;

for i=1:Np,
    
    maxwellElectronKineticEnergy = maxwellElectronKineticEnergy + FelectronMaxwell(i)*me*energyElectron(i)*deElectron(i);
    maxwellElectronTotalEnergy = maxwellElectronTotalEnergy + FelectronMaxwell(i)*(me*energyElectron(i) + me)*deElectron(i);
    maxwellProtonKineticEnergy = maxwellProtonKineticEnergy + FprotonMaxwell(i)*me*energyProton(i)*deProton(i);
    maxwellProtonTotalEnergy = maxwellProtonTotalEnergy + FprotonMaxwell(i)*(me*energyProton(i) + mp)*deProton(i);
    for j=startx:endx,
        protonRestEnergy = protonRestEnergy + fp2(i,j)*mp;
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

massFactor2 = (mprotonreal/massFactor)/melectronreal;

totalEnergy = (magneticEnergy + electronKineticEnergy/massFactor2 + protonKineticEnergy);

epsilonB = magneticEnergy/totalEnergy;
epsilonE = electricEnergy/fluxTotalEnergy;
epsilon_e = (electronKineticEnergy/totalEnergy)/massFactor2;
epsilon_e_accelerated = electronAcceleratedTotalEnergy/fluxTotalEnergy;
epsilon_p = protonKineticEnergy/totalEnergy;
epsilon_p_accelerated = protonAcceleratedTotalEnergy/fluxTotalEnergy;

epsilonBnonrel = magneticEnergy/fluxKineticEnergy;
epsilonEnonrel = electricEnergy/fluxKineticEnergy;
epsilon_enonrel = electronTotalEnergy/fluxKineticEnergy;
epsilon_e_acceleratednonrel = electronAcceleratedTotalEnergy/fluxKineticEnergy;
epsilon_pnonrel = protonTotalEnergy/fluxKineticEnergy;
epsilon_p_acceleratednonrel = protonAcceleratedTotalEnergy/fluxKineticEnergy;

epsilon_e_maxwell_acceleratednonrel = (electronTotalEnergy - maxwellElectronTotalEnergy)/fluxKineticEnergy;
epsilon_e_maxwell_accelerated = (electronTotalEnergy - maxwellElectronTotalEnergy)/fluxTotalEnergy;

epsilon_p_maxwell_acceleratednonrel = (protonTotalEnergy - maxwellProtonTotalEnergy)/fluxKineticEnergy;
epsilon_p_maxwell_accelerated = (protonTotalEnergy - maxwellProtonTotalEnergy)/fluxTotalEnergy;

epsilon_e_kinetic_maxwell_acceleratednonrel = (electronKineticEnergy - maxwellElectronKineticEnergy)/fluxKineticEnergy;
epsilon_e_kinetic_maxwell_accelerated = (electronKineticEnergy - maxwellElectronKineticEnergy)/fluxTotalEnergy;

epsilon_p_kinetic_maxwell_acceleratednonrel = (protonKineticEnergy - maxwellProtonKineticEnergy)/fluxKineticEnergy;
epsilon_p_kinetic_maxwell_accelerated = (protonKineticEnergy - maxwellProtonKineticEnergy)/fluxTotalEnergy;

total_energy = (magneticEnergy + electricEnergy + electronTotalEnergy + protonTotalEnergy);
total_kinetic_energy = (magneticEnergy + electricEnergy + electronKineticEnergy + protonKineticEnergy);

fraction_e = electronTotalEnergy/total_energy;
fraction_e_kinetic = electronKineticEnergy/total_energy;
fraction_e_accelerated = (electronTotalEnergy - maxwellElectronTotalEnergy)/total_energy;

fraction_p = protonTotalEnergy/total_energy;
fraction_p_kinetic = protonKineticEnergy/total_energy;
fraction_p_accelerated = (protonTotalEnergy - maxwellProtonTotalEnergy)/total_energy;

fractionB = magneticEnergy/total_energy;
fractionBrest = magneticEnergy/protonRestEnergy;
fractionE = electricEnergy/total_energy;

fraction_e_accelerated_kinetic = (electronKineticEnergy - maxwellElectronKineticEnergy)/total_kinetic_energy;

fraction_e_accelerated_kinetic_scaled = (electronKineticEnergy - maxwellElectronKineticEnergy)/total_kinetic_energy/sqrt(18);

fractionB_kinetic = magneticEnergy/total_kinetic_energy;





startPowerP = 95;
endPowerP = 115;

startPowerE = 140;
endPowerE = 150;

startPower = startPowerE;
endPower = endPowerE;

Fea(1:Np) = 0;

Fea(startPower) = Felectron(startPower);
Fea(endPower) = Felectron(endPower);

gammap = log(Fea(startPower)/Fea(endPower))/log((energyElectron(startPower)+1)/(energyElectron(endPower)+1));


ap = exp(log(Fea(startPower)) - gammap*log((energyElectron(startPower)+1)));

for i = 1:startPower-1,
    Fea(i) = Felectron(i);
end;
for i = startPower:Np,
    Fea(i) = ap*((energyElectron(i)+1)^gammap);
end;

extrapolatedElectronTotalEnergy = 0;
extrapolatedElectronKineticEnergy = 0;
for i = 1:Np,
    extrapolatedElectronKineticEnergy = extrapolatedElectronKineticEnergy + Fea(i)*(me*energyElectron(i))*deElectron(i);
    extrapolatedElectronTotalEnergy = extrapolatedElectronTotalEnergy + Fea(i)*(me*energyElectron(i) + me)*deElectron(i);
end;

epsilon_e_accelerated_extrapolated = (extrapolatedElectronTotalEnergy - maxwellElectronTotalEnergy)/fluxTotalEnergy;
epsilon_e_accelerated_extrapolatednonrel = (extrapolatedElectronTotalEnergy - maxwellElectronTotalEnergy)/fluxKineticEnergy;

epsilon_e_chev = (electronKineticEnergy - maxwellElectronKineticEnergy)/fluxKineticEnergy2/sqrt(18);
epsilon_B_chev = magneticEnergy/fluxKineticEnergy2;