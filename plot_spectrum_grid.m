clear;
directory_name = './output/';
file_name = 'ParticleBinning6';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 30;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata)).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);

Np=size(fp1,1);
Nx=size(fp1,2);

minEe = 0.001;
maxEe = 1000;
minEp = 0.1;
maxEp = 5000;
minE = minEe;
maxE = maxEe;
factor = (maxE/minE)^(1.0/(Np-1));

me = 1;
mp = 100;
m = me;

energy(1:Np) = 0;
de(1:Np) = 0;
energy(1) = minE;
for i = 2:Np,
    energy(i) = energy(i-1)*factor;
end;
de(1) = energy(2) - energy(1);
for i = 2:Np,
    de(i) = energy(i) - energy(i-1);
end;

Fp1(1:Np)=0;
Fp2(1:Np)=0;

samplingFactor = 20;

startx = fix(7000/samplingFactor)+1;
endx = fix(12000/samplingFactor);

for i=1:Np,
    for j=startx:endx,
        Fp1(i)=Fp1(i)+fp1(i,j)*(energy(i) + m)*(energy(i) + m)/de(i);
        Fp2(i)=Fp2(i)+fp2(i,j)*(energy(i) + m)*(energy(i) + m)/de(i);
    end;
end;


startPowerP = 160;
endPowerP = 170;

Fpa(1:Np) = 0;

Fpa(startPowerP) = Fp2(startPowerP);
Fpa(endPowerP) = Fp2(endPowerP);

gammap = log(Fpa(startPowerP)/Fpa(endPowerP))/log(energy(startPowerP)/energy(endPowerP));


ap = exp(log(Fpa(startPowerP)) - gammap*log(energy(startPowerP)));

for i = startPowerP:endPowerP,
    Fpa(i) = ap*(energy(i)^gammap);
end;


figure(1);
%plot(energy(1:Np),Fp2(1:Np),'red', energy(startPowerP:endPowerP), Fpa(startPowerP:endPowerP),'blue');
loglog(energy(1:Np),Fp2(1:Np),'red');
title('F(E)');
xlabel('Ekin/me c^2');
ylabel('F(E)*E^2');
name = strcat('approximation gamma = ',num2str(gammap-2));
%legend('Fe', name,'Location','southeast');
grid;

dlmwrite('Ep0.dat',energy,'delimiter',' ');
dlmwrite('Fp0.dat',Fp2,'delimiter',' ');