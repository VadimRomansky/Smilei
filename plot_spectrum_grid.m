clear;
directory_name = './output/';
file_name = 'ParticleBinning6';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata)).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);

Np=size(fp1,1);
Nx=size(fp1,2);

minE = 0.1;
maxE = 1000;
factor = (maxE/minE)^(1.0/(Np-1));

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

startx = fix(60000/samplingFactor)+1;
endx = fix(70000/samplingFactor);

for i=1:Np,
    for j=startx:endx,
        Fp1(i)=Fp1(i)+fp1(i,j)/de(i);
        Fp2(i)=Fp2(i)+fp2(i,j)/de(i);
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
plot(energy(1:Np),Fp2(1:Np),'red');
title('F(E)');
xlabel('Ekin/me c^2');
ylabel('F(E)');
name = strcat('approximation gamma = ',num2str(gammap-2));
%legend('Fe', name,'Location','southeast');
grid;