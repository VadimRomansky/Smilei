clear;
%directory_name = './output_theta0-90_gamma1.5_sigma0.004/';
directory_name = './output/';
file_name1 = 'ParticleBinning6';
file_name2 = 'ParticleBinning7';
file_number = '.h5';
full_name1 = strcat(directory_name, file_name1, file_number);
full_name2 = strcat(directory_name, file_name2, file_number);
info = h5info(full_name1);
Ndata = size(info.Datasets,1);
Ndata = 4;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata/2)+1).Name;
name3 = info.Datasets(Ndata).Name;

info = h5info(full_name1);
fe1= hdf5read(full_name1, name1);
fe2 = hdf5read(full_name1, name2);
fe3 = hdf5read(full_name1, name3);

name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata/2)+1).Name;
name3 = info.Datasets(Ndata).Name;

fp1= hdf5read(full_name2, name1);
fp2 = hdf5read(full_name2, name2);
fp3 = hdf5read(full_name2, name3);

Np=size(fe1,1);
Nx=size(fe1,2);

minEe = 0.1;
maxEe = 10000;
minEp = 0.1;
maxEp = 50000;
minE1 = minEe;
maxE1 = maxEe;
factor1 = (maxE1/minE1)^(1.0/(Np-1));
minE2 = minEp;
maxE2 = maxEp;
factor2 = (maxE2/minE2)^(1.0/(Np-1));

me = 1;
mp = 100;
m1 = me;
m2 = mp;

energy1(1:Np) = 0;
de1(1:Np) = 0;
energy1(1) = minE1;

energy2(1:Np) = 0;
de2(1:Np) = 0;
energy2(1) = minE2;
for i = 2:Np,
    energy1(i) = energy1(i-1)*factor1;
    energy2(i) = energy2(i-1)*factor2;
end;
de1(1) = energy1(2) - energy1(1);
de2(1) = energy2(2) - energy2(1);
for i = 2:Np,
    de1(i) = energy1(i) - energy1(i-1);
    de2(i) = energy2(i) - energy2(i-1);
end;

Fe1(1:Np)=0;
Fe2(1:Np)=0;
Fe3(1:Np)=0;

Fp1(1:Np)=0;
Fp2(1:Np)=0;
Fp3(1:Np)=0;

samplingFactor = 20;

startx = fix(10000/samplingFactor)+1;
endx = fix(80000/samplingFactor);

for i=1:Np,
    for j=startx:endx,
        Fe1(i)=Fe1(i)+fe1(i,j)*energy1(i)*energy1(i)/de1(i);
        Fe2(i)=Fe2(i)+fe2(i,j)*energy1(i)*energy1(i)/de1(i);
        Fe3(i)=Fe3(i)+fe3(i,j)*energy1(i)*energy1(i)/de1(i);
        
        Fp1(i)=Fp1(i)+fp1(i,j)*energy2(i)*energy2(i)/de2(i);
        Fp2(i)=Fp2(i)+fp2(i,j)*energy2(i)*energy2(i)/de2(i);
        Fp3(i)=Fp3(i)+fp3(i,j)*energy2(i)*energy2(i)/de2(i);
    end;
end;

set(0, 'DefaultLineLineWidth', 2);
figure(1);
hold on;
%plot(energy(1:Np),Fp2(1:Np),'red', energy(startPowerP:endPowerP), Fpa(startPowerP:endPowerP),'blue');
%loglog(energy(1:Np),Fp2(1:Np),'red');
%loglog(me*energy(1:Np)+m,Fp1(1:Np),'red', me*energy(1:Np)+m,Fp2(1:Np),'green',me*energy(1:Np)+m,Fp3(1:Np),'blue');
%loglog(me*energy(1:Np)/m+1,Fp1(1:Np),'red', me*energy(1:Np)/m+1,Fp2(1:Np),'green',me*energy(1:Np)/m+1,Fp3(1:Np),'blue');
%loglog(energy(1:Np)*m/me,Fp1(1:Np),'red', energy(1:Np)*m/me,Fp2(1:Np),'green',energy(1:Np)*m/me,Fp3(1:Np),'blue');
plot(energy1(1:Np),Fe1(1:Np),'--','Color','red');
plot(energy1(1:Np),Fe2(1:Np),'--','Color','green');
plot(energy1(1:Np),Fe3(1:Np),'--','Color','blue');
plot(energy2(1:Np),Fp1(1:Np),'red');
plot(energy2(1:Np),Fp2(1:Np),'green');
plot(energy2(1:Np),Fp3(1:Np),'blue');
title('F(E)');
%xlabel('Ekin/me c^2');
xlabel('E/me c^2');
ylabel('F(E)E^2');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
%legend('Fe', name,'Location','southeast');
legend('electrons t=0','electrons t=T/2','electrons t=T','protons t=0','protons t=T/2','protons t=T');
%legend('t=0','t=500','t=1000')
grid;
