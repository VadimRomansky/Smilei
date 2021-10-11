clear;
directory_name = './output/';
file_name = 'ParticleBinning6';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name1 = info.Datasets(1).Name;
name2 = info.Datasets(Ndata).Name;
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

mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;

m = me;

gam = 1.048;
beta = sqrt(1 - 1/(gam*gam));
c = 2.99792458*10^10;
Te = 2.6*10^9;
Temin = 10^9;
Temax = 2*10^11;
Tp = 2*10^11;
Tpmin = 10^9;
Tpmax = 10^13;

T = Te;
Tmax = Temax;
Tmin = Temin;
kB = 1.3806488*10^-16;
theta = kB*T/(m*c*c);

energy(1:Np) = 0;
de(1:Np) = 0;
energy(1) = minE;
Fjuttner(1:Np) = 0;
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
        Fp1(i)=Fp1(i)+fp1(i,j)/de(i);
        Fp2(i)=Fp2(i)+fp2(i,j)/de(i);
    end;
end;

norm = 1.0;
normp = 0.0;
for i = 1:Np,
    normp = normp + Fp2(i)*de(i)*me/m;
end;

for i = 1:Np,
    Fp2(i) = Fp2(i)*norm/normp;
end;

index1 = 50;
index2 = 100;

Tleft = Tmin;
Tright = Tmax;

for j = 1:20,
    T1 = Tleft + (Tright - Tleft)/3;
    T2 = Tleft + (Tright - Tleft)*2/3;
    s1 = 0;
    s2 = 0;
    theta = kB*T1/(m*c*c);
    bes = besselk(2, 1/theta);
    for i = index1:index2,
        gam = energy(i)*me/m + 1;
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner(i) = (1.0/(theta*bes))*exp1*gam*gam*beta;
        s1 = s1 + ((Fjuttner(i) - Fp2(i))^2)*de(i);
    end;
    theta = kB*T2/(m*c*c);
    bes = besselk(2, 1/theta);
    for i = index1:index2,
        gam = energy(i)*me/m + 1;
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner(i) = (1.0/(theta*bes))*exp1*gam*gam*beta;
        s2 = s2 + ((Fjuttner(i) - Fp2(i))^2)*de(i);
    end;
    if(s1 < s2)
        Tright = T2;
    else 
        Tleft = T1;
    end;
end;
T = (Tleft + Tright)/2;
theta = kB*T/(m*c*c);
bes = besselk(2, 1/theta);
for i = 1:Np,   
    gam = energy(i)*me/m + 1;
    beta = sqrt(1.0 - 1.0/(gam*gam));
    exp1 = exp(-gam/theta);       
    Fjuttner(i) = (1.0/(theta*bes))*exp1*gam*gam*beta;
end;

startPowerP = 95;
endPowerP = 115;

startPowerE = 140;
endPowerE = 150;

startPower = startPowerE;
endPower = endPowerE;

Fpa(1:Np) = 0;

Fpa(startPower) = Fp2(startPower);
Fpa(endPower) = Fp2(endPower);

gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));


ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

for i = startPower-20:endPower+20,
    Fpa(i) = ap*((me*energy(i)+m)^gammap);
end;



figure(1);
hold on;
plot(energy(1:Np)+m/me,Fp2(1:Np),'red','LineWidth',2);
plot(energy(1:Np)+m/me, Fjuttner(1:Np),'blue','LineWidth',2);
plot(energy(1:Np)+m/me,Fpa(1:Np),'green','LineWidth',2);
title('F(E)');
xlabel('E/me c^2');
ylabel('F(E)');
name = strcat('powerlaw \gamma = ',num2str(gammap));
legend('Fe', 'maxwell-juttner',name,'Location','southeast');
grid;

dlmwrite('Ee9.dat',energy,'delimiter','\n');
dlmwrite('Fe9.dat',Fp2,'delimiter','\n');