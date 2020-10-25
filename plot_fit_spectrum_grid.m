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


mp = 1.67*10^-24;
mass_ratio = 64;
me = mp/mass_ratio;
gam = 1.5;
beta = sqrt(1 - 1/(gam*gam));
c = 2.99792458*10^10;
Te = 2.6*10^9;
Temin = 10^9;
Temax = 10^12;
Tp = 2*10^9;
Tpmin = 10^8;
Tpmax = 10^11;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);

minE = 0.1;
maxE = 1000;
factor = (maxE/minE)^(1.0/(Np-1));

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

startx = fix(10000/samplingFactor)+1;
endx = fix(17000/samplingFactor);

for i=1:Np,
    for j=startx:endx,
        Fp1(i)=Fp1(i)+fp1(i,j)/de(i);
        Fp2(i)=Fp2(i)+fp2(i,j)/de(i);
    end;
end;

norm = 1.0;
normp = 0;
for i = 1:Np,
    normp = normp + Fp2(i)*de(i);
end;

for i = 1:Np,
    Fp2(i) = Fp2(i)*norm/normp;
end;

index1 = 50;
index2 = 110;

Teleft = Temin;
Teright = Temax;

for j = 1:20,
    Te1 = Teleft + (Teright - Teleft)/3;
    Te2 = Teleft + (Teright - Teleft)*2/3;
    s1 = 0;
    s2 = 0;
    thetae = kB*Te1/(me*c*c);
    bes = besselk(2, 1/thetae);
    for i = index1:index2,
        gam = energy(i) + 1;
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-(energy(i) + 1.0)/thetae);       
        Fjuttner(i) = (1.0/(thetae*bes))*exp1*gam*gam*beta;
        s1 = s1 + ((Fjuttner(i) - Fp2(i))^2)*de(i);
    end;
    thetae = kB*Te2/(me*c*c);
    bes = besselk(2, 1/thetae);
    for i = index1:index2,
        gam = energy(i) + 1;
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-(energy(i) + 1.0)/thetae);       
        Fjuttner(i) = (1.0/(thetae*bes))*exp1*gam*gam*beta;
        s2 = s2 + ((Fjuttner(i) - Fp2(i))^2)*de(i);
    end;
    if(s1 < s2)
        Teright = Te2;
    else 
        Teleft = Te1;
    end;
end;
Te = (Teleft + Teright)/2;
thetae = kB*Te/(me*c*c);
bes = besselk(2, 1/thetae);
for i = 1:Np,   
    gam = energy(i) + 1;
    beta = sqrt(1.0 - 1.0/(gam*gam));
    exp1 = exp(-(energy(i) + 1.0)/thetae);       
    Fjuttner(i) = (1.0/(thetae*bes))*exp1*gam*gam*beta;
end;


figure(1);
plot(energy(1:Np),Fp2(1:Np),'red', energy(1:Np), Fjuttner(1:Np),'blue');
title('F(E)');
xlabel('Ekin/me c^2');
ylabel('F(E)');
legend('Fe', 'maxwell-juttner','Location','southeast');
grid;