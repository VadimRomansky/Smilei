clear;
directory_name = './output/';
file_name = 'ParticleBinning7';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 11;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(Ndata).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);

Np=size(fp1,1);
Nx=size(fp1,2);
Ns = 5;

minEe = 0.001;
maxEe = 1000;
minEp = 0.1;
maxEp = 5000;
minE = minEp;
maxE = maxEp;
factor = (maxE/minE)^(1.0/(Np-1));

mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;

m = mp;

gam = 1.048;
beta = sqrt(1 - 1/(gam*gam));
c = 2.99792458*10^10;
Te = 2.6*10^9;
Temin = 10^8;
Temax = 2*10^12;
Tp = 2*10^11;
Tpmin = 10^9;
Tpmax = 10^13;

T = Tp;
Tmax = Tpmax;
Tmin = Tpmin;
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

Fp(1:Np,1:Ns)=0;

samplingFactor = 20;

startx(1:Ns) = 0; 
endx(1:Ns) = 0; 

startx(1) = 5000/samplingFactor+1;
endx(1) = 12000/samplingFactor;

startx(2) = 19000/samplingFactor+1;
endx(2) = 26000/samplingFactor;

startx(3) = 30000/samplingFactor+1;
endx(3) = 37000/samplingFactor;

startx(4) = 43000/samplingFactor+1;
endx(4) = 50000/samplingFactor;

startx(5) = 55000/samplingFactor+1;
endx(5) = 62000/samplingFactor;


startx(1)= 100000/samplingFactor;
endx(1) = fix(150000/samplingFactor);

for i =2:Ns,
   startx(i) = startx(1);
   endx(i) = endx(1);
end;

Nt = [1,2,5,11,24];
for k = 1:Ns,
    name = info.Datasets(Nt(k)).Name;
    fp= hdf5read(full_name, name);
    norm = 0;
    for i=1:Np,
        for j=startx:endx,
            Fp(i,k)=Fp(i,k)+fp(i,j)/de(i);
            norm = norm + fp(i,j);
        end;
    end;
    for i=1:Np,
        Fp(i,k) = Fp(i,k)/norm;
    end;
   

    index1 = 10;
    index2 = 80;

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
            s1 = s1 + ((Fjuttner(i) - Fp(i,k))^2)*de(i);
        end;
        theta = kB*T2/(m*c*c);
        bes = besselk(2, 1/theta);
        for i = index1:index2,
            gam = energy(i)*me/m + 1;
            beta = sqrt(1.0 - 1.0/(gam*gam));
            exp1 = exp(-gam/theta);       
            Fjuttner(i) = (1.0/(theta*bes))*exp1*gam*gam*beta;
            s2 = s2 + ((Fjuttner(i) - Fp(i,k))^2)*de(i);
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

    startPowerP = 115;
    endPowerP = 135;

    startPowerE = 138;
    endPowerE = 147;

    startPower = startPowerE;
    endPower = endPowerE;

    Fpa(1:Np) = 0;

    Fpa(startPower) = Fp(startPower,k);
    Fpa(endPower) = Fp(endPower,k);

    %gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));

    polyfitx(1:endPower-startPower + 1) = 0;
    polyfity(1:endPower-startPower + 1) = 0;

    for i = 1:endPower-startPower + 1,
        polyfitx(i) = log((me*energy(i+startPower - 1)+m));
        polyfity(i) = log(Fp(i+startPower - 1,k));
    end;
    p = polyfit(polyfitx, polyfity, 1);

    %ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

    for i = startPower-5:endPower+5,
        %Fpa(i) = ap*((me*energy(i)+m)^gammap);
        Fpa(i) = exp(polyval(p, log(me*energy(i)+m)));
    end;


    figure(k);
    hold on;
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlim([1.0 1000]);
    ylim([10^-10 10]);
    plot(energy(1:Np)+m/me, Fp(1:Np,k),'red','LineWidth',2);
    plot(energy(1:Np)+m/me, Fjuttner(1:Np),'blue','LineWidth',2);
    plot(energy(1:Np)+m/me, Fpa(1:Np),'green','LineWidth',2);
    plot(energy(startPower) + m/me, Fp(startPower,k),'o','Color','red');
    plot(energy(endPower) + m/me, Fp(endPower,k),'o','Color','red');
    title('F(E)');
    xlabel('E/me c^2');
    ylabel('F(E)');
    name = strcat('powerlaw \gamma = ',num2str(p(1)));
    legend('Fe', 'maxwell-juttner',name,'Location','southeast');
    grid;
end;