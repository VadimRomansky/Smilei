clear;
directory_name = './output/';
file_name = 'ParticleBinning6';
field_file_name = 'Fields0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
field_full_name = strcat(directory_name, field_file_name, file_number);
info = h5info(full_name);
field_info = h5info(field_full_name);
Ndata = size(info.Datasets,1);
%Ndata = 7;
name1 = info.Datasets(1).Name;
name2 = info.Datasets(Ndata).Name;
fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);

name2x = strcat(field_info.Groups.Groups(Ndata).Name, '/Bx');
name2y = strcat(field_info.Groups.Groups(Ndata).Name, '/By');
name2z = strcat(field_info.Groups.Groups(Ndata).Name, '/Bz');

Bx= hdf5read(field_full_name, name2x);
By= hdf5read(field_full_name, name2y);
Bz= hdf5read(field_full_name, name2z);

Np=size(fp1,1);
Nx=size(fp1,2);

Ns =10;
shockX = 60000;
boxSize = shockX/Ns;

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

startPowerP = 95;
endPowerP = 115;

startPowerE = 138;
endPowerE = 148;

startPower = startPowerE;
endPower = endPowerE;
Fpa(1:Np,1:Ns) = 0;

energy(1:Np) = 0;
de(1:Np) = 0;
energy(1) = minE;
Fjuttner(1:Np,1:Ns) = 0;

for i = 2:Np,
    energy(i) = energy(i-1)*factor;
end;
de(1) = energy(2) - energy(1);
for i = 2:Np,
    de(i) = energy(i) - energy(i-1);
end;

Fp1(1:Np,1:Ns)=0;
Fp2(1:Np,1:Ns)=0;

samplingFactor = 20;
fieldSamplingFactor = 4;
fieldNy = size(Bx,1);

Ts(1:Ns) = 0;
gammap(1:Ns) = 0;
Bxa(1:Ns) = 0;
Bya(1:Ns) = 0;
Bza(1:Ns) = 0;
Ba(1:Ns) = 0;


for k = 1:Ns,
    startx = fix((k-1)*boxSize/samplingFactor)+1;
    endx = fix(k*boxSize/samplingFactor);
    
    fieldstartx = fix((k-1)*boxSize/fieldSamplingFactor)+1;
    fieldendx = fix(k*boxSize/fieldSamplingFactor);
    
    B2 = 0;
    
    for j=fieldstartx:fieldendx,
        for l = 1:fieldNy,
            Bxa(k) = Bxa(k) + Bx(l,j)/(fieldNy*(fieldendx - fieldstartx + 1));
            Bya(k) = Bya(k) + By(l,j)/(fieldNy*(fieldendx - fieldstartx + 1));
            Bza(k) = Bza(k) + Bz(l,j)/(fieldNy*(fieldendx - fieldstartx + 1));
            B2 = B2 + (Bx(l,j)*Bx(l,j) + By(l,j)*By(l,j) + Bz(l,j)*Bz(l,j))/(fieldNy*(fieldendx - fieldstartx + 1));
        end;
    end;
    Ba(k) = sqrt(B2);
    
    for i=1:Np,
        for j=startx:endx,
            Fp1(i,k)=Fp1(i,k)+fp1(i,j)/de(i);
            Fp2(i,k)=Fp2(i,k)+fp2(i,j)/de(i);
        end;
    end;

    norm = 1.0;
    normp = 0.0;
    for i = 1:Np,
        normp = normp + Fp2(i,k)*de(i)*me/m;
    end;

    for i = 1:Np,
        Fp2(i,k) = Fp2(i,k)*norm/normp;
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
            Fjuttner(i,k) = (1.0/(theta*bes))*exp1*gam*gam*beta;
            s1 = s1 + ((Fjuttner(i,k) - Fp2(i,k))^2)*de(i);
        end;
        theta = kB*T2/(m*c*c);
        bes = besselk(2, 1/theta);
        for i = index1:index2,
            gam = energy(i)*me/m + 1;
            beta = sqrt(1.0 - 1.0/(gam*gam));
            exp1 = exp(-gam/theta);       
            Fjuttner(i,k) = (1.0/(theta*bes))*exp1*gam*gam*beta;
            s2 = s2 + ((Fjuttner(i,k) - Fp2(i,k))^2)*de(i);
        end;
        if(s1 < s2)
            Tright = T2;
        else 
            Tleft = T1;
        end;
    end;
    T = (Tleft + Tright)/2;
    Ts(k) = T;
    theta = kB*T/(m*c*c);
    bes = besselk(2, 1/theta);
    for i = 1:Np,   
        gam = energy(i)*me/m + 1;
        beta = sqrt(1.0 - 1.0/(gam*gam));
        exp1 = exp(-gam/theta);       
        Fjuttner(i,k) = (1.0/(theta*bes))*exp1*gam*gam*beta;
    end;

    startPower = startPowerE;
    endPower = endPowerE;

    %gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));

    polyfitx(1:endPower-startPower + 1) = 0;
    polyfity(1:endPower-startPower + 1) = 0;

    for i = 1:endPower-startPower + 1,
        polyfitx(i) = log((me*energy(i+startPower - 1)+m));
        polyfity(i) = log(Fp2(i+startPower - 1,k));
    end;
    p = polyfit(polyfitx, polyfity, 1);
    gammap(k) = p(1);

    %ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

    for i = startPower-20:endPower+20,
        %Fpa(i) = ap*((me*energy(i)+m)^gammap);
        Fpa(i,k) = exp(polyval(p, log(me*energy(i)+m)));
    end;



    figure(k);
    hold on;
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    plot(energy(1:Np)+m/me,Fp2(1:Np,k),'red','LineWidth',2);
    plot(energy(1:Np)+m/me, Fjuttner(1:Np,k),'blue','LineWidth',2);
    plot(energy(1:Np)+m/me,Fpa(1:Np,k),'green','LineWidth',2);
    plot(energy(startPower) + m/me,Fp2(startPower,k),'o','Color','red');
    plot(energy(endPower) + m/me,Fp2(endPower,k),'o','Color','red');
    title('F(E)');
    xlabel('E/me c^2');
    ylabel('F(E)');
    name = strcat('powerlaw \gamma = ',num2str(gammap(k)));
    legend('Fe', 'maxwell-juttner',name,'Location','southeast');
end;

figure(Ns+1);
hold on;
plot(1:Ns, Ts(1:Ns),'red','LineWidth',2);
title('T(x)');
xlabel('x');
ylabel('T');
grid;

figure(Ns+2);
hold on;
plot(1:Ns, gammap(1:Ns),'red','LineWidth',2);
title('\gamma');
xlabel('x');
ylabel('\gamma');
grid;

figure(Ns+3);
hold on;
plot(1:Ns, Bxa(1:Ns),'red','LineWidth',2);
plot(1:Ns, Bya(1:Ns),'green','LineWidth',2);
plot(1:Ns, Bza(1:Ns),'blue','LineWidth',2);
plot(1:Ns, Ba(1:Ns),'black','LineWidth',2);
legend('Bx','By','Bz','B_{ms}');
title('B');
xlabel('x');
ylabel('B');
grid;

