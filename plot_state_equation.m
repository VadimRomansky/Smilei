clear;
%directory_name = './output_gamma0.3_sigma0.0002_theta30/';
directory_name = './output/';
file_name = 'ParticleBinning6';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 5;
name = info.Datasets(Ndata).Name;
fp2 = hdf5read(full_name, name);

Np=size(fp2,1);
Nx=size(fp2,2);

file_name = 'ParticleBinning0';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
name3 = info.Datasets(Ndata).Name;
fp3 = hdf5read(full_name, name3);

dx = 0.2;
factor = 1.0/(dx*dx);
Ny = 200;
N=size(fp3,1);


samplingFactor = fix(N/Nx);

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

startPowerP = 125;
endPowerP = 145;

startPowerE = 137;
endPowerE = 148;

startPower = startPowerE;
endPower = endPowerE;

gam = 1.048;
beta = sqrt(1 - 1/(gam*gam));
c = 2.99792458*10^10;
Te = 2.6*10^9;
Temin = 10^7;
Temax = 2*10^12;
Tp = 2*10^11;
Tpmin = 10^8;
Tpmax = 10^12;

Tmax = Temax;
Tmin = Temin;
kB = 1.3806488*10^-16;

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

Nd = 10;
T(1:Nd) = Te;
n(1:Nd) = 0;

startx(1:Nd) = 0;
endx(1:Nd) = 0;
length = 1500;

startx(1) = 1000;
endx(1) = startx(1) + length;
for i=2:Nd,
    startx(i) = endx(i-1);
    endx(i) = startx(i) + length;
end;

for k=1:Nd,
    for j = startx(k):endx(k),
        n(k) = n(k) + fp3(j)*factor*samplingFactor/(Ny*length);
    end;

    for i=1:Np,
        for j=startx(k)/samplingFactor:endx(k)/samplingFactor,
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

    %for debug
    %Fp2(Np)=10^-10;

    %remove zeros
    for i = 2:Np-1,
        if(Fp2(i) <= 0) & (Fp2(i-1) > 0) & (Fp2(i+1) > 0)
            Fp2(i) = 0.5*(Fp2(i-1)+Fp2(i+1));
        end;
    end;

    maxNonZero = 0;
    for j = Np:-1:1,
        if(Fp2(j) > 0)
            maxNonZero = j;
            break;
        end;
    end;
    minNonZero = Np;
    for i = 1:Np
        if(Fp2(i) > 0)
            minNonZero = i;
            break;
        end;
    end;

%for i = minNonZero:maxNonZero-1,
%    if(Fp2(i) <= 0)
%        nextNonZero = i+1;
%        while(Fp2(nextNonZero)<=0)
%            nextNonZero = nextNonZero + 1;
%        end;
%        s = (log(Fp2(i-1))- log(Fp2(nextNonZero)))/(log(me*energy(i-1)+m) - log(me*energy(nextNonZero)+m));
%        for j = i:nextNonZero-1,
%            Fp2(j) = Fp2(i-1)*power((me*energy(j)+m)/(me*energy(i-1)+m),s);
%        end;
%    end;
%end;

    index1 = 40;
    index2 = 70;

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
    T(k) = (Tleft + Tright)/2;
    theta = kB*T(k)/(m*c*c);
    bes = besselk(2, 1/theta);
    Fshifted(1:Np) = 0;
    gam(1:Np) = 1.0;
    for i = 1:Np,   
        gam(i) = energy(i)*me/m + 1;
        beta = sqrt(1.0 - 1.0/(gam(i)*gam(i)));
        exp1 = exp(-gam(i)/theta);       
        Fjuttner(i) = (1.0/(theta*bes))*exp1*gam(i)*gam(i)*beta;
        Fshifted(i) = juttner_shifted_integrated(gam(i), 0.127, sqrt(1.0 - 1.0/2.25));
    end;
 
    normShifted = 0;
    for i = 1:Np,
        normShifted = normShifted + Fshifted(i)*de(i)*me/m;
    end;
 
    x0(1:2)=[theta,0.3];
    %error = evaluate_error(Fp2, 0.127, sqrt(1.0 - 1.0/2.25), gam, Np, index1, index2);
    fun = @(x)evaluate_error(Fp2,x(1),x(2), gam, Np, index1, index2);
    [x,fval] = fminunc(fun,x0);
    for i = 1:Np,   
        Fshifted(i) = juttner_shifted_integrated(gam(i), x(1), x(2));
        %Fshifted(i) = juttner_shifted_integrated(gam(i), 0.02, 0.3);
    end;
    T(k)=x(1)*m*c*c/kB;

    figure(k);
    hold on;
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlim([1.0 10000]);
    ylim([10^-10 1]);
    plot(energy(1:Np)+m/me,Fp2(1:Np),'red','LineWidth',2);
    plot(energy(1:Np)+m/me, Fshifted(1:Np), 'black', 'LineWidth',2);
end;

figure(Nd+1);
colormap Jet;
c = linspace(1,10,Nd);
scatter(n(1:Nd),T(1:Nd),[],c,'filled');
title('T(n)');
xlabel('n/n_0');
ylabel('T');