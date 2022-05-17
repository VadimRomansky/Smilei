clear;
directory_name = './output/';
file_name = 'ParticleBinning7';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 6;
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
minE = minEp;
maxE = maxEp;
factor = (maxE/minE)^(1.0/(Np-1));

mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;

m = mp;

startPowerP = 125;
endPowerP = 145;

startPowerE = 149;
endPowerE = 159;

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

shockx = 38000;

startx = fix((shockx - 2560)/samplingFactor)+1;
endx = fix((shockx - 320)/samplingFactor);

startx = fix(1000/samplingFactor)+1;
endx = fix(20000/samplingFactor);

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

index1 = 15;
index2 = 30;

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
%T = 3*10^10;
theta = kB*T/(m*c*c);
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

x0(1:2)=[0.1,0.1];
error = evaluate_error(Fp2, 0.127, sqrt(1.0 - 1.0/2.25), gam, Np, index1, index2);
fun = @(x)evaluate_error(Fp2,x(1),x(2), gam, Np, index1, index2);
[x,fval] = fminunc(fun,x0);
for i = 1:Np,   
    Fshifted(i) = juttner_shifted_integrated(gam(i), x(1), x(2));
end;
T3=x(1)*m*c*c/kB;

Fpa(1:Np) = 0;

Fpa(startPower) = Fp2(startPower);
Fpa(endPower) = Fp2(endPower);

%gammap = log(Fpa(startPower)/Fpa(endPower))/log((me*energy(startPower)+m)/(me*energy(endPower)+m));

polyfitx(1:endPower-startPower + 1) = 0;
polyfity(1:endPower-startPower + 1) = 0;

for i = 1:endPower-startPower + 1,
    polyfitx(i) = log((me*energy(i+startPower - 1)+m));
    %polyfitx(i) = log((me*energy(i+startPower - 1)));
    polyfity(i) = log(Fp2(i+startPower - 1));
end;
p = polyfit(polyfitx, polyfity, 1);

%ap = exp(log(Fpa(startPower)) - gammap*log((me*energy(startPower)+m)));

for i = startPower-5:endPower+5,
    Fpa(i) = exp(polyval(p, log(me*energy(i)+m)));
    %Fpa(i) = exp(polyval(p, log(me*energy(i))));
end;



figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlim([1.0 10000]);
ylim([10^-10 1]);
plot(energy(1:Np)+m/me,Fp2(1:Np),'red','LineWidth',2);
plot(energy(1:Np)+m/me, Fjuttner(1:Np),'blue','LineWidth',2);
plot(energy(1:Np)+m/me,Fpa(1:Np),'green','LineWidth',2);
plot(energy(1:Np)+m/me, Fshifted(1:Np), 'black', 'LineWidth',2);
plot(energy(startPower) + m/me,Fp2(startPower),'o','Color','red');
plot(energy(endPower) + m/me,Fp2(endPower),'o','Color','red');

%plot(energy(1:Np),Fp2(1:Np),'red','LineWidth',2);
%plot(energy(1:Np), Fjuttner(1:Np),'blue','LineWidth',2);
%plot(energy(1:Np),Fpa(1:Np),'green','LineWidth',2);
%plot(energy(startPower),Fp2(startPower),'o','Color','red');
%plot(energy(endPower),Fp2(endPower),'o','Color','red');
title('F(E)');
xlabel('E/me c^2');
ylabel('F(E)');
name = strcat('powerlaw \gamma = ',num2str(p(1)));
%legend('Fe', 'maxwell-juttner',name,'Location','southeast');
name2 = strcat('juttner \theta = ', num2str(x(1)),' \beta = ', num2str(x(2)));
legend('Fe', 'maxwell-juttner',name, name2,'Location','southeast');
grid;

output(1:167,1:4) = 0;
for i = 1:167,
    output(i,1) = log10(energy(i) + m/me);
    output(i,2) = Fp2(i);
    output(i,3) = Fjuttner(i);
    output(i,4) = Fpa(i);
end;

dlmwrite('electrons.dat',output,'delimiter',' ');