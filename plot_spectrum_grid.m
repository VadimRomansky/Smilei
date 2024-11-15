clear;
%directory_name = './output_gamma0.3_sigma0.0002_dx0.2_theta0-90/';
directory_name = './output/';
file_name = 'ParticleBinning6';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
<<<<<<< HEAD
%Ndata = 10;
=======
Ndata = 10;
>>>>>>> a20194666d8315ed0adc456b2063cc0475b72fc6
name1 = info.Datasets(1).Name;
name2 = info.Datasets(fix(Ndata/2)+1).Name;
name3 = info.Datasets(Ndata).Name;

fp1= hdf5read(full_name, name1);
fp2 = hdf5read(full_name, name2);
fp3 = hdf5read(full_name, name3);

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
Fp3(1:Np)=0;

samplingFactor = 20;

<<<<<<< HEAD
startx = fix(100/samplingFactor)+1;
endx = fix(Nx);
=======
startx = fix(70000/samplingFactor)+1;
endx = fix(80000/samplingFactor);
>>>>>>> a20194666d8315ed0adc456b2063cc0475b72fc6

for i=1:Np,
    for j=startx:endx,
        Fp1(i)=Fp1(i)+fp1(i,j)/de(i);
        Fp2(i)=Fp2(i)+fp2(i,j)/de(i);
        Fp3(i)=Fp3(i)+fp3(i,j)/de(i);
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

for i=1:Np,
        Fp1(i)=Fp1(i);
        Fp2(i)=Fp2(i);
        Fp3(i)=Fp3(i);
end;

set(0, 'DefaultLineLineWidth', 2);
figure(2);
%plot(energy(1:Np),Fp2(1:Np),'red', energy(startPowerP:endPowerP), Fpa(startPowerP:endPowerP),'blue');
%loglog(energy(1:Np),Fp2(1:Np),'red');
%loglog(me*energy(1:Np)+m,Fp1(1:Np),'red', me*energy(1:Np)+m,Fp2(1:Np),'green',me*energy(1:Np)+m,Fp3(1:Np),'blue');
%loglog(me*energy(1:Np)/m+1,Fp1(1:Np),'red', me*energy(1:Np)/m+1,Fp2(1:Np),'green',me*energy(1:Np)/m+1,Fp3(1:Np),'blue');
%loglog(energy(1:Np)*m/me,Fp1(1:Np),'red', energy(1:Np)*m/me,Fp2(1:Np),'green',energy(1:Np)*m/me,Fp3(1:Np),'blue');
loglog(energy(1:Np),Fp1(1:Np),'red', energy(1:Np),Fp2(1:Np),'green',energy(1:Np),Fp3(1:Np),'blue');
title('F(E)');
%xlabel('Ekin/me c^2');
xlabel('E/me c^2');
ylabel('F(E)');
name = strcat('approximation gamma = ',num2str(gammap-2));
%legend('Fe', name,'Location','southeast');
legend('t=0','t=T/2','t=T')
%legend('t=0','t=500','t=1000')
grid;

gammaKin(1:Np)=0;
for i=1:Np,
    gammaKin(i) = energy(i)*me/m;
end;

dlmwrite('Ee0.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs0.dat',Fp3,'delimiter',' ');
dlmwrite('Ee1.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs1.dat',Fp3,'delimiter',' ');
dlmwrite('Ee2.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs2.dat',Fp3,'delimiter',' ');
dlmwrite('Ee3.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs3.dat',Fp3,'delimiter',' ');
dlmwrite('Ee4.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs4.dat',Fp3,'delimiter',' ');
dlmwrite('Ee5.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs5.dat',Fp3,'delimiter',' ');
dlmwrite('Ee6.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs6.dat',Fp3,'delimiter',' ');
dlmwrite('Ee7.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs7.dat',Fp3,'delimiter',' ');
dlmwrite('Ee8.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs8.dat',Fp3,'delimiter',' ');
dlmwrite('Ee9.dat',gammaKin,'delimiter',' ');
dlmwrite('Fs9.dat',Fp3,'delimiter',' ');