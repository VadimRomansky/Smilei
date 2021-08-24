clear;
directory_name = './output_theta0-90_gamma1.5sigma0.004/';
file_name = 'ParticleBinning7';
file_ending = '.h5';

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Nd = 10;
start = 0;
fileNumber = 6;
full_name = strcat(directory_name, file_name, num2str(fileNumber), file_ending);

info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 130;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

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

startx(1:Nd) = 0;
endx(1:Nd) = 0;
Fp(1:Nd,1:Np)=0;

samplingFactor = 20;

for i = 1:Nd,
    startx(i) = fix(20000/samplingFactor)+1;
    endx(i) = fix(30000/samplingFactor);
end;

for k = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(k+start-1), file_ending);
    info = h5info(full_name);
    name = info.Datasets(Ndata).Name;
    fp = hdf5read(full_name, name);
    for i=1:Np,
        for j=startx:endx,
            Fp(k,i)=Fp(k, i)+fp(i,j)/de(i);
        end;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1.5);

figure(1);
hold on;
title ('F_p(E)');
xlabel ('E_{kin}/{m_p c^2}');
ylabel ('F_p(E)');
for j=1:Nd,
    plot (energy(1:Np)/100,100*Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');
grid ;

tempOutput(1:Np, 1:2) = 0;
for i = 1:Np,
    %tempOutput(i) = Fp(1,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(1,i)*100;
end;

dlmwrite('Ap0.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep0.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp0.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(2,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(2,i)*100;
end;
dlmwrite('Ap1.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep1.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp1.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(3,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(3,i)*100;
end;
dlmwrite('Ap2.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep2.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp2.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(4,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(4,i)*100;
end;
dlmwrite('Ap3.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep3.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp3.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(5,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(5,i)*100;
end;
dlmwrite('Ap4.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep4.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp4.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(6,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(6,i)*100;
end;
dlmwrite('Ap5.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep5.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp5.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(7,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(7,i)*100;
end;
dlmwrite('Ap6.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep6.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp6.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(8,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(8,i)*100;
end;
dlmwrite('Ap7.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep7.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp7.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(9,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(9,i)*100;
end;
dlmwrite('Ap8.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep8.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp8.dat',tempOutput*100,'delimiter',' ');

for i = 1:Np,
    %tempOutput(i) = Fp(10,i);
    tempOutput(i,1) = energy(i)/100;
    tempOutput(i,2) = Fp(10,i)*100;
end;
dlmwrite('Ap9.dat',tempOutput,'delimiter',' ');
%dlmwrite('Ep9.dat',energy/100,'delimiter',' ');
%dlmwrite('Fp9.dat',tempOutput*100,'delimiter',' ');