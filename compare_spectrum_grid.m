clear;
directory_name = './output/';
file_name = 'ParticleBinning6';
file_ending = '.h5';

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Nd = 10;
start = 0;
fileNumber = start;
full_name = strcat(directory_name, file_name, num2str(fileNumber), file_ending);

info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 130;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minE = 0.1;
maxE = 1000;
factor = (maxE/minE)^(1.0/(Np-1));

energy(1:Np) = 0;
energy(1) = minE;
for i = 2:Np,
    energy(i) = energy(i-1)*factor;
end;

startx(1:Nd) = 0;
endx(1:Nd) = 0;
Fp(1:Nd,1:Np)=0;

samplingFactor = 20;

for i = 1:Nd,
    startx(i) = fix(10000/samplingFactor)+1;
    endx(i) = fix(15000/samplingFactor);
end;

for k = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(k+start-1), file_ending);
    name = info.Datasets(Ndata).Name;
    fp = hdf5read(full_name, name);
    for i=1:Np,
        for j=startx:endx,
            Fp(k,i)=Fp(k, i)+fp(i,j);
        end;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1.5);

figure(1);
hold on;
title ('F_e(E)');
xlabel ('E/{m_e c^2}');
ylabel ('F_e(E)*E^2');
for j=1:Nd,
    plot (energy(1:Np),Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');
grid ;
dlmwrite('Ee0.dat',energy,'delimiter','\n');
dlmwrite('Fs0.dat',Fp(1),'delimiter','\n');
dlmwrite('Ee1.dat',energy,'delimiter','\n');
dlmwrite('Fs1.dat',Fp(2),'delimiter','\n');
dlmwrite('Ee2.dat',energy,'delimiter','\n');
dlmwrite('Fs2.dat',Fp(3),'delimiter','\n');
dlmwrite('Ee3.dat',energy,'delimiter','\n');
dlmwrite('Fs3.dat',Fp(4),'delimiter','\n');
dlmwrite('Ee4.dat',energy,'delimiter','\n');
dlmwrite('Fs4.dat',Fp(5),'delimiter','\n');
dlmwrite('Ee5.dat',energy,'delimiter','\n');
dlmwrite('Fs5.dat',Fp(6),'delimiter','\n');
dlmwrite('Ee6.dat',energy,'delimiter','\n');
dlmwrite('Fs6.dat',Fp(7),'delimiter','\n');
dlmwrite('Ee7.dat',energy,'delimiter','\n');
dlmwrite('Fs7.dat',Fp(8),'delimiter','\n');
dlmwrite('Ee8.dat',energy,'delimiter','\n');
dlmwrite('Fs8.dat',Fp(9),'delimiter','\n');
dlmwrite('Ee9.dat',energy,'delimiter','\n');
dlmwrite('Fs9.dat',Fp(10),'delimiter','\n');