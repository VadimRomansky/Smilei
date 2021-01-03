clear;
directory_name = './output_theta0-90/';
file_name = 'ParticleBinning6';
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
    startx(i) = fix(7000/samplingFactor)+1;
    endx(i) = fix(15000/samplingFactor);
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
title ('F_e(E)');
xlabel ('E/{m_e c^2}');
ylabel ('F_e(E)');
for j=1:Nd,
    plot (energy(1:Np),Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');
grid ;

tempOutput(1:Np) = 0;
for i = 1:Np,
    tempOutput(i) = Fp(1,i);
end;

dlmwrite('Ee0.dat',energy,'delimiter',' ');
dlmwrite('Fs0.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(2,i);
end;
dlmwrite('Ee1.dat',energy,'delimiter',' ');
dlmwrite('Fs1.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(3,i);
end;
dlmwrite('Ee2.dat',energy,'delimiter',' ');
dlmwrite('Fs2.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(4,i);
end;
dlmwrite('Ee3.dat',energy,'delimiter',' ');
dlmwrite('Fs3.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(5,i);
end;
dlmwrite('Ee4.dat',energy,'delimiter',' ');
dlmwrite('Fs4.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(6,i);
end;
dlmwrite('Ee5.dat',energy,'delimiter',' ');
dlmwrite('Fs5.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(7,i);
end;
dlmwrite('Ee6.dat',energy,'delimiter',' ');
dlmwrite('Fs6.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(8,i);
end;
dlmwrite('Ee7.dat',energy,'delimiter',' ');
dlmwrite('Fs7.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(9,i);
end;
dlmwrite('Ee8.dat',energy,'delimiter',' ');
dlmwrite('Fs8.dat',tempOutput,'delimiter',' ');

for i = 1:Np,
    tempOutput(i) = Fp(10,i);
end;
dlmwrite('Ee9.dat',energy,'delimiter',' ');
dlmwrite('Fs9.dat',tempOutput,'delimiter',' ');