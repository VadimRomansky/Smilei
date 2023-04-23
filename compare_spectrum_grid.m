clear;
directory_name = './output_gamma0.3_sigma0.0002_dx0.2_theta0-90/';
file_name = 'ParticleBinning6';
file_ending = '.h5';

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Nd = 10;
start = 0;
fileNumber = 1;
full_name = strcat(directory_name, file_name, num2str(fileNumber), file_ending);

info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 15;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

mass_ratio = 100;
me=1;
mp=me*mass_ratio;

minEe = 0.001;
maxEe = 5000;
minEp = 0.1;
maxEp = 5000;
minE = minEe;
maxE = maxEe;
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
    startx(i) = fix(2000/samplingFactor)+1;
    endx(i) = fix(50000/samplingFactor);
end;
for i = 5:Nd,
    startx(i) = fix(35000/samplingFactor)+1;
    endx(i) = fix(45000/samplingFactor);
end;
% endx(1) = fix(51000/samplingFactor);
% endx(2) = fix(50000/samplingFactor);
% endx(3) = fix(49000/samplingFactor);
% endx(4) = fix(48000/samplingFactor);
% endx(5) = fix(49000/samplingFactor);
% endx(6) = fix(50000/samplingFactor);
% endx(7) = fix(52000/samplingFactor);
% endx(8) = fix(53000/samplingFactor);
% endx(9) = fix(54000/samplingFactor);
% endx(10) = fix(55000/samplingFactor);
% for i = 1:Nd,
%     startx(i) = endx(i) - fix(10000/samplingFactor);
% end;

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
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
for j=1:Nd,
    plot (energy(1:Np)/100,100*Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');
grid ;

for k = 1:Nd,
    %tempOutput(1:Np, 1:2) = 0;
    tempOutput(1:Np)=0;
    for i = 1:Np,
         tempOutput(i) = Fp(k,i)*m/me;
         %tempOutput(i,1) = me*energy(i)/m;
         %tempOutput(i,2) = Fp(k,i)*m/me;
     end;
 
     %dlmwrite('Ap' + num2str(k) + 'dat',tempOutput,'delimiter',' ');
     dlmwrite('Ee' + num2str(k) + '.dat',energy/100,'delimiter',' ');
     dlmwrite('Fs' + num2str(k) + '.dat',tempOutput,'delimiter',' ');
end;