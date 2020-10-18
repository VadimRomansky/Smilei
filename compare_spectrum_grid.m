clear;
directory_name = './output1/';
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
Ndata = 90;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

startx(1:Nd) = 0;
endx(1:Nd) = 0;
Fp(1:Nd,1:Np)=0;

samplingFactor = 20;

for i = 1:Nd,
    startx(i) = fix(30000/samplingFactor)+1;
    endx(i) = fix(40000/samplingFactor);
end;

for k = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(k+start-1), file_ending);
    name = info.Datasets(Ndata).Name;
    fp = hdf5read(full_name, name);
    for i=1:Np,
        for j=startx:endx,
            Fp(k,i)=Fp(k, i)+fp(i,j)*i*i;
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
    plot (1:Np,Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');
grid ;