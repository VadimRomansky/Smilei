clear;
directory_name = './outputs/';
file_name = 'ParticleBinning0';
file_ending = '.h5';

Number = {0, 2, 3, 4, 8};
Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Nd = 5;
start = 0;
fileNumber = start;
full_name = strcat(directory_name, file_name, num2str(Number{1}), file_ending);

info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 100;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Nx=size(fp,1);


Fp(1:Nd,1:Nx)=0;


for k = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(Number{k}), file_ending);
    name = info.Datasets(Ndata).Name;
    fp = hdf5read(full_name, name);
    for i=1:Nx,
        Fp(k,i)= fp(i);
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1.5);

figure(1);
hold on;
title ('n(x)');
xlabel ('x');
ylabel ('n');
for j=1:Nd,
    plot (1:Nx,Fp(j, 1:Nx),'color',Color{j});
end;

legend(LegendTitle{Number{1}+1}, LegendTitle{Number{2}+1}, LegendTitle{Number{3}+1}, LegendTitle{Number{4}+1}, LegendTitle{Number{5}+1},'Location','northwest');
grid ;