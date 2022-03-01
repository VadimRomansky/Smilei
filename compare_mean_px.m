clear;
directory_name = './output4/';
file_name = 'ParticleBinning10';
file_ending = '.h5';

Number = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
%LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};
%LegendTitle = {'{\gamma} = 1.05','{\gamma} = 1.5','{\gamma} = 3', '{\gamma} = 7','{\gamma} = 10','{\gamma} = 20'};
LegendTitle = {'T = 0.00001','T = 0.0001','T = 0.001', 'T = 0.01','T = 0.1','T = 0.5'};


Nd = 6;

full_name = strcat(directory_name, file_name, num2str(Number{1}), file_ending);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name = info.Datasets(1).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minP = -100.0;
maxP = 100.0;
dp = (maxP - minP)/Np;
P(1:Np) = 0;
for i = 1:Np,
    P(i) = minP + (i-0.5)*dp;
end;

dt = 4000*0.5*0.1;

figure(3);
hold on;
meanP(1:Ndata,1:Nd) = 0;
for l = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(Number{l}), file_ending);
    info = h5info(full_name);
    %Ndata = size(info.Datasets,1);
    name = info.Datasets(Ndata).Name;
    fp= hdf5read(full_name, name);
    plot(P(1:Np),fp(1:Np),Color{l});
    for k = 1:Ndata,
        name = info.Datasets(k).Name;
        fp= hdf5read(full_name, name);
        particles = 0;
        for j = 1:Np,
            px = minP + (j-0.5)*dp;
            meanP(k,l) = meanP(k,l) + px*fp(j);
            particles = particles+fp(j);
        end;
        meanP(k,l) = meanP(k,l)/particles;
    end;
end;

time(1:Ndata) = 0;
for i = 1:Ndata,
    time(i) = (i-1)*dt;
end;

set(0, 'DefaultLineLineWidth', 2);

figure(1);
hold on;
for l = 1:Nd,
    plot(time(1:Ndata),meanP(1:Ndata,l)/meanP(1,l),Color{l});
end;
legend(LegendTitle{Number{1}+1}, LegendTitle{Number{2}+1}, LegendTitle{Number{3}+1}, LegendTitle{Number{4}+1},LegendTitle{Number{5}+1}, LegendTitle{Number{6}+1},'Location','northwest');
title('Px/Px0');
xlabel('t');
ylabel('Px/Px0');
grid;
