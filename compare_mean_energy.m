clear;
directory_name = './output1/';
file_name = 'ParticleBinning2';
file_ending = '.h5';

Number = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
%LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};
LegendTitle = {'{\gamma} = 1.05', '{\gamma} = 1.5','{\gamma} = 3', '{\gamma} = 7','{\gamma} = 10','{\gamma} = 20'};

Nd = 6;

full_name = strcat(directory_name, file_name, num2str(Number{1}), file_ending);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name = info.Datasets(1).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);

minE = 0.001;
maxE = 5000.0;
factor = (maxE/minE)^(1.0/(Np-1));
energy(1:Np)=0;
energy(1) = minE;
for i = 2:Np,
    energy(i) = energy(i-1)*factor;
end;
de(1:Np)=0;
de(1) = energy(2) - energy(1);
for i = 2:Np,
    de(i) = energy(i) - energy(i-1);
end;

dt = 4000*0.5*0.1;

meanE(1:Ndata,1:Nd) = 0;
set(0, 'DefaultLineLineWidth', 2);
figure(3);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
for l = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(Number{l}), file_ending);
    info = h5info(full_name);
    %Ndata = size(info.Datasets,1);
    name = info.Datasets(Ndata).Name;
    fp= hdf5read(full_name, name);
    Fp(1:Np)=0;
    for i = 1:Np,
        Fp(i)=fp(i)/de(i);
    end;
    plot(energy(1:Np),Fp(1:Np),Color{l});
    for k = 1:Ndata,
        name = info.Datasets(k).Name;
        fp= hdf5read(full_name, name);
        particles = 0;
        for j = 1:Np,
            meanE(k,l) = meanE(k,l) + energy(j)*fp(j);
            particles = particles+fp(j);
        end;
        %meanE(k,l) = meanE(k,l)/particles;
    end;
end;

time(1:Ndata) = 0;
for i = 1:Ndata,
    time(i) = (i-1)*dt;
end;

figure(1);
hold on;
for l = 1:Nd,
    plot(time(1:Ndata),meanE(1:Ndata,l)/meanE(1,l),Color{l});
end;
legend(LegendTitle{Number{1}+1}, LegendTitle{Number{2}+1}, LegendTitle{Number{3}+1}, LegendTitle{Number{4}+1}, LegendTitle{Number{5}+1}, LegendTitle{Number{5}+1},'Location','northwest');
title('E_{kin}');
xlabel('t');
ylabel('E_{kin}/E0');
grid;
