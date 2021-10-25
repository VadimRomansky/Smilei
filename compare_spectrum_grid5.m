clear;
directory_name = './output_theta0-90_gamma1.5sigma0.004/';
file_name = 'ParticleBinning6';
file_ending = '.h5';

Number = {1, 3, 4, 6, 9};
Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Nd = 5;
start = 0;
fileNumber = start;
full_name = strcat(directory_name, file_name, num2str(Number{1}), file_ending);

info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 10;
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
    startx(i) = fix(15000/samplingFactor)+1;
    endx(i) = fix(20000/samplingFactor);
end;
%startx(2) = fix(27000/samplingFactor)+1;
%endx(2) = fix(32000/samplingFactor);

%startx(5) = fix(33000/samplingFactor)+1;
%endx(5) = fix(38000/samplingFactor);

for k = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(Number{k}), file_ending);
    %faaake
    info = h5info(full_name);
    Ndata = size(info.Datasets,1);
    name = info.Datasets(Ndata).Name;
    fp = hdf5read(full_name, name);
    for i=1:Np,
        for j=startx:endx,
            %Fp(k,i)=Fp(k, i)+fp(i,j)*(energy(i)+1)*(energy(i)+1)/de(i);
            Fp(k,i)=Fp(k, i)+fp(i,j)/de(i);
        end;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1.5);

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_e(E)');
xlabel ('E/{m_e c^2}');
ylabel ('F_e(E)*E^2');
for j=1:Nd,
    plot (energy(1:Np)+1,Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{Number{1}+1}, LegendTitle{Number{2}+1}, LegendTitle{Number{3}+1}, LegendTitle{Number{4}+1}, LegendTitle{Number{5}+1},'Location','northwest');
grid ;

for j = 1:Nd,
    norm = 0;
    for i = 1:Np,
        norm = norm + Fp(j,i)*de(i);
    end;
    for i = 1:Np,
        Fp(j,i) = Fp(j,i)/norm;
    end;
end;

output(1:180,1:Nd+1)=0;
for i = 1:180, 
    output(i,1) = log10(energy(i) + 1);
    for j = 1:Nd,
        output(i, j +1) = Fp(j,i);
    end;
end;
dlmwrite('compareelectrons.dat',output,'delimiter',' ');