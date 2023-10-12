clear;
directory_name = './output_theta80_gamma0.5_sigma0.0002_mass25-400/';
file_name = 'ParticleBinning6';
file_ending = '.h5';

Suffix = {'_25','_64','_100','_144','_400'};
Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'m = 25', 'm = 64','m = 100', 'm = 144', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Nd = 4;
full_name = strcat(directory_name, file_name, Suffix{1}, file_ending);

info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 5;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minE = 0.001;
maxE = 1000;
factor = (maxE/minE)^(1.0/(Np-1));

mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;

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

startx(1:Nd) = 0;
endx(1:Nd) = 0;
Fp(1:Nd,1:Np)=0;

samplingFactor = 20;

for i = 1:Nd,
    startx(i) = fix(15000/samplingFactor)+1;
    endx(i) = fix(25000/samplingFactor);
end;
%startx(2) = fix(23000/samplingFactor)+1;
%endx(2) = fix(28000/samplingFactor);

for k = 1:Nd,
    full_name = strcat(directory_name, file_name, Suffix{k}, file_ending);
    %faaake
    info = h5info(full_name);
    %Ndata = size(info.Datasets,1);
    name = info.Datasets(Ndata).Name;
    fp = hdf5read(full_name, name);
    for i=1:Np,
        for j=startx(k):endx(k),
            %Fp(k,i)=Fp(k, i)+fp(i,j)*(energy(i)+1)*(energy(i)+1)/de(i);
            Fp(k,i)=Fp(k, i)+fp(i,j)/de(i);
        end;
    end;
end;

for k = 1:Nd,
    norm = 1.0;
    normp = 0.0;
    for i = 1:Np,
        normp = normp + Fp(k,i)*de(i)*me/m;
    end;

    for i = 1:Np,
        Fp(k,i) = Fp(k,i)*norm/normp;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1.5);

for j = 1:Nd,
    norm = 0;
    for i = 1:Np,
        norm = norm + Fp(j,i)*de(i);
    end;
    for i = 1:Np,
        Fp(j,i) = Fp(j,i)/norm;
    end;
end;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_e(E)');
xlabel ('E/{m_e c^2}');
ylabel ('F_e(E)');
for j=1:Nd,
    plot (energy(1:Np)+m/me,Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4},'Location','northwest');
grid ;

output(1:180,1:Nd+1)=0;
for i = 1:180, 
    output(i,1) = log10(energy(i) + 1);
    for j = 1:Nd,
        output(i, j +1) = Fp(j,i);
    end;
end;
dlmwrite('compareelectrons.dat',output,'delimiter',' ');