clear;
directory_name = './output1/';
file_name = 'ParticleBinning6';
file_ending = '.h5';

Number = {1,2};
Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
%LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};
LegendTitle = {'T = 0.2', 'T = 0.00001','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};

Nd = 2;
full_name = strcat(directory_name, file_name, num2str(Number{1}), file_ending);

info = h5info(full_name);
Ndata(1:Nd) = 0;
Ndata(1) = 36;
Ndata(2) = 10;
name = info.Datasets(Ndata(1)).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minE(1:Nd) = 0;
minE(1) = 0.001;
minE(2) = 0.001;
maxE(1:Nd) = 0;
maxE(1)= 5000;
maxE(2) = 1000;

mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;

m = me;

energy(1:Np,1:Nd) = 0;
de(1:Np,1:Nd) = 0;
energy(1,1) = minE(1);
energy(1,2) = minE(2);
for j = 1:Nd,
    factor = (maxE(j)/minE(j))^(1.0/(Np-1));
    for i = 2:Np,
    	energy(i,j) = energy(i-1,j)*factor;
    end;
    de(1,j) = energy(2,j) - energy(1,j);
    for i = 2:Np,
        de(i,j) = energy(i,j) - energy(i-1,j);
    end;
end;

startx(1:Nd) = 0;
endx(1:Nd) = 0;
Fp(1:Nd,1:Np)=0;

samplingFactor = 20;

for i = 1:Nd,
    startx(i) = fix(23000/samplingFactor)+1;
    endx(i) = fix(28000/samplingFactor);
end;
startx(1) = fix(1000/samplingFactor)+1;
endx(1) = fix(65000/samplingFactor);
startx(2) = fix(500/samplingFactor)+1;
endx(2) = fix(27000/samplingFactor);

for k = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(Number{k}), file_ending);
    %faaake
    info = h5info(full_name);
    %Ndata = size(info.Datasets,1);
    name = info.Datasets(Ndata(k)).Name;
    fp = hdf5read(full_name, name);
    for i=1:Np,
        for j=startx(k):endx(k),
            %Fp(k,i)=Fp(k, i)+fp(i,j)*(energy(i)+1)*(energy(i)+1)/de(i);
            Fp(k,i)=Fp(k, i)+fp(i,j)/de(i,k);
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

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F_e(E)');
xlabel ('E/{m_e c^2}');
ylabel ('F_e(E)');
for j=1:Nd,
    plot (energy(1:Np,j)+m/me,Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2},'Location','northwest');
grid ;