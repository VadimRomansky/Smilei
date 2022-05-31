clear;
directory_name = './output/';
file_name = 'ParticleBinning';
file_ending = '.h5';

Number = {6,7,8};
Color = {'blue','red','magenta'};
LegendTitle = {'electrons', 'protons','positrons'};

Nd = 3;
full_name = strcat(directory_name, file_name, num2str(Number{1}), file_ending);

info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 30;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minE(1:Nd) = 0;
minE(1) = 0.001;
minE(2) = 0.1;
minE(3) = 0.001;
maxE(1:Nd) = 0;
maxE(1)= 50000;
maxE(2) = 50000;
maxE(3) = 50000;

mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;

m(1:Nd) = [me, mp, me];

energy(1:Np,1:Nd) = 0;
de(1:Np,1:Nd) = 0;
energy(1,1) = minE(1);
energy(1,2) = minE(2);
energy(1,3) = minE(3);
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
    startx(i) = fix(300000/samplingFactor)+1;
    endx(i) = fix(400000/samplingFactor);
end;
%startx(2) = fix(23000/samplingFactor)+1;
%endx(2) = fix(28000/samplingFactor);

for k = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(Number{k}), file_ending);
    %faaake
    info = h5info(full_name);
    %Ndata = size(info.Datasets,1);
    name = info.Datasets(Ndata).Name;
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
        normp = normp + Fp(k,i)*de(i,k)*me/m(k);
    end;

    for i = 1:Np,
        Fp(k,i) = Fp(k,i)*norm/normp;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1.5);

% for j = 1:Nd,
%     norm = 0;
%     for i = 1:Np,
%         norm = norm + Fp(j,i)*de(i,j);
%     end;
%     for i = 1:Np,
%         Fp(j,i) = Fp(j,i)/norm;
%     end;
% end;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title ('F({\gamma})');
xlabel ('{\gamma}*m/m_e');
ylabel ('F({\gamma})');
for j=1:Nd,
    %plot (energy(1:Np,j)+m(j)/me,Fp(j, 1:Np),'color',Color{j});
    plot (energy(1:Np,j) + m(j)/mp,Fp(j, 1:Np),'color',Color{j});
end;

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3},'Location','northwest');
grid ;

output(1:180,1:Nd+1)=0;
for i = 1:180, 
    output(i,1) = log10(energy(i) + 1);
    for j = 1:Nd,
        output(i, j +1) = Fp(j,i);
    end;
end;
dlmwrite('compareelectrons.dat',output,'delimiter',' ');