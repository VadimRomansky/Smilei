clear;
directory_name = './output/';
file_name = 'ParticleBinning6';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
%Ndata = 10;
%Ndata = 100;
%Ndata = 1;
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Color = {'red','blue','green','black','magenta', [1.0,0.6,0]};

me=1;
mp=100;
m=me;

Np=size(fp,1);
Nx=size(fp,2);

Ns=6;

minEe = 0.00001;
maxEe = 1000;
minEp = 0.1;
maxEp = 5000;
minE = minEe;
maxE = maxEe;
factor = (maxE/minE)^(1.0/(Np-1));

Fp(1:Np,1:Ns)=0;

samplingFactor = 20;

startx(1:Ns) = 0; 
endx(1:Ns) = 0; 

startx(1) = 100/samplingFactor+1;
endx(1) = startx(1) + 1000/samplingFactor;

startx(2) = 1000/samplingFactor;
endx(2) = startx(2) + 1000/samplingFactor;

startx(3) = 2000/samplingFactor;
endx(3) = startx(3) + 5000/samplingFactor;

startx(4) = 10000/samplingFactor;
endx(4) = startx(4) + 5000/samplingFactor;

startx(5) = 15000/samplingFactor;
endx(5) = startx(5) + 5000/samplingFactor;

startx(6) = 50000/samplingFactor;
endx(6) = startx(6) + 2000/samplingFactor;

%startx(1)= 5000/samplingFactor;
%endx(1) = fix(10000/samplingFactor);
%for i =2:Ns,
%    startx(i) = endx(i-1);
%    endx(i) = startx(i)+fix(5000/samplingFactor);
%end;

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

for k=1:Ns,
    norm = 0;
    for i=1:Np,
        for j=startx(k):endx(k),
            Fp(i,k)=Fp(i,k)+fp(i,j)/de(i);
            norm = norm + fp(i,j);
        end;
    end;
    %comment for compare
    %for i =1:Np,
    %    Fp(i,k) = Fp(i,k)/norm;
    %end;
end;

figure(1);
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
xlim([0.00001 1000]);
%ylim([10^-10 10]);
for k = 1:Ns,
    plot(energy(1:Np),Fp(1:Np,k),'color',Color{k},'LineWidth',2);
end;
title('F(E)');
xlabel('E/me c^2');
ylabel('F(E)');
%legend('1','2','3','4','5','6');
legend('far far downstream', ' far downstream', 'downstream', 'front', 'upstream', 'far upstream','Location','northeast');
grid;