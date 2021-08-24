clear;
directory_name = './output/';
file_name = 'ParticleBinning7';
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name = info.Datasets(Ndata).Name;
fp= hdf5read(full_name, name);

Color = {'red','blue','green','black','magenta', [1.0,0.6,0]};

Np=size(fp,1);
Nx=size(fp,2);

Ns=5;

Fp(1:Np,1:Ns)=0;

samplingFactor = 20;

startx(1:Ns) = 0; 
endx(1:Ns) = 0; 

startx(1) = 0/samplingFactor+1;
endx(1) = 5000/samplingFactor;

startx(2) = 7000/samplingFactor+1;
endx(2) = 12000/samplingFactor;

startx(3) = 16000/samplingFactor+1;
endx(3) = 21000/samplingFactor;

startx(4) = 25000/samplingFactor+1;
endx(4) = 30000/samplingFactor;

startx(5) = 25000/samplingFactor+1;
endx(5) = 30000/samplingFactor;


%startx(1)= 5000/samplingFactor;
%endx(1) = fix(10000/samplingFactor);
%for i =2:Ns,
%    startx(i) = endx(i-1);
%    endx(i) = startx(i)+fix(5000/samplingFactor);
%end;

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

for k=1:Ns,
    norm = 0;
    Nt = fix(Ndata*k/Ns);
    name = info.Datasets(Nt).Name;
    fp= hdf5read(full_name, name);
    for i=1:Np,
        for j=startx(k):endx(k),
            Fp(i,k)=Fp(i,k)+fp(i,j)*energy(i)*energy(i)/de(i);
            norm = norm + fp(i,j);
        end;
    end;
    for i=1:Np,
        Fp(i,k) = Fp(i,k)/norm;
    end;
end;

figure(2);
hold on;
for k = 2:Ns,
    plot(energy(1:Np),Fp(1:Np,k),'color',Color{k});
end;
title('F(E)');
xlabel('E/m_e c^2');
ylabel('F(E)E^2');
%legend('far downstream', 'downstream', 'front', 'upstream', 'far upstream', 'far far upstream','Location','northwest');
%legend('1', '2', '3', '4', '5', '6','Location','northwest');
%legend('700 {\omega}_{pi}^{-1}', '1400 {\omega}_{pi}^{-1}', '2100 {\omega}_{pi}^{-1}', '2800 {\omega}_{pi}^{-1}', '3500 {\omega}_{pi}^{-1}','Location','northwest');
legend('1400 {\omega}_{pi}^{-1}', '2100 {\omega}_{pi}^{-1}', '2800 {\omega}_{pi}^{-1}', '3500 {\omega}_{pi}^{-1}','Location','northwest');
grid;