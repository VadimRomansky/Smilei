clear;
directory_name = './output/';
file_name = 'ParticleBinning7';
%need distribution in the rest frame?
file_number = '.h5';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
Ndata = 6;
name = info.Datasets(Ndata).Name;

fp = hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minEe = 0.001;
maxEe = 1000;
minEp = 0.1;
maxEp = 5000;
minE = minEp;
maxE = maxEp;
factor = (maxE/minE)^(1.0/(Np-1));

me = 1;
mp = 100;
m = mp;

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

Fp(1:Np)=0;

samplingFactor = 20;

startx = fix(5000/samplingFactor)+1;
endx = fix(15000/samplingFactor);

for i=1:Np,
    for j=startx:endx,
        Fp(i)=Fp(i)+fp(i,j)/de(i);
    end;
end;

%int 2 Pi cos x cos x sin x dx = 4 Pi/3
pconst = 4*pi/3;

P = 0;
E = 0;
Ek = 0;

for i=1:Np,
    en = me*energy(i) + m;
    gam = en/m;
    p = sqrt(en*en - m*m);
    v = p/en;
    P = P + Fp(i)*pconst*p*v*de(i)/(4*pi);%because F(p) already has 4 pi
    E = E + Fp(i)*en*de(i);
    Ek = Ek + Fp(i)*(en - m)*de(i);
end;

index = (P+E)/E;
index1 = (P+Ek)/Ek;


