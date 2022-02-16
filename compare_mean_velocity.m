clear;
directory_name = './output1/';
file_name = 'ParticleBinning8';
file_ending = '.h5';

Number = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
%LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30', '{\theta} = 40', '{\theta} = 50','{\theta} = 60', '{\theta} = 70', '{\theta} = 80', '{\theta} = 90'};
LegendTitle = {'{\gamma} = 1.05','{\gamma} = 1.5','{\gamma} = 3', '{\gamma} = 7','{\gamma} = 10','{\gamma} = 20'};

Nd = 6;

full_name = strcat(directory_name, file_name, num2str(Number{1}), file_ending);
info = h5info(full_name);
Ndata = size(info.Datasets,1);
name = info.Datasets(1).Name;
fp= hdf5read(full_name, name);

Np=size(fp,1);
Nx=size(fp,2);

minV = -1.0;
maxV = 1.0;
dv = (maxV - minV)/Np;

dt = 4000*0.5*0.1;

meanV(1:Ndata,1:Nd) = 0;
n(1:Ndata,1:Nd) = 0;
for l = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(Number{l}), file_ending);
    info = h5info(full_name);
    %Ndata = size(info.Datasets,1);
    fp= hdf5read(full_name, name);
    for k = 1:Ndata,
        name = info.Datasets(k).Name;
        fp= hdf5read(full_name, name);
        particles = 0;
        for i = 1:Nx,
            for j = 1:Np,
                vx = minV + (j-0.5)*dv;
                meanV(k,l) = meanV(k,l) + vx*fp(j,i);
                particles = particles+fp(j,i);
                n(k,l) = n(k,l) + fp(j,i);
            end;
        end;
        meanV(k,l) = meanV(k,l)/particles;
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
    plot(time(1:Ndata),meanV(1:Ndata,l),Color{l});
end;
legend(LegendTitle{Number{1}+1}, LegendTitle{Number{2}+1}, LegendTitle{Number{3}+1}, LegendTitle{Number{4}+1}, LegendTitle{Number{5}+1}, LegendTitle{Number{6}+1},'Location','northwest');
title('V_x');
xlabel('t');
ylabel('V_x');
grid;

figure(2);
hold on;
for l = 1:Nd,
    plot(time(1:Ndata),n(1:Ndata,l),Color{l});
end;
legend(LegendTitle{Number{1}+1}, LegendTitle{Number{2}+1}, LegendTitle{Number{3}+1}, LegendTitle{Number{4}+1}, LegendTitle{Number{5}+1}, LegendTitle{Number{6}+1},'Location','northwest');
title('N');
xlabel('t');
ylabel('N');
grid;
