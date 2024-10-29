clear all, close all, clc

load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_3_10;
velocidad_up=Vel_Up_3_10;
velocidad_este=Vel_East_3_10;

dia_1=3;
dia_2=10;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

figura=1;

lgd={'Bin 10','Bin 1','Bin 5'};
%bin=[1,3,6,9,10]; %elijo el bin para calcular
bin=[10,1,5]; %elijo el bin para calcular


fs=8;%es la frecuencia de muestreo del equipo
dt=1/fs;


% con estos datos subdivido la serie completa

N_muestras=[8*60*1];% frecuencia de muestreo * segundos * minutos
N_max=length(velocidad_norte); % La longitud máxima para coger los datos

%strcat('Vel_Up_',dia_1,'_',dia_2)




% VELOCIDADES ESTE NORTE ARRIBA

for i=1:length(bin)
vel=velocidad_up(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
figure (figura)
loglog(f_total,pxx_total)
xlabel(['Frequency (Hz)'])
ylabel('S(m^2s^-^2Hz^-^1)')
title('Power Spectra Density Velocity Up')
hold on
end
legend(lgd{:})
hold off

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',dia_1,'_',dia_2);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()

for i=1:length(bin)
vel=velocidad_norte(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
figure (figura)
loglog(f_total,pxx_total)
xlabel(['Frequency (Hz)'])
ylabel('S(m^2s^-^2Hz^-^1)')
title('Power Spectra Density Velocity North')
hold on
end
legend(lgd{:})
hold off

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',dia_1,'_',dia_2);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()

for i=1:length(bin)
    
vel=velocidad_este(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
figure (figura)
loglog(f_total,pxx_total)
xlabel(['Frequency (Hz)'])
ylabel('S(m^2s^-^2Hz^-^1)')
title('Power Spectra Density Velocity East')
hold on
end
legend(lgd{:})
hold off

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',dia_1,'_',dia_2);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()
%Hago la media después de calcular el periodograma
%% Fourier segmentado

clear all, close all, clc

load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_16_23;
velocidad_up=Vel_Up_16_23;
velocidad_este=Vel_East_16_23;

dia_1=16;
dia_2=23;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

%tipo_serie='Corta';
tipo_serie='30_MINUTOS';


figura=1;

lgd={'Bin 10','Bin 1','Bin 5'};
%bin=[1,3,6,9,10]; %elijo el bin para calcular
bin=[10,1,5]; %elijo el bin para calcular


fs=8;%es la frecuencia de muestreo del equipo
dt=1/fs;


% con estos datos subdivido la serie completa

N_muestras=[8*60*30];% frecuencia de muestreo * segundos * minutos
N_max=length(velocidad_norte); % La longitud máxima para coger los datos

%strcat('Vel_Up_',dia_1,'_',dia_2)


figure (figura)
for b=1:length(bin)
    vel=velocidad_up(:,bin(b)); %CAMBIAR VELOCIDAD
for l=1:length(N_muestras)
    N_min=N_muestras(1,l);
    subseries=N_max/N_min;
    posicion_total=1;            
    vel_subserie=[];
        for j=1:subseries
                for i=1:N_min
                    vel_subserie(i,j)=vel(posicion_total,1);
                    posicion_total=posicion_total+1;
                end
            [pxx,f]=periodogram(vel_subserie(:,j));
            frecuencia(:,j)=f;
            periodograma(:,j)=pxx;
        end
    S(:,l)=mean(periodograma,2);%hago la media de cada fila y las guardo en cada columna para cada subserie
    F(:,l)=mean(frecuencia,2);
    S(1,l)=NaN; %borro el primer valor de la serie porque no sirve
    F(1,l)=NaN;
    S(end,l)=NaN;
    F(end,l)=NaN;

    
    loglog(F(:,l),S(:,l))
    xlabel(['Frequency (Hz)'])
    ylabel('S(m^2s^-^2Hz^-^1)')
    title('Power Spectra Density Velocity Up')
    hold on

end
end

m=-5/3;
x1=3;
x2=8;
y1=0.01;
y2=10^(m*log10(x2/x1)+log10(y1));
x=[x1,x2];
y=[y1,y2];
plot(x,y)
legend(lgd{:})
hold off


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',tipo_serie,'_',dia_1,'_',dia_2);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

figure (figura)
for b=1:length(bin)
    vel=velocidad_norte(:,bin(b)); %CAMBIAR VELOCIDAD
for l=1:length(N_muestras)
    N_min=N_muestras(1,l);
    subseries=N_max/N_min;
    posicion_total=1;            
    vel_subserie=[];
        for j=1:subseries
                for i=1:N_min
                    vel_subserie(i,j)=vel(posicion_total,1);
                    posicion_total=posicion_total+1;
                end
            [pxx,f]=periodogram(vel_subserie(:,j));
            frecuencia(:,j)=f;
            periodograma(:,j)=pxx;
        end
    S(:,l)=mean(periodograma,2);%hago la media de cada fila y las guardo en cada columna para cada subserie
    F(:,l)=mean(frecuencia,2);
    S(1,l)=NaN; %borro el primer valor de la serie porque no sirve
    F(1,l)=NaN;
    S(end,l)=NaN;
    F(end,l)=NaN;

    
    loglog(F(:,l),S(:,l))
    xlabel(['Frequency (Hz)'])
    ylabel('S(m^2s^-^2Hz^-^1)')
    title('Power Spectra Density Velocity North')
    hold on

end
end

m=-5/3;
x1=3;
x2=8;
y1=0.01;
y2=10^(m*log10(x2/x1)+log10(y1));
x=[x1,x2];
y=[y1,y2];
plot(x,y)
legend(lgd{:})
hold off


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',tipo_serie,'_',dia_1,'_',dia_2);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

figure (figura)
for b=1:length(bin)
    vel=velocidad_este(:,bin(b)); %CAMBIAR VELOCIDAD
for l=1:length(N_muestras)
    N_min=N_muestras(1,l);
    subseries=N_max/N_min;
    posicion_total=1;            
    vel_subserie=[];
        for j=1:subseries
                for i=1:N_min
                    vel_subserie(i,j)=vel(posicion_total,1);
                    posicion_total=posicion_total+1;
                end
            [pxx,f]=periodogram(vel_subserie(:,j));
            frecuencia(:,j)=f;
            periodograma(:,j)=pxx;
        end
    S(:,l)=mean(periodograma,2);%hago la media de cada fila y las guardo en cada columna para cada subserie
    F(:,l)=mean(frecuencia,2);
    S(1,l)=NaN; %borro el primer valor de la serie porque no sirve
    F(1,l)=NaN;
    S(end,l)=NaN;
    F(end,l)=NaN;

    
    loglog(F(:,l),S(:,l))
    xlabel(['Frequency (Hz)'])
    ylabel('S(m^2s^-^2Hz^-^1)')
    title('Power Spectra Density Velocity East')
    hold on

end
end

m=-5/3;
x1=3;
x2=8;
y1=0.01;
y2=10^(m*log10(x2/x1)+log10(y1));
x=[x1,x2];
y=[y1,y2];
plot(x,y)
legend(lgd{:})
hold off


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',tipo_serie,'_',dia_1,'_',dia_2);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

%% WAVELETS

clear all, close all, clc

%load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_16_23;
velocidad_up=Vel_Up_16_23;
velocidad_este=Vel_East_16_23;
longitud=length(velocidad_norte);

clear Vel_North_16_23 Vel_Up_16_23 Vel_East_16_23

dia_1=16;
dia_2=23;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

bin=[10,1,5];

fs=8;

figura=1;


for i=1:length(bin)
    
vel=velocidad_este(:,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)

c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(),nombre_figura,'fig')

saveas(figure(),nombre_figura,'jpeg')

close()

end

for i=1:length(bin)
    
vel=velocidad_up(:,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)

c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(),nombre_figura,'fig')

saveas(figure(),nombre_figura,'jpeg')

close()

end

for i=1:length(bin)
    
vel=velocidad_norte(:,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)

c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(),nombre_figura,'fig')

saveas(figure(),nombre_figura,'jpeg')

close()


end


for i=1:length(bin)
    
vel=velocidad_este(:,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)

c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(),nombre_figura,'fig')

saveas(figure(),nombre_figura,'jpeg')

close()

end

for i=1:length(bin)
    
vel=velocidad_up(:,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)

c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(),nombre_figura,'fig')

saveas(figure(),nombre_figura,'jpeg')

close()


end

for i=1:length(bin)
    
vel=velocidad_norte(:,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)

c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_WSST_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(),nombre_figura,'fig')

saveas(figure(),nombre_figura,'jpeg')

close()

end

%% WAVELETS PRUEBA 2 EN FUNCIONAMIENTO

%module load matlab
%matlab -nodisplay -nodesktop

clear all, close all, clc

%load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_3_10;
velocidad_up=Vel_Up_3_10;
velocidad_este=Vel_East_3_10;

longitud=length(velocidad_norte);


clear Vel_North_3_10 Vel_Up_3_10 Vel_East_3_10

dia_1=3;
dia_2=10;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

bin=[10,1,5];

fs=8;

for i=1:length(bin)
    
vel=velocidad_este(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/ESTE/Vel_East_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel


vel=velocidad_este(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/ESTE/Vel_East_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel


end





for i=1:length(bin)
    
vel=velocidad_up(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/UP/Vel_Up_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel


vel=velocidad_up(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/UP/Vel_Up_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel


end



for i=1:length(bin)
    
vel=velocidad_norte(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/NORTE/Vel_North_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_norte(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

cwt(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/NORTE/Vel_North_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

end

%% wsst
 
clear all, close all, clc

%load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_16_23;
velocidad_up=Vel_Up_16_23;
velocidad_este=Vel_East_16_23;
longitud=length(velocidad_norte);

clear Vel_North_16_23 Vel_Up_16_23 Vel_East_16_23

dia_1=16;
dia_2=23;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

bin=[10,1,5];
tramos=7;
fs=8;

for i=1:length(bin)
for j=0:tramos-1   

vel=velocidad_norte((1+j*longitud/tramos):((j+1)*longitud/tramos),bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_WSST_',int2str(j),'_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_este((1+j*longitud/tramos):((j+1)*longitud/tramos),bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_',int2str(j),'_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_up((1+j*longitud/tramos):((j+1)*longitud/tramos),bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_',int2str(j),'_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel


end
end


%% TROZO DE PRUEBA
for i=1:length(bin)
vel=velocidad_norte(longitud/4:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_WSST_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel


vel=velocidad_norte(longitud/2:3*longitud/4,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_WSST_3_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_norte(3*longitud/4:3*longitud,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_WSST_4_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

end

for i=1:length(bin)
    
vel=velocidad_este(1:longitud/4,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_este(longitud/4:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')
clear vel

vel=velocidad_este(longitud/2:3*longitud/4,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_3_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')
clear vel

vel=velocidad_este(3*longitud/4:longitud,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_4_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')
clear vel

end



for i=1:length(bin)
    
vel=velocidad_up(1:longitud/4,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_up(longitud/4:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_up(longitud/2:3*longitud/4,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_3_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

vel=velocidad_up(3*longitud/4:longitud,bin(i)); % CAMBIAR VELOCIDAD

wsst(vel,fs)
fig=gcf;
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_4_',dia_1,'_',dia_2,'_',int2str(bin(i)));

%saveas(fig,nombre_figura,'fig')

saveas(fig,nombre_figura,'jpeg')

close()

clear vel

end

%[sst,f]=wsst(vel,fs);

%[wt,f]=cwt(velocidad_up(:,10),fs);





%% TRAMOS MÁS CORTOS WAVELETS

for i=1:length(bin)
    
vel=velocidad_este(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

[wt,f]=cwt(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(wt))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear wt f tms 

figura=figura+1;

vel=velocidad_este(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

[wt,f]=cwt(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(wt))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear wt f tms 

figura=figura+1;

end

for i=1:length(bin)
    
vel=velocidad_este(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

[sst,f]=wsst(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(sst))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear sst f tms

figura=figura+1;

vel=velocidad_este(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

[sst,f]=wsst(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(sst))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_East_WSST_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear sst f tms

figura=figura+1;

end






for i=1:length(bin)
    
vel=velocidad_up(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

[wt,f]=cwt(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(wt))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear wt f tms

figura=figura+1;

vel=velocidad_up(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

[wt,f]=cwt(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(wt))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear wt f tms

figura=figura+1;


end



for i=1:length(bin)
    
vel=velocidad_up(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

[sst,f]=wsst(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(sst))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear sst f tms


figura=figura+1;

vel=velocidad_up(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

[sst,f]=wsst(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(sst))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_Up_WSST_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear sst f tms

figura=figura+1;

end


for i=1:length(bin)
    
vel=velocidad_norte(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

[wt,f]=cwt(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(wt))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear wt f tms

figura=figura+1;

vel=velocidad_norte(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

[wt,f]=cwt(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(wt))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear wt f tms

figura=figura+1;

end



for i=1:length(bin)
    
vel=velocidad_norte(1:longitud/2,bin(i)); % CAMBIAR VELOCIDAD

[sst,f]=wsst(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(sst))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_WSST_1_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear sst f tms

figura=figura+1;

vel=velocidad_norte(longitud/2:longitud,bin(i)); % CAMBIAR VELOCIDAD

[sst,f]=wsst(vel,fs);

figure (figura)
tms = (0:numel(vel)-1)/fs;
surface(tms,f,abs(sst))
axis tight
shading flat
c=colorbar;
c.Label.String='S(m^2s^-^2Hz^-^1)';
title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/WAVELETS/Vel_North_WSST_2_',dia_1,'_',dia_2,'_',int2str(bin(i)));

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

close()

clear sst f tms

figura=figura+1;

end

%[sst,f]=wsst(vel,fs);

%[wt,f]=cwt(velocidad_up(:,10),fs);

%% COHERENCIA FOURIER Y WAVELETS

clear all, close all, clc

load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_3_10;
velocidad_up=Vel_Up_3_10;
velocidad_este=Vel_East_3_10;

dia_1=3;
dia_2=10;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

figura=1;

lgd={'Bin 10','Bin 1','Bin 5'};
%bin=[1,3,6,9,10]; %elijo el bin para calcular
bin=[10,1,5]; %elijo el bin para calcular


fs=8;%es la frecuencia de muestreo del equipo
dt=1/fs;
