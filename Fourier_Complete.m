load('Arriba.mat')
load('Este.mat')
load('Norte.mat')

figura=1;

lgd={'Bin 10','Bin 1','Bin 5'};
bin = [10,1,5];
size_letra=16;

fs=8; 

%FOURIER COMPLETO UP

figure (figura)
for i=1:length(bin)
vel=Vel_Up(:,bin(i)); % CAMBIAR VELOCIDAD
[pxx_total,f_total]=periodogram(vel,[],[],fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
loglog(f_total,pxx_total)
hold on
end

hold off
xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/GRAFICOS/fourier_up');

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel pxx_total f_total

close()

% FOURIER COMPLETO ESTE

figure (figura)
for i=1:length(bin)
vel=Vel_East(:,bin(i)); % CAMBIAR VELOCIDAD
[pxx_total,f_total]=periodogram(vel,[],[],fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
loglog(f_total,pxx_total)
hold on
end

hold off
xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/GRAFICOS/fourier_east');

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel pxx_total f_total

close()

% FOURIER COMPLETO NORTE

figure (figura)
for i=1:length(bin)
vel=Vel_North(:,bin(i)); % CAMBIAR VELOCIDAD
[pxx_total,f_total]=periodogram(vel,[],[],fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
loglog(f_total,pxx_total)
hold on
end

hold off
xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/GRAFICOS/fourier_north');

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel pxx_total f_total

close()

%% FOURIER DE 3 AL 10 DE ENERO DE 2024

load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

size_letra=16;

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

%N_muestras=[8*60*1];% frecuencia de muestreo * segundos * minutos
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
hold on
end
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold on
end
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold on
end
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel pxx_total f_total

close()

% FOURIER TRAMOS DE 30 MINUTOS 

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')




nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')



nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')



nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()


% FOURIER TRAMOS DE 10 MINUTOS

% Fourier segmentado

%tipo_serie='Corta';
tipo_serie='100_SEGUNDOS';


figura=1;

lgd={'Bin 10','Bin 1','Bin 5'};
%bin=[1,3,6,9,10]; %elijo el bin para calcular
bin=[10,1,5]; %elijo el bin para calcular


fs=8;%es la frecuencia de muestreo del equipo
dt=1/fs;


% con estos datos subdivido la serie completa

N_muestras=[8*100*1];% frecuencia de muestreo * segundos * minutos
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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')




nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')
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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')




nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')



nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

%% FOURIER DE 16 AL 23 DE ENERO DE 2024

load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

size_letra=16;

velocidad_norte=Vel_North_16_23;
velocidad_up=Vel_Up_16_23;
velocidad_este=Vel_East_16_23;

dia_1=16;
dia_2=23;

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
hold on
end
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')
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
hold on
end
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold on
end
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')


nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel pxx_total f_total

close()



% FOURIER TRAMOS DE 30 MINUTOS 

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')




nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')



nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')



nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()


% FOURIER TRAMOS DE 10 MINUTOS

% Fourier segmentado

%tipo_serie='Corta';
tipo_serie='100_SEGUNDOS';


figura=1;

lgd={'Bin 10','Bin 1','Bin 5'};
%bin=[1,3,6,9,10]; %elijo el bin para calcular
bin=[10,1,5]; %elijo el bin para calcular


fs=8;%es la frecuencia de muestreo del equipo
dt=1/fs;


% con estos datos subdivido la serie completa

N_muestras=[8*100*1];% frecuencia de muestreo * segundos * minutos
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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')




nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_Up_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')
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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')




nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_North_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

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
hold off

xlabel(['Frequency (Hz)'],'FontSize', size_letra, 'FontWeight', 'bold')
ylabel('S(m^2s^-^2Hz^-^1)','FontSize', size_letra, 'FontWeight', 'bold')
%title('Power Spectra Density Velocity Up','FontSize', size_letra, 'FontWeight', 'bold')
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes
legend(lgd{:},'FontSize', size_letra, 'FontWeight', 'bold')



nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Vel_East_',tipo_serie,'_',dia_1,'_',dia_2);

%saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

