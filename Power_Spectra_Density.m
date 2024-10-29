clear all, close all, clc

numero_archivos=12;
subarchivos=2;
figura=1;
lgd={'Bin 10','Bin 9','Bin 6','Bin 3','Bin 1'};

for bloque=0:numero_archivos

for subbloque=1:subarchivos


fpath = '/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS'; 

nombre='VELOCIDADES_';
fnumber=int2str(bloque); %Number of files 
barra='_';
codigo=int2str(subbloque);
archivo='.mat';
fname=strcat(nombre,fnumber,barra,codigo,archivo)

%load([ fpath '/' fname ]);

Variables={'Burst_VelEast','Burst_VelNorth','Burst_VelUp'};
Data2=load([ fpath '/' fname ],Variables{:});

Data2.fs=8;%es la frecuencia de muestreo del equipo
dt=1/Data2.fs;

%Inicio=57601;
%Final=129601;

% con estos datos subdivido la serie completa

N_muestras=[8*30*1];% frecuencia de muestreo * segundos * minutos
N_max=length(Data2.Burst_VelUp); % La longitud máxima para coger los datos



%bin=[1,3,6,9,10]; %elijo el bin para calcular
bin=[10,9,6,3,1]; %elijo el bin para calcular


%% VELOCIDADES ESTE NORTE ARRIBA

for i=1:length(bin)
vel=Data2.Burst_VelUp(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],Data2.fs); %obtengo el periodograma de la velocidad
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

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/COMPLETO/VEL_UP/Vel_Up_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()

for i=1:length(bin)
vel=Data2.Burst_VelNorth(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],Data2.fs); %obtengo el periodograma de la velocidad
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

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/COMPLETO/VEL_NORTH/Vel_North_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()

for i=1:length(bin)
    
vel=Data2.Burst_VelEast(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],Data2.fs); %obtengo el periodograma de la velocidad
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

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/COMPLETO/VEL_EAST/Vel_East_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()
%Hago la media después de calcular el periodograma

figure (figura)
for b=1:length(bin)
    vel=Data2.Burst_VelUp(:,bin(b)); %CAMBIAR VELOCIDAD
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


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/VEL_UP/Vel_Up_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

figure (figura)
for b=1:length(bin)
    vel=Data2.Burst_VelNorth(:,bin(b)); %CAMBIAR VELOCIDAD
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


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/VEL_NORTH/Vel_North_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

figure (figura)
for b=1:length(bin)
    vel=Data2.Burst_VelEast(:,bin(b)); %CAMBIAR VELOCIDAD
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


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/VEL_EAST/Vel_East_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

%clearvars -except figura l n numero_archivos subarchivos

end

end



%% VELOCIDADES X Y Z

for i=1:length(bin)
vel=Data2.Burst_VelZ(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],Data2.fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
figure (figura)
loglog(f_total,pxx_total)
xlabel(['Frequency (Hz)'])
ylabel('S(m^2s^-^2Hz^-^1)')
title('Power Spectra Density Velocity Z')
hold on
end
legend('Bin 1','Bin 3','Bin 6','Bin 9','Bin 10')
hold off

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/COMPLETO/VEL_Z/Vel_Z_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()

for i=1:length(bin)
vel=Data2.Burst_VelY(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],Data2.fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
figure (figura)
loglog(f_total,pxx_total)
xlabel(['Frequency (Hz)'])
ylabel('S(m^2s^-^2Hz^-^1)')
title('Power Spectra Density Velocity Y')
hold on
end
legend('Bin 1','Bin 3','Bin 6','Bin 9','Bin 10')
hold off

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/COMPLETO/VEL_Y/Vel_Y_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()

for i=1:length(bin)
vel=Data2.Burst_VelX(:,bin(i)); % CAMBIAR VELOCIDAD

[pxx_total,f_total]=periodogram(vel,[],[],Data2.fs); %obtengo el periodograma de la velocidad
% estoy aplicando una ventana rectangular debido a que es la que viene por
% defecto
figure (figura)
loglog(f_total,pxx_total)
xlabel(['Frequency (Hz)'])
ylabel('S(m^2s^-^2Hz^-^1)')
title('Power Spectra Density Velocity X')
hold on
end
legend('Bin 1','Bin 3','Bin 6','Bin 9','Bin 10')
hold off

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/COMPLETO/VEL_X/Vel_X_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel pxx_total f_total

close()
%Hago la media después de calcular el periodograma

figure (figura)
for b=1:length(bin)
    vel=Data2.Burst_VelZ(:,bin(b)); %CAMBIAR VELOCIDAD
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
    title('Power Spectra Density Velocity Z')
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
legend('Bin 1','Bin 3','Bin 6','Bin 9','Bin 10')
hold off


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/VEL_Z/Vel_Z_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

figure (figura)
for b=1:length(bin)
    vel=Data2.Burst_VelY(:,bin(b)); %CAMBIAR VELOCIDAD
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
    title('Power Spectra Density Velocity Y')
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
legend('Bin 1','Bin 3','Bin 6','Bin 9','Bin 10')
hold off


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/VEL_Y/Vel_Y_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

figure (figura)
for b=1:length(bin)
    vel=Data2.Burst_VelX(:,bin(b)); %CAMBIAR VELOCIDAD
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
    title('Power Spectra Density Velocity X')
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
legend('Bin 1','Bin 3','Bin 6','Bin 9','Bin 10')
hold off


nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/ESPECTROS/VEL_X/Vel_X_',fnumber,barra,codigo);

saveas(figure(figura),nombre_figura,'fig')

saveas(figure(figura),nombre_figura,'jpeg')

figura=figura+1;

clear vel vel_subserie S F pxx f frecuencia periodograma posicion_total

close()

%clearvars -except figura l n numero_archivos subarchivos

end

end


% 
figure(3)
hago la media antes de calcular el periodograma
figure(3)
for l=1:length(N_muestras)
    N_min=N_muestras(1,l);
    subseries=N_max/N_min;
    posicion_total=1;
    transformada=[];
    frecuencia=[];
    vel_subserie=[];
    T=[];
    V=[];
        for j=1:subseries
                for i=1:N_min
                    vel_subserie(i,j)=vel(posicion_total,1);
                    posicion_total=posicion_total+1;
                end
            transformada(:,j)=fft(vel_subserie(:,j));
        end
    T(:,l)=mean(transformada,2);%hago la media de cada fila 
    V(:,l)=ifft(T(:,l)); %hago la transformada inversa
    [pxx,f]=periodogram(V(:,l));
    F(:,l)=f;
    S(:,l)=pxx;
    S(1,l)=NaN; %borro el primer valor de la serie porque no sirve
    F(1,l)=NaN;

    loglog(F(:,l),S(:,l))
    xlabel(['Frequency (Hz)'])
    ylabel('S(m^2s^-^2Hz^-^1)')
    title('Periodograma Subseries FFT')
    hold on
end
legend('25','50','100')
hold off

clear S F T V
% para sacar las frecuencias
% Ns=length(vel);
% 
% ks=1:1:Ns;
% 
% dfs=1/(Ns*dt);
% 
% fks=ks*dfs;
% 
% % el módulo de la transformada de Fourier al cuadrado
% 
%  S_l(:,j)=dt*(abs(XHs1_ft(:,j))).^2;

figure (4)
for l=1:length(N_muestras)
    N_min=N_muestras(1,l);
    subseries=N_max/N_min;
    posicion_total=1;            
    vel_subserie=[];
    transformada=[];
    frecuencia=[];
    FFT=[];
    S=[];
    F=[];
        for j=1:subseries
            ks=[];
            dfs=[];
            fks=[];
                for i=1:N_min
                    vel_subserie(i,j)=vel(posicion_total,1);
                    posicion_total=posicion_total+1;
                end
            transformada(:,j)=fft(vel_subserie(:,j));
            Ns=length(vel_subserie(:,j));
            ks=1:1:Ns;
            dfs=1/(Ns*dt);            
            fks=ks*dfs;
            frecuencia(:,j)=fks;
        end
    FFT=mean(transformada,2);
    S=2/(dt*N_muestras(l)).*(abs(FFT)).^2;
    F=mean(frecuencia,2);
    %S(1,l)=NaN; 
    %F(1,l)=NaN;
    
    loglog(F(2:round(length(FFT)/2),1),S(2:round(length(FFT)/2),1))
    xlabel(['Frequency (Hz)'])
    ylabel('S(m^2s^-^2Hz^-^1)')
    title('Periodograma Subseries FFT 2')
    hold on
end
legend('25','50','100')
hold off

    
