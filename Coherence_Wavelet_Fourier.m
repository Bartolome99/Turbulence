clear all, close all, clc

load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_3_10;
velocidad_up=Vel_Up_3_10;
velocidad_este=Vel_East_3_10;

clear Vel_North_3_10 Vel_Up_3_10 Vel_East_3_10

dia_1=3;
dia_2=10;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

figura=1;

lgd={'Coherence Layer 1-10','Coherence Layer 1-5','Coherence Layer 5-10'};
%bin=[1,3,6,9,10]; %elijo el bin para calcular
bin=[10,1,5]; %elijo el bin para calcular

s=10;
i=1;
m=5;

fs=8;%es la frecuencia de muestreo del equipo

filtro=8*60*30; %frecuencia*segundos*minutos

longitud=length(velocidad_norte);
tramos=7;

inicio=1;
final=fs*100; %segundos
semana=3.5*24*60*60*fs;
% Fourier

% este
[cxy_m f]=mscohere(velocidad_este(:,i),velocidad_este(:,s),[],[],[],fs);
[cxy_i f]=mscohere(velocidad_este(:,i),velocidad_este(:,m),[],[],[],fs);
[cxy_s f]=mscohere(velocidad_este(:,m),velocidad_este(:,s),[],[],[],fs);

cxy_smooth_i=smoothdata(cxy_i,'movmean',filtro);
cxy_smooth_s=smoothdata(cxy_s,'movmean',filtro);
cxy_smooth_m=smoothdata(cxy_m,'movmean',filtro);

plot(f,cxy_smooth_m)
hold on
plot(f,cxy_smooth_i)
plot(f,cxy_smooth_s)
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend(lgd{:})
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Este_',dia_1,'_',dia_2);

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


close(figura)

% norte

[cxy_m f]=mscohere(velocidad_norte(:,i),velocidad_norte(:,s),[],[],[],fs);
[cxy_i f]=mscohere(velocidad_norte(:,i),velocidad_norte(:,m),[],[],[],fs);
[cxy_s f]=mscohere(velocidad_norte(:,m),velocidad_norte(:,s),[],[],[],fs);

cxy_smooth_i=smoothdata(cxy_i,'movmean',filtro);
cxy_smooth_s=smoothdata(cxy_s,'movmean',filtro);
cxy_smooth_m=smoothdata(cxy_m,'movmean',filtro);

plot(f,cxy_smooth_m)
hold on
plot(f,cxy_smooth_i)
plot(f,cxy_smooth_s)
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend(lgd{:})
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Norte_',dia_1,'_',dia_2);

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


close(figura)

% arriba

[cxy_m f]=mscohere(velocidad_up(:,i),velocidad_up(:,s),[],[],[],fs);
[cxy_i f]=mscohere(velocidad_up(:,i),velocidad_este(:,m),[],[],[],fs);
[cxy_s f]=mscohere(velocidad_up(:,m),velocidad_up(:,s),[],[],[],fs);

cxy_smooth_i=smoothdata(cxy_i,'movmean',filtro);
cxy_smooth_s=smoothdata(cxy_s,'movmean',filtro);
cxy_smooth_m=smoothdata(cxy_m,'movmean',filtro);

plot(f,cxy_smooth_m)
hold on
plot(f,cxy_smooth_i)
plot(f,cxy_smooth_s)
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend(lgd{:})
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Up_',dia_1,'_',dia_2);

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


close(figura)


%% FOURIER PARA TODOS LOS BIN

clear all, close all, clc

load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_3_10;
velocidad_up=Vel_Up_3_10;
velocidad_este=Vel_East_3_10;

clear Vel_North_3_10 Vel_Up_3_10 Vel_East_3_10

dia_1=3;
dia_2=10;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

figura=1;

fs=8;%es la frecuencia de muestreo del equipo

filtro=8*60*30; %frecuencia*segundos*minutos

longitud=length(velocidad_norte);
tramos=7;

inicio=1;
final=fs*100; %segundos
semana=3.5*24*60*60*fs;
% Fourier
figure(1)
% este
for i=1:9
[cxy f]=mscohere(velocidad_este(:,i),velocidad_este(:,10),[],[],[],fs);
cxy_smooth=smoothdata(cxy,'movmean',filtro);
semilogx(f,cxy_smooth)
hold on
end
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend()
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Este_Completo_',dia_1,'_',dia_2);

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

figure(2)
% norte
for i=1:9
[cxy f]=mscohere(velocidad_norte(:,i),velocidad_norte(:,10),[],[],[],fs);
cxy_smooth=smoothdata(cxy,'movmean',filtro);
semilogx(f,cxy_smooth)
hold on
end
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend()
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Norte_Completo_',dia_1,'_',dia_2);

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

figure(3)

% arriba
for i=1:9
[cxy f]=mscohere(velocidad_up(:,i),velocidad_up(:,10),[],[],[],fs);
cxy_smooth=smoothdata(cxy,'movmean',filtro);
semilogx(f,cxy_smooth)
hold on
end
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend()
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Up_Completo_',dia_1,'_',dia_2);

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')
%%
%close(figura)

% norte

[cxy_m f]=mscohere(velocidad_norte(:,i),velocidad_norte(:,s),[],[],[],fs);
[cxy_i f]=mscohere(velocidad_norte(:,i),velocidad_norte(:,m),[],[],[],fs);
[cxy_s f]=mscohere(velocidad_norte(:,m),velocidad_norte(:,s),[],[],[],fs);

cxy_smooth_i=smoothdata(cxy_i,'movmean',filtro);
cxy_smooth_s=smoothdata(cxy_s,'movmean',filtro);
cxy_smooth_m=smoothdata(cxy_m,'movmean',filtro);

plot(f,cxy_smooth_m)
hold on
plot(f,cxy_smooth_i)
plot(f,cxy_smooth_s)
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend(lgd{:})
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Norte_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

%saveas(figura,nombre_figura,'jpeg')


%close(figura)

% arriba

[cxy_m f]=mscohere(velocidad_up(:,i),velocidad_up(:,s),[],[],[],fs);
[cxy_i f]=mscohere(velocidad_up(:,i),velocidad_este(:,m),[],[],[],fs);
[cxy_s f]=mscohere(velocidad_up(:,m),velocidad_up(:,s),[],[],[],fs);

cxy_smooth_i=smoothdata(cxy_i,'movmean',filtro);
cxy_smooth_s=smoothdata(cxy_s,'movmean',filtro);
cxy_smooth_m=smoothdata(cxy_m,'movmean',filtro);

plot(f,cxy_smooth_m)
hold on
plot(f,cxy_smooth_i)
plot(f,cxy_smooth_s)
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend(lgd{:})
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Fourier_Up_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

%saveas(figura,nombre_figura,'jpeg')


%close(figura)

%% Wavelet Corto

clear all, close all, clc

load('TIME.mat')
load('VEL_EAST.mat')
load('VEL_NORTH.mat')
load('VEL_UP.mat')

velocidad_norte=Vel_North_3_10;
velocidad_up=Vel_Up_3_10;
velocidad_este=Vel_East_3_10;

clear Vel_North_3_10 Vel_Up_3_10 Vel_East_3_10

dia_1=3;
dia_2=10;

dia_1=int2str(dia_1);
dia_2=int2str(dia_2);

fs=8;%es la frecuencia de muestreo del equipo
inicio=1;
final=fs*1000; %segundos
s=10;
i=1;
m=5;


% este
figure(1)

t=tiledlayout(3,1)

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_este(inicio:final,i),velocidad_este(inicio:final,s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 1-10')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_este(inicio:final,i),velocidad_este(inicio:final,m),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")


title('Wavelet Coherence Layer 1-5')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_este(inicio:final,m),velocidad_este(inicio:final,s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 5-10')

figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Wavelet_Este_',int2str(inicio),'_',int2str(final),'_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


%close(figura)



% norte

figure(2)

t=tiledlayout(3,1)

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_norte(inicio:final,i),velocidad_norte(inicio:final,s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 1-10')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_norte(inicio:final,i),velocidad_norte(inicio:final,m),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")


title('Wavelet Coherence Layer 1-5')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_norte(inicio:final,m),velocidad_norte(inicio:final,s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 5-10')

figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Wavelet_Norte_',int2str(inicio),'_',int2str(final),'_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


%close(figura)


% arriba

figure(3)

t=tiledlayout(3,1)

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_up(inicio:final,i),velocidad_up(inicio:final,s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 1-10')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_up(inicio:final,i),velocidad_up(inicio:final,m),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")


title('Wavelet Coherence Layer 1-5')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_up(inicio:final,m),velocidad_up(inicio:final,s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 5-10')

figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Wavelet_Up_',int2str(inicio),'_',int2str(final),'_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


%close(figura)


%% Wavelet tramos

for j=0:tramos-1   

% este

t=tiledlayout(3,1)

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_este((1+j*longitud/tramos):((j+1)*longitud/tramos),i),velocidad_este((1+j*longitud/tramos):((j+1)*longitud/tramos),s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 1-10')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_este((1+j*longitud/tramos):((j+1)*longitud/tramos),i),velocidad_este((1+j*longitud/tramos):((j+1)*longitud/tramos),m),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")


title('Wavelet Coherence Layer 1-5')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_este((1+j*longitud/tramos):((j+1)*longitud/tramos),m),velocidad_este((1+j*longitud/tramos):((j+1)*longitud/tramos),s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 5-10')

figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Wavelet_Este_',int2str(j),'_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


close(figura)

end

% norte
for j=0:tramos-1  
t=tiledlayout(3,1)

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_norte((1+j*longitud/tramos):((j+1)*longitud/tramos),i),velocidad_norte((1+j*longitud/tramos):((j+1)*longitud/tramos),s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 1-10')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_norte((1+j*longitud/tramos):((j+1)*longitud/tramos),i),velocidad_norte((1+j*longitud/tramos):((j+1)*longitud/tramos),m),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")


title('Wavelet Coherence Layer 1-5')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_norte((1+j*longitud/tramos):((j+1)*longitud/tramos),m),velocidad_norte((1+j*longitud/tramos):((j+1)*longitud/tramos),s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 5-10')

figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Wavelet_Norte_',int2str(j),'_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


close(figura)
end

% arriba
for j=0:tramos-1  

t=tiledlayout(3,1)

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_up((1+j*longitud/tramos):((j+1)*longitud/tramos),i),velocidad_up((1+j*longitud/tramos):((j+1)*longitud/tramos),s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 1-10')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_up((1+j*longitud/tramos):((j+1)*longitud/tramos),i),velocidad_up((1+j*longitud/tramos):((j+1)*longitud/tramos),m),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")


title('Wavelet Coherence Layer 1-5')

nexttile

[wcoh,wcs,f]=wcoherence(velocidad_up((1+j*longitud/tramos):((j+1)*longitud/tramos),m),velocidad_up((1+j*longitud/tramos):((j+1)*longitud/tramos),s),fs);
t=0:1/fs:(length(wcoh)-1)/fs; 
h=surface(t,f,wcoh);
h.EdgeColor = "none";
axis tight;
shading flat;
clim([0 1])
c=colorbar('Ticks',[0,0.25,0.5,0.75,1]);
c.Label.String='Coherence';
%c.Label.String='Magnitude-Squared Coherence';
%title("")
xlabel("Time (s)")
ylabel("Frequency (Hz)")

title('Wavelet Coherence Layer 5-10')

figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/',dia_1,'_',dia_2,'_ENERO/Coherencia_Wavelet_Up_',int2str(j),'_',dia_1,'_',dia_2);

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')


close(figura)

end


%% FOURIER PARA TODA LA SERIE

clear all, close all, clc

load Este.mat

fs=8;%es la frecuencia de muestreo del equipo

filtro=8*60*30; %frecuencia*segundos*minutos

longitud=length(Vel_East);
tramos=7;

inicio=1;
final=fs*100; %segundos
semana=3.5*24*60*60*fs;
% Fourier


figure(1)
% este
for i=1:9
[cxy f]=mscohere(Vel_East(:,i),Vel_East(:,10),[],[],[],fs);
cxy_smooth=smoothdata(cxy,'movmean',filtro);
semilogx(f,cxy_smooth)
hold on
end
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend()
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/Coherencia_Fourier_Este_Completo');

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

clear Vel_East

figure(2)

load Norte.mat

% norte
for i=1:9
[cxy f]=mscohere(Vel_North(:,i),Vel_North(:,10),[],[],[],fs);
cxy_smooth=smoothdata(cxy,'movmean',filtro);
semilogx(f,cxy_smooth)
hold on
end
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend()
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/Coherencia_Fourier_Norte_Completo');

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

clear Vel_North


% arriba

load Arriba.mat
figure(3)
for i=1:9
[cxy f]=mscohere(Vel_Up(:,i),Vel_Up(:,10),[],[],[],fs);
cxy_smooth=smoothdata(cxy,'movmean',filtro);
semilogx(f,cxy_smooth)
hold on
end
hold off
xlabel("Frequency (Hz)")
ylabel("Coherence")
legend()
figura=gcf;

nombre_figura=strcat('/Volumes/ECOSISTEMAS/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/FILTRADOS/TRANSFORMADOS/Coherencia_Fourier_Up_Completo');

saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

clear Vel_Up