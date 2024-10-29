%% Histogramas
clear all

load Arriba.mat

figura=figure(1)

[xx nn] = hist(Vel_Up(:,:),100);
j=1;
for i=[1,5,10]
    cum_xx(:,j) = xx(:,i)/sum(xx(:,i));
    j=j+1;
end
plot(nn,cum_xx,'LineWidth',2)

size_letra=16;
xlim([-4 4])
ylim([0 0.5])
legend({'Bin 1','Bin 5','Bin 10'},'FontSize', size_letra, 'FontWeight', 'bold')
xlabel('Velocity upward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Probability density function', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

%title('Histogram Velocity Up')

nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/GRAFICOS/','Histogram_Velocity_Up');

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')


clear all

load Este.mat


figura=figure(2)

[xx nn] = hist(Vel_East(:,:),100);
j=1;
for i=[1,5,10]
    cum_xx(:,j) = xx(:,i)/sum(xx(:,i));
    j=j+1;
end
plot(nn,cum_xx,'LineWidth',2)

size_letra=16;
xlim([-4 4])
ylim([0 0.5])
legend({'Bin 1','Bin 5','Bin 10'},'FontSize', size_letra, 'FontWeight', 'bold')
xlabel('Velocity eastward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Probability density function', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

%title('Histogram Velocity East')

nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/GRAFICOS/','Histogram_Velocity_East');

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')


clear all

load Norte.mat


figura=figure(3)

[xx nn] = hist(Vel_North(:,:),100);
j=1;
for i=[1,5,10]
    cum_xx(:,j) = xx(:,i)/sum(xx(:,i));
    j=j+1;
end
plot(nn,cum_xx,'LineWidth',2)

size_letra=16;
xlim([-4 4])
ylim([0 0.5])
legend({'Bin 1','Bin 5','Bin 10'},'FontSize', size_letra, 'FontWeight', 'bold')
xlabel('Velocity northward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Probability density function', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

%title('Histogram Velocity North')

nombre_figura=strcat('/Volumes/NO NAME/TURBULENCIA/ADCP/103062/ARCHIVOS_MATLAB/GRAFICOS/','Histogram_Velocity_North');

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

