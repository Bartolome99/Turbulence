clear

load("Este.mat")

data =abs(Vel_East) ; % Matriz de datos 
clear Vel_East
% Crear un boxplot horizontal
figure(1)
figura = figure(1);
%boxplot(data, 'Orientation', 'horizontal', 'MarkerEdgeColor', [0.5 0.5 0.5]); % Puntos grises
boxplot(data, 'Orientation', 'horizontal', 'Symbol', 'kx'); % 'ko' para puntos negros); punto x r rojo g verde

% Etiquetas opcionales para los ejes
% Aumentar el tamaño de las letras en los ejes y leyendas
size_letra=16;
xlabel('Velocity eastward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Bins', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
xlim([0 4])

% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/','Box_plot_East');

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

close()

clear data

load("Norte.mat")

data =abs(Vel_North) ; % Matriz de datos 
clear Vel_North
% Crear un boxplot horizontal
figure(2)
figura = figure(2);
boxplot(data, 'Orientation', 'horizontal', 'Symbol', 'kx'); % 'ko' para puntos negros); punto x r rojo g verde

% Etiquetas opcionales para los ejes
% Aumentar el tamaño de las letras en los ejes y leyendas
size_letra=16;
xlabel('Velocity northward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Bins', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
xlim([0 4])

% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/','Box_plot_North');

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

close()

clear data

load("Arriba.mat")

data =abs(Vel_Up) ; % Matriz de datos 
clear Vel_Up
% Crear un boxplot horizontal
figure(3)
figura = figure(3);
boxplot(data, 'Orientation', 'horizontal', 'Symbol', 'kx'); % 'ko' para puntos negros); punto x r rojo g verde

% Etiquetas opcionales para los ejes
% Aumentar el tamaño de las letras en los ejes y leyendas
size_letra=16;
xlabel('Velocity upward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Bins', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
xlim([-4 4])

% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/','Box_plot_Up');

%saveas(figura,nombre_figura,'fig')

saveas(figura,nombre_figura,'jpeg')

saveas(figura,nombre_figura,'eps')

saveas(figura,nombre_figura,'pdf')

close()

clear data

%% Box-Plot semana

clear

load("VEL_EAST.mat")

data =Vel_East_3_10 ; % Matriz de datos 
clear Vel_East_3_10
% Crear un boxplot horizontal
figure(1)
figura = figure(1);
%boxplot(data, 'Orientation', 'horizontal', 'MarkerEdgeColor', [0.5 0.5 0.5]); % Puntos grises
boxplot(data, 'Orientation', 'horizontal', 'Symbol', 'kx'); % 'ko' para puntos negros); punto x r rojo g verde

% Etiquetas opcionales para los ejes
% Aumentar el tamaño de las letras en los ejes y leyendas
size_letra=16;
xlabel('Velocity eastward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Bins', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
xlim([-4 4])

% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

% nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/','Box_plot_East');
% 
% %saveas(figura,nombre_figura,'fig')
% 
% saveas(figura,nombre_figura,'jpeg')
% 
% saveas(figura,nombre_figura,'eps')
% 
% saveas(figura,nombre_figura,'pdf')
% 
% close()

clear data

load("VEL_NORTH.mat")

data =Vel_North_3_10 ; % Matriz de datos 
clear Vel_North_3_10
% Crear un boxplot horizontal
figure(2)
figura = figure(2);
boxplot(data, 'Orientation', 'horizontal', 'Symbol', 'kx'); % 'ko' para puntos negros); punto x r rojo g verde

% Etiquetas opcionales para los ejes
% Aumentar el tamaño de las letras en los ejes y leyendas
size_letra=16;
xlabel('Velocity northward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Bins', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
xlim([-4 4])

% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

% nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/','Box_plot_North');
% 
% %saveas(figura,nombre_figura,'fig')
% 
% saveas(figura,nombre_figura,'jpeg')
% 
% saveas(figura,nombre_figura,'eps')
% 
% saveas(figura,nombre_figura,'pdf')
% 
% close()

clear data

load("VEL_UP.mat")

data =Vel_Up_3_10 ; % Matriz de datos 
clear Vel_Up_3_10
% Crear un boxplot horizontal
figure(3)
figura = figure(3);
boxplot(data, 'Orientation', 'horizontal', 'Symbol', 'kx'); % 'ko' para puntos negros); punto x r rojo g verde

% Etiquetas opcionales para los ejes
% Aumentar el tamaño de las letras en los ejes y leyendas
size_letra=16;
xlabel('Velocity upward component (m/s)', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje X
ylabel('Bins', 'FontSize', size_letra, 'FontWeight', 'bold'); % Etiqueta del eje Y
title('', 'FontSize', size_letra, 'FontWeight', 'bold'); % Título
xlim([-4 4])

% Aumentar el tamaño de las etiquetas del eje Y (nombres de las columnas)
set(gca, 'FontSize', size_letra); % Cambia el tamaño de fuente de los números en los ejes

% nombre_figura=strcat('/Users/bartolome/Documents/UPCT/ARTÍCULOS/TURBULENCIA/','Box_plot_Up');
% 
% %saveas(figura,nombre_figura,'fig')
% 
% saveas(figura,nombre_figura,'jpeg')
% 
% saveas(figura,nombre_figura,'eps')
% 
% saveas(figura,nombre_figura,'pdf')
% 
% close()

clear data
