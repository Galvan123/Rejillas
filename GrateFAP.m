clc
clear all
%-------------------------------------------------------------------------
g = 9.80665; % Aceleración debida a la gravedad (m/s^2)
% Barras en trasnversales (si, bt=1, No, bt=0);
bt=1;
L = 1; % Ancho de la rejilla en metros (eje y)
ey = 0.060; % Espaciamiento entre barras en la dirección x (longitud) en metros
ty = 0.03; % Espesor de las barras longitudinales en metros


% Barras en longitudinales (si, bl=1, No, bl=0);
bl=1;
ex = 0.48; % Espaciamiento entre barras en la dirección y (ancho) en metros
tx = 0.12; % Espesor de las barras transversales en metros
%Coefciiente de obstruccion
C=0.25;
% Barras diagonales
nd=0;
c1=[0.5020    0.5020    0.5020];
%Cv  es el coeficiente de descarga del vertedero, que generalmente tiene un valor aproximado de 1.66m1/2/s
Cv=1.66;
%Co es el coeficiente de descarga del orificio, típicamente entre 0.6 y 0.65.
Co=0.625;
%------------------------------------------------------------------------
file='GratesFAP03';
% Tirantes
Depths= readtable(file,'Sheet','Tirante','ReadVariableNames',true);
columnNames = Depths.Properties.VariableNames;
tiempo=Depths(:,1);
tiempo=table2array(tiempo(:, columnNames{1}))
h=hour(tiempo);
m=minute(tiempo);
hora=h+m/60;
time=datenum(tiempo);
dt=0.75*(time(end)-time(1));
dt=time(1)+dt;
Depths=Depths(:,2:end);
Depths=table2array(Depths);
% Caudales
Flows= readtable(file,'Sheet','Caudal','ReadVariableNames',true);
Flows=Flows(:,2:end);
Flows=table2array(Flows);
% Velocidades
Speeds= readtable(file,'Sheet','Velocidad','ReadVariableNames',true);
Speeds=Speeds(:,2:end);
Speeds=table2array(Speeds);
%calles
streets= readtable(file,'Sheet','Calles','ReadVariableNames',true);
namestreet=streets(:,1:2);
namestreet=table2array(namestreet);
width=streets(:,3:end);
width=table2array(width);
%--------------------------------------------------------------------------
c1=[0.4667    0.6745    0.1882];
c2=[0.8510    0.3255    0.0980];
c3=[0.4941    0.1843    0.5569];
c4=[0.4667    0.6745    0.1882];
ns=length(width(:,1));
nq=length(Flows(:,1));
for i=1:ns
    sector=namestreet{i,1};
    name=namestreet{i,2};
    W=width(i,1);% Ancho de la calle
    W=W-0.6;% Ancho de la rejilla
    W=round(10*W)/10
    % Cálculo de la cantidad de barras
    nx = round((W + ex) / (ex + tx))
    %dx = L - nx * (ex + tx) + ex;
    ex = (W-tx*nx)/(nx - 1)
    
    ny = round((L + ey) / (ey + ty))
    %dy = W - ny * (ey + ty) + ey;
    ey = (L-ty*ny)/(ny - 1)
    
    % Area vacia o efectiva
    AH=ex*ey*(nx-1)*(ny-1);
    Ah=AH*10^4
    % Longitud efectiva
    LE=ey*(ny-1);
    Le=100*LE
    % Ancho efectivo
    Be=ex*(nx-1)
    % Parámetros geométricos característicos de cada rejilla
    alfa = -1.924 * (Le^0.631/Ah^0.279) * (1 + nd)^-0.089 * (nx + 1)^-0.238 * (ny + 1)^-0.045
    beta = -26.803* Le^-4.953+1.213
    dimensions(i,1)=L;
    dimensions(i,2)=W;
    dimensions(i,3)=nx;
    dimensions(i,4)=ex;
    dimensions(i,5)=tx;
    dimensions(i,6)=ny;
    dimensions(i,7)=ey;
    dimensions(i,8)=ty;
    % Mostrar los resultados
    %     fprintf('Cantidad de barras en la dirección longitudinal (x): %d\n', nx);
    %     fprintf('Cantidad de barras en la dirección transversal (y): %d\n', ny);
    Q=Flows(:,i);
    Y=Depths(:,i);
    V=Speeds(:,i);
    % Eficiencia variqble en el tiempo
    %C = randi([24 26])/100  
    C=0.25
    for j=1:nq
        q=Q(j);
        if q<0
            q=0;
            Q(j)=q;
        end
        y=Y(j);
        v=V(j);
        Qm=AH*max(V);
        F=v/sqrt(g*y);
        FF(j,i)=F;        
        E=alfa*F*(100*y/Le)^0.812+beta;
        if E>1
            E=1;
        end  
        BB(j,i)=E;
        E=(1-C)*E;
        EE(j,i)=E;
        FY(j,i)=(E-beta)/alfa;
        qu=q/Be;
        QU(j,i)=qu;
        QI(j,i)=E*qu*Be;
%       qq(j,i)=E*qu*Be;
%       qq(j,i+1)=E*q;
        DQ(j,i)=q-E*qu*Be;
    end
    Qmax=max(Q);
    Qin=max(QI(:,i))
    Qv=Cv*Be*max(Y)^1.5;
    Qo=Co*AH*sqrt(2*g*y);
    Emin=min(EE(:,i));    
    SW=[hora,QI(:,i)];
    SW = SW(~isnan(SW(:,2)), :);
    nam=sprintf(strcat(name,'.dat'));  
    save(nam, 'SW', '-ascii');
    %hold on
    %plot(time,QYV(:,1),'--')
    %plot(time,QYV(:,5),'-')
    %plot(QYV(:,4),QYV(:,5),'-')
    % Crear una figura para los subgráficos
    %figure;
    % Crear una nueva figura
    fig = figure;
    % Establecer el tamaño y la posición de la figura
    fig.Position = [200, 200, 1000, 600]; % [left, bottom, width, height]
    % Primer subgráfico:
    subplot(2, 2, 1); % 2 filas, 2 columnas, posición 1
    plot(time, Y,'-','color',c1,'linewidth',2);
    title('Tirante');
    xlabel('Tiempo (h)','FontSize',10,'FontWeight','bold','Color','k');
    ylabel('y (m)','FontSize',10,'FontWeight','bold','Color','k');
    grid on
    box on
    datetick('x', 'HH:MM', 'keepticks');
    % Segundo subgráfico:
    subplot(2, 2, 2); % 2 filas, 2 columnas, posición 2
    plot(time, V,'-','color',c2,'linewidth',2);
    title('Velocidad');
    xlabel('Tiempo (h)','FontSize',10,'FontWeight','bold','Color','k');
    ylabel('V (m/s)','FontSize',10,'FontWeight','bold','Color','k');
    grid on
    box on
    datetick('x', 'HH:MM', 'keepticks');
    
    % Tercer subgráfico:
    subplot(2, 2, 3); % 2 filas, 2 columnas, posición 3
    plot(FY(:,i), EE(:,i),'-','color',c3,'linewidth',2)
    title('Eficiencia & Froude');
    xlabel('F*(y/L)^{0.812}','FontSize',10,'FontWeight','bold','Color','k');
    ylabel('Eficiencia','FontSize',10,'FontWeight','bold','Color','k');
    grid on
    box on
    %datetick('x', 'HH:MM', 'keepticks');
    
    % Cuarto subgráfico:    
    subplot(2, 2, 4); % 2 filas, 2 columnas, posición 4
    hold on
    %plot(time, QYV(:,1), 'DisplayName', 'Caudal superficial'); % Primera serie de datos
    plot(time, Q,'-','color','b','linewidth',2,'DisplayName', 'Caudal Superficial');
    plot(time, QU(:,i),'-.','color','m','linewidth',2,'DisplayName', 'Caudal unitario');
    plot(time, QI(:,i),'-.','color','g','linewidth',2,'DisplayName', 'Caudal Interceptado');
    plot(time, DQ(:,i),'-','color','r','linewidth',2,'DisplayName', 'Caudal No Interceptado');
    
    text(dt,Qmax/3,strcat('E=',num2str(Emin)),'Color','r','FontWeight','bold','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','bottom');
    text(dt,Qmax/6,strcat('Qin=',num2str(Qin)),'Color','r','FontWeight','bold','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','bottom');
    title('Caudales');
    xlabel('Tiempo (h)','FontSize',9,'FontWeight','bold','Color','k');
    ylabel('Q (m3/s)','FontSize',9,'FontWeight','bold','Color','k');
    datetick('x', 'HH:MM', 'keepticks');
    grid on
    box on
    legend('show', 'Location', 'best'); % Muestra la leyenda con los nombres especificados
    name4=strcat('Eficiencia-',name);
    sgtitle(name);
    pdffigurec(name4);
    % Crear la rejilla
    figure;
    H=4*L; % Largo deseado, tres veces el ancho
    set(gcf, 'Position', [100, 100, 250*W, 200*H]);
    
    % Dibujar las barras en la dirección x (longitudinales)
    %----------------------------------------------------------------------
    % Generar archivo DXF
    %name2 = street{i,2};
    name3 = strrep(name, '_', ' ');
    name3 = regexprep(name3, '\s+', ' ');
    name3 = strtrim(name3);
    dxf_filename = strcat(name,'.dxf');
    fileID = fopen(dxf_filename, 'w');
    % Escribir el encabezado del archivo DXF
    fprintf(fileID, '0\nSECTION\n2\nHEADER\n0\nENDSEC\n');
    fprintf(fileID, '0\nSECTION\n2\nTABLES\n0\nENDSEC\n');
    fprintf(fileID, '0\nSECTION\n2\nBLOCKS\n0\nENDSEC\n');
    fprintf(fileID, '0\nSECTION\n2\nENTITIES\n');
    for i = 0:nx-1
        x = i * (ex + tx);
        rectangle('Position', [x, 0, tx, L], 'FaceColor', 'b');
        % Escribir un rectángulo en la dirección x
        fprintf(fileID, '0\nPOLYLINE\n8\n0\n66\n1\n70\n1\n');
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n0\n30\n0\n', x);
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', x, L);
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', x+tx, L);
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n0\n30\n0\n', x+tx);
        fprintf(fileID, '0\nSEQEND\n');
    end
    % Dibujar las barras en la dirección y (transversales)
    for j = 0:ny-1
        y = j * (ey + ty);
        rectangle('Position', [0, y, W, ty], 'FaceColor', 'r');
        % Escribir un rectángulo en la dirección y
        fprintf(fileID, '0\nPOLYLINE\n8\n0\n66\n1\n70\n1\n');
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n0\n20\n%.3f\n30\n0\n', y);
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', W, y);
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', W, y+ty);
        fprintf(fileID, '0\nVERTEX\n8\n0\n10\n0\n20\n%.3f\n30\n0\n', y+ty);
        fprintf(fileID, '0\nSEQEND\n');
    end
    rectangle('Position', [0, 0, W, L],'EdgeColor',c1,'LineWidth',2);
    % Dibujar el rectángulo en AutoCAD como una polilínea cerrada
    fprintf(fileID, '0\nPOLYLINE\n8\n0\n66\n1\n70\n1\n');
    fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', 0, 0);
    fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', W, 0);
    fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', W, L);
    fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', 0, L);
    fprintf(fileID, '0\nVERTEX\n8\n0\n10\n%.3f\n20\n%.3f\n30\n0\n', 0, 0);
    fprintf(fileID, '0\nSEQEND\n');
    
    % Cerrar la sección de ENTITIES y el archivo DXF
    fprintf(fileID, '0\nENDSEC\n0\nEOF\n');
    fclose(fileID);
    %--------------------------------------------------------------------------
    d=sqrt(L^2+W^2);
    d=0.04*d;
    % Añadir anotaciones para la longitud y el ancho de la rejilla
    %annotation('doublearrow', [0, 0.9], [0.05, 0.05], 'HeadStyle', 'vback1', 'LineStyle', '-', 'Color', 'g');
    
    text(W/2, -d/2, sprintf('W = %.1f m',W), 'HorizontalAlignment', 'center','FontSize',15,'FontWeight','bold');
    
    %annotation('doublearrow', [0.05, 0.05], [0.1, 0.9], 'HeadStyle', 'vback1', 'LineStyle', '--', 'Color', 'r');
    text(-d/2, L/2, sprintf('L = %.1f m',L), 'Rotation',90,'HorizontalAlignment','center','FontSize',15,'FontWeight','bold');
    
    e=W/6;
    % Añadir texto con información adicional
    text(0, L+d/2, sprintf('n_x=%d', nx), 'HorizontalAlignment', 'left','FontSize',13,'FontWeight','bold');
    text(e + 0.2, L+d/2, sprintf('e_x=%.2f m', ex), 'HorizontalAlignment', 'left','FontSize',13,'FontWeight','bold');
    text(2*e + 0.2, L+d/2, sprintf('t_x=%.2f m', tx), 'HorizontalAlignment', 'left','FontSize',13,'FontWeight','bold');
    
    text(3*e + 0.2, L+d/2, sprintf('n_y=%d', ny), 'HorizontalAlignment', 'left','FontSize',13,'FontWeight','bold');
    text(4*e + 0.2, L+d/2, sprintf('e_y=%.2f m', ey), 'HorizontalAlignment', 'left','FontSize',13,'FontWeight','bold');
    text(5*e + 0.2, L+d/2, sprintf('t_y=%.2f m', ty), 'HorizontalAlignment', 'left','FontSize',13,'FontWeight','bold');
    
    % Configurar la gráfica
    xlabel('Longitud (m)');
    ylabel('Ancho (m)');
    title(strcat('Rejilla', ' - ',name3));
    %axis equal;
    xlim([-d, W+d]);
    ylim([-d, L+d]);
    grid on;
    box on;
    axis on;
    hold on;
    % Abrir el archivo DXF en AutoCAD (opcional)
    %system(['start ', dxf_filename]);
    pdffigurec(name); 
end
% Convertir la matriz de valores a una tabla
Rejillas = array2table(dimensions, 'VariableNames', {'L', 'W', 'nx', 'ex','tx','ny', 'ey','ty'});
% Concatenar las tablas
VRejillas = [streets(:,1:3), Rejillas];
%Exportar la tabla a un archivo Excel
writetable(VRejillas, 'Dimensiones_Rejillas.xlsx');
% Exportar numero de froude, eficiencia y Qin
% Nombre del archivo Excel de salida
filename = 'Eficiencia_Caudal.xlsx';
% Escribir FF en la primera hoja
FR= array2table([time FF], 'VariableNames', columnNames)
writetable(FR, filename, 'Sheet', 'Froude');
% Escribir EE en la segunda hoja
EF= array2table([time EE], 'VariableNames', columnNames)
writetable(EF, filename, 'Sheet', 'Eficiencia');
% Escribir QI en la tercera hoja
QQ= array2table([time QI], 'VariableNames', columnNames)
writetable(QQ, filename, 'Sheet', 'Caudal_Interceptado');
%Caudal no interceptado
QN= array2table([time DQ], 'VariableNames', columnNames)
writetable(QN, filename, 'Sheet', 'Caudal_No_Interceptado');
%--------------------------------------------------------------------------
function pdffigurec(namef)
nf=strcat(namef,'.png');
%set(gcf,'Units','inches')
%axis(gca,'equal');
set(gca,'FontSize',10,'FontName','Arial', 'FontUnits','normalized');
exportgraphics(gcf,nf,'Resolution',1500);
end

