%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Del4_OptimalPlacement_Enhanced.m
%
% MATLAB-skript för Del 4 – Optimal placering av ljudkälla för att minimera ljud i sovhörnan
%
% Projektet handlar om att hitta optimal placering av en TV (ljudkälla)
% för att minimera ljudnivån i en sovhörna i ett L-format rum.
%
% Projekt: Ljudvågor i ett slutet rum (SF1550 VT 2025)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;  % Rensa alla variabler och figurer

%% ------------------------- Parametrar ----------------------------------
% Grundläggande parametrar för projektet
omega = 30;         % Frekvens för ljudvågen
aTV = 1;            % Källstyrka (satt till 1 enligt uppgiften)
yTV_fixed = 0.6;    % Y-koordinat för TV:n när den sitter på väggen vid x=0.6

% Upplösningsparametrar för beräkningar
N_coarse = 200;     % Låg upplösning för snabbare beräkningar i grova steg
N_fine = 1000;      % Hög upplösning för noggrannare resultat i finare steg

% Ljudkällans form - en funktion som avtar med avståndet
S0_handle = @(r) cos(24*r) .* exp(-900*r.^2);

% Definition av sovhörnan: Ω_sov = [0, 0.25] x [0.5, 0.75]
% Sovhörnan är området som vi vill minimera ljudet i
sov_x_max = 0.25;   % Högsta x-värde för sovhörnan
sov_y_min = 0.5;    % Lägsta y-värde för sovhörnan
sov_y_max = 0.75;   % Högsta y-värde för sovhörnan

%% ------------------------- Del 4.1: Grovlokalisering ---------------------
% I denna del testar vi olika TV-positioner på väggen för att få en grov
% uppfattning om var ljudet i sovhörnan är som lägst
fprintf('=== Del 4.1: Grovlokalisering av optimal TV-placering på väggen ===\n');

% Testa 20 jämnt fördelade positioner på väggen mellan x=0.6 och x=1
xs_values = linspace(0.6, 1, 20);   % 20 x-positioner på väggen
A_values = zeros(size(xs_values));   % Förbered array för relativa ljudstyrkor

% Beräkna relativa ljudstyrkan A för varje position
for i = 1:length(xs_values)
    xs = xs_values(i);
    % Anropa funktionen som beräknar relativa ljudstyrkan i sovhörnan
    A_values(i) = evaluateA_x(xs, yTV_fixed, omega, aTV, S0_handle, N_coarse, sov_x_max, sov_y_min, sov_y_max);
    fprintf('xs = %.3f, A = %.4f\n', xs, A_values(i));
end

% Hitta bästa och sämsta position för visualisering senare
[~, best_idx] = min(A_values);      % Index för lägsta ljudnivå (bäst)
[~, worst_idx] = max(A_values);     % Index för högsta ljudnivå (sämst)
best_xs = xs_values(best_idx);      % x-koordinat för bästa position
worst_xs = xs_values(worst_idx);    % x-koordinat för sämsta position

% Skapa en figur som visar resultaten
figure;
% Rita grundläggande kurvan för A(xs)
plot(xs_values, A_values, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
% Markera bästa (lägsta) och sämsta (högsta) värden
plot(best_xs, A_values(best_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(worst_xs, A_values(worst_idx), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
% Lägg till axelanvändning och titel
xlabel('TV placering xs (på väggen)', 'FontSize', 12);
ylabel('Relativ ljudstyrka A', 'FontSize', 12);
title('Grovlokalisering: A(xs) med N = 200', 'FontSize', 14);
grid on;
legend('A(xs)', 'Min A', 'Max A', 'Location', 'best');

%% -------------------- Del 4.2: Gyllene snitt-sökning ---------------------
% I denna del använder vi gyllene snitt-metoden för att noggrannare hitta 
% den optimala positionen på väggen
fprintf('\n=== Del 4.2: Gyllene snitt-sökning för optimalt xs med N = %d ===\n', N_fine);

% Sökintervall för xs på väggen
a_int = 0.6;    % Vänstra väggen
b_int = 1.0;    % Högra änden av rummet
tol = 1e-4;     % Tolerans för sökningen (hur nära vi vill komma det optimala)

% Använder gyllene snitt-sökning för att hitta optimalt xs
[opt_xs, opt_A] = goldenSectionSearch(@(xs) evaluateA_x(xs, yTV_fixed, omega, aTV, S0_handle, N_fine, sov_x_max, sov_y_min, sov_y_max), a_int, b_int, tol);

fprintf('Optimalt xs = %.5f med A = %.5f\n', opt_xs, opt_A);

% Skapa en mer detaljerad graf för resultatet med finare upplösning
xs_fine = linspace(a_int, b_int, 50);  % 50 punkter i intervallet
% Beräkna A för varje punkt
A_fine = arrayfun(@(xs) evaluateA_x(xs, yTV_fixed, omega, aTV, S0_handle, N_fine, sov_x_max, sov_y_min, sov_y_max), xs_fine);
figure;
% Rita kurvan med högre noggrannhet
plot(xs_fine, A_fine, 'r.-', 'LineWidth', 1.5);
hold on;
% Markera den optimala punkten
plot(opt_xs, opt_A, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
% Lägg till axelanvändning och titel
xlabel('TV placering xs (på väggen)', 'FontSize', 12);
ylabel('Relativ ljudstyrka A', 'FontSize', 12);
title('Gyllene snitt-sökning: A(xs) med N = 1000', 'FontSize', 14);
grid on;
legend('A(xs)', 'Optimalt xs', 'Location', 'best');

%% --------------- Del 4.3: Tvådimensionell optimering ----------------------
% I denna del optimerar vi TV-positionen i både x- och y-led,
% vilket betyder att den kan placeras var som helst i rummet, inte bara på väggen
fprintf('\n=== Del 4.3: Tvådimensionell optimering ===\n');

% Definiera objektfunktionen - vad vi vill minimera
% Vi använder en robust version som hanterar begränsningar för rummets form
objFun = @(x) evaluateA_xy(x, omega, aTV, S0_handle, N_coarse, sov_x_max, sov_y_min, sov_y_max);

% Välj flera startgissningar för att undvika att fastna i lokala minima
start_guesses = [0.7, 0.59;    ... % Nära väggen, x-koordinat nära optimala från Del 4.2
                 0.8, 0.58;    ... % Längre höger, nära väggen
                 0.9, 0.57;    ... % Ännu längre höger, nära väggen
                 0.65, 0.55];      % Vänster, lite längre från väggen

% Körparametrar för optimering
options = optimset('Display', 'iter', 'TolX', 1e-4, 'MaxFunEvals', 200, 'MaxIter', 100);
results = struct('x_opt',{}, 'A_opt',{});

% Steg 1: Grov optimering med lägre upplösning
fprintf('Steg 1: Grov optimering med N = %d\n', N_coarse);
for j = 1:size(start_guesses,1)
    x0 = start_guesses(j,:);
    
    % Direkt kontroll om startpositionen är giltig
    xs = x0(1);
    ys = x0(2);

% Kontrollera om positionen är giltig
    if xs < 0 || xs > 1 || ys < 0 || ys > 0.75 || (xs > 0.25 && xs < 0.6 && ys > 0.4 && ys < 0.6) || ys >= 0.6
        fprintf('Startpunkt [%.2f, %.2f] är ogiltig, hoppar över\n', xs, ys);
        continue;
    end
    
    % Utför optimering
    [xj, Aj] = fminsearch(objFun, x0, options);
    
    % Projicera resultat om utanför rummet
    xs = xj(1);
    ys = xj(2);
    
    % Kontrollera om resultatet är utanför rummet och projicera in det
    outsideRoom = xs < 0 || xs > 1 || ys < 0 || ys > 0.75 || (xs > 0.25 && xs < 0.6 && ys > 0.4 && ys < 0.6) || ys >= 0.6;
    
    if outsideRoom
        fprintf('  Resultat [%.2f, %.2f] är utanför giltigt område - justerar tillbaka\n', xs, ys);
        
        % Projicera in i rummet
        xs = max(0, min(1, xs));
        ys = max(0, min(0.75, ys));
        
        % Hantera väggen - allt ovanför går till under
        if ys >= 0.6
            ys = 0.59;
        end
        
        % Hantera L-utskärningen - simpel projection till höger
        if xs > 0.25 && xs < 0.6 && ys > 0.4 && ys < 0.6
            xs = 0.61; % projicera höger om L-utskärning
        end
        
        % Uppdatera xj och Aj
        xj = [xs, ys];
        Aj = evaluateA_xy_simplified(xj, omega, aTV, S0_handle, N_coarse, sov_x_max, sov_y_min, sov_y_max);
    end
    
    results(j).x_opt = xj;
    results(j).A_opt = Aj;
    fprintf('Start [%0.2f, %0.2f] -> xs = %.5f, ys = %.5f, A = %.5f\n', x0(1), x0(2), xj(1), xj(2), Aj);
end

% Välj bästa av resultaten
temp_A = [results.A_opt];
[~, best_coarse] = min(temp_A);
x_best_coarse = results(best_coarse).x_opt;
A_best_coarse = results(best_coarse).A_opt;

fprintf('\nBästa grova resultat: xs = %.5f, ys = %.5f, A = %.5f\n', x_best_coarse(1), x_best_coarse(2), A_best_coarse);

% Steg 2: Fin optimering med högre upplösning
fprintf('\nSteg 2: Fin optimering med N = %d\n', N_fine);
objFun_fine = @(x) evaluateA_xy(x, omega, aTV, S0_handle, N_fine, sov_x_max, sov_y_min, sov_y_max);

[x_opt, A_opt] = fminsearch(objFun_fine, x_best_coarse, options);

% Kontrollera och projicera resultat om utanför rummet
xs = x_opt(1);
ys = x_opt(2);

% Kontrollera om resultatet är utanför rummet
outsideRoom = xs < 0 || xs > 1 || ys < 0 || ys > 0.75 || (xs > 0.25 && xs < 0.6 && ys > 0.4 && ys < 0.6) || ys >= 0.6;

if outsideRoom
    fprintf('  Slutligt resultat [%.2f, %.2f] är utanför giltigt område - justerar tillbaka\n', xs, ys);
    
    % Projicera in i rummet
    xs = max(0, min(1, xs));
    ys = max(0, min(0.75, ys));
    
    % Hantera väggen - allt ovanför går till under
    if ys >= 0.6
        ys = 0.59;
    end
    
    % Hantera L-utskärningen - simpel projection till höger
    if xs > 0.25 && xs < 0.6 && ys > 0.4 && ys < 0.6
        xs = 0.61; % projicera höger om L-utskärning
    end
    
    % Uppdatera x_opt och A_opt
    x_opt = [xs, ys];
    A_opt = objFun_fine(x_opt);
end

fprintf('Optimal 2D-placering: xs = %.5f, ys = %.5f, A = %.5f\n', x_opt(1), x_opt(2), A_opt);
fprintf('Förbättring jämfört med väggplacering: %.2f%%\n', 100*(opt_A - A_opt)/opt_A);

%% ------------------ Visualisering av resultat ----------------------------
% I denna del visualiserar vi resultaten genom att beräkna och visa ljudfälten 
% för de olika placeringarna

% Beräkna ljudfälten för optimala (2D och vägg) och dåliga placeringar

% Optimal 2D placering - ljudkälla placerad på optimala koordinater
S_opt_2D = @(x,y) aTV * S0_handle(sqrt((x-x_opt(1)).^2 + (y-x_opt(2)).^2));
[~, Sol_opt_2D] = hhsolver(omega, S_opt_2D, N_fine);

% Optimal väggplacering - ljudkälla placerad på väggen
S_opt_wall = @(x,y) aTV * S0_handle(sqrt((x-opt_xs).^2 + (y-yTV_fixed).^2));
[~, Sol_opt_wall] = hhsolver(omega, S_opt_wall, N_fine);

% Dålig placering - ljudkälla placerad på sämsta stället på väggen
S_bad = @(x,y) aTV * S0_handle(sqrt((x-worst_xs).^2 + (y-yTV_fixed).^2));
[~, Sol_bad] = hhsolver(omega, S_bad, N_fine);

%% Visualisering 1: 2D Contourplot av lösningarna
% Skapa konturer som visar ljudets amplitud i rummet

% Definiera koordinater för rummet och sovhörnan för att rita dem
% Rummet har en L-form
room_x = [0, 0, 0.25, 0.25, 0.6, 0.6, 1, 1, 0];  % x-koordinater för rummets hörn
room_y = [0, 0.75, 0.75, 0.4, 0.4, 0.6, 0.6, 0, 0];  % y-koordinater för rummets hörn
% Sovhörnan är ett rektangulärt område i övre vänstra hörnet
sov_x = [0, sov_x_max, sov_x_max, 0, 0];  % x-koordinater för sovhörnans hörn
sov_y = [sov_y_min, sov_y_min, sov_y_max, sov_y_max, sov_y_min];  % y-koordinater

figure('Position', [100, 100, 1200, 400]);

% Rita optimal 2D-placering
subplot(1,3,1);
contour(Sol_opt_2D.x, Sol_opt_2D.y, abs(Sol_opt_2D.u), 20);  % 20 konturlinjer
axis equal;  % Samma skala i x och y
hold on;
% Rita rummet och sovhörnan
plot(room_x, room_y, 'k-', 'LineWidth', 2);  % Rummets gränser
patch(sov_x, sov_y, 'r', 'FaceAlpha', 0.2);  % Sovhörnan (röd transparent)
% Markera källans position
plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
title(sprintf('Optimal 2D: (%.2f, %.2f), A = %.4f', x_opt(1), x_opt(2), A_opt), 'FontSize', 12);
axis off;

% Rita optimal väggplacering
subplot(1,3,2);
contour(Sol_opt_wall.x, Sol_opt_wall.y, abs(Sol_opt_wall.u), 20);
axis equal;
hold on;
plot(room_x, room_y, 'k-', 'LineWidth', 2);
patch(sov_x, sov_y, 'r', 'FaceAlpha', 0.2);
plot(opt_xs, yTV_fixed, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
title(sprintf('Optimal på vägg: (%.2f, %.2f), A = %.4f', opt_xs, yTV_fixed, opt_A), 'FontSize', 12);
axis off;

% Rita dålig placering
subplot(1,3,3);
contour(Sol_bad.x, Sol_bad.y, abs(Sol_bad.u), 20);
axis equal;
hold on;
plot(room_x, room_y, 'k-', 'LineWidth', 2);
patch(sov_x, sov_y, 'r', 'FaceAlpha', 0.2);
plot(worst_xs, yTV_fixed, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
title(sprintf('Dålig placering: (%.2f, %.2f), A = %.4f', worst_xs, yTV_fixed, A_values(worst_idx)), 'FontSize', 12);
axis off;

%% Visualisering 2: Färgkarta med shading
% Skapa färgkartor som visar ljudvågorna med jämn skuggning

figure('Position', [100, 100, 1200, 400]);

% Optimal 2D
subplot(1,3,1);
surf(Sol_opt_2D.x, Sol_opt_2D.y, abs(Sol_opt_2D.u));  % Skapa 3D-yta för ljudamplitud
view(0, 90);  % Visa uppifrån (2D-vy)
shading interp;  % Jämn interpolerad skuggning
colorbar;  % Visa färgskala
axis equal;
hold on;
% Rita rummet och sovhörnan ovanpå 3D-ytan
plot3(room_x, room_y, ones(size(room_x))*max(abs(Sol_opt_2D.u(:))), 'k-', 'LineWidth', 2);
plot3(sov_x, sov_y, ones(size(sov_x))*max(abs(Sol_opt_2D.u(:))), 'r-', 'LineWidth', 2);
plot3(x_opt(1), x_opt(2), max(abs(Sol_opt_2D.u(:))), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
title(sprintf('Optimal 2D: (%.2f, %.2f), A = %.4f', x_opt(1), x_opt(2), A_opt), 'FontSize', 12);

% Optimal vägg
subplot(1,3,2);
surf(Sol_opt_wall.x, Sol_opt_wall.y, abs(Sol_opt_wall.u));
view(0, 90);
shading interp;
colorbar;
axis equal;
hold on;
plot3(room_x, room_y, ones(size(room_x))*max(abs(Sol_opt_wall.u(:))), 'k-', 'LineWidth', 2);
plot3(sov_x, sov_y, ones(size(sov_x))*max(abs(Sol_opt_wall.u(:))), 'r-', 'LineWidth', 2);
plot3(opt_xs, yTV_fixed, max(abs(Sol_opt_wall.u(:))), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
title(sprintf('Optimal på vägg: (%.2f, %.2f), A = %.4f', opt_xs, yTV_fixed, opt_A), 'FontSize', 12);

% Dålig placering
subplot(1,3,3);
surf(Sol_bad.x, Sol_bad.y, abs(Sol_bad.u));
view(0, 90);
shading interp;
colorbar;
axis equal;
hold on;
plot3(room_x, room_y, ones(size(room_x))*max(abs(Sol_bad.u(:))), 'k-', 'LineWidth', 2);
plot3(sov_x, sov_y, ones(size(sov_x))*max(abs(Sol_bad.u(:))), 'r-', 'LineWidth', 2);
plot3(worst_xs, yTV_fixed, max(abs(Sol_bad.u(:))), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
title(sprintf('Dålig placering: (%.2f, %.2f), A = %.4f', worst_xs, yTV_fixed, A_values(worst_idx)), 'FontSize', 12);

%% Visualisering 3: 3D perspektiv med belysning
% Skapa 3D-visualiseringar med belysning för att tydligt visa ljudvågornas form

figure('Position', [100, 100, 1200, 400]);

% Optimal 2D
subplot(1,3,1);
surf(Sol_opt_2D.x, Sol_opt_2D.y, abs(Sol_opt_2D.u));
shading interp;  % Jämn skuggning
light;  % Lägg till ljuskälla
lighting gouraud;  % Jämn belysning
view(-60, 60);  % Snett perspektiv
colormap('jet');
material shiny;
hold on;
plot3(x_opt(1), x_opt(2), max(abs(Sol_opt_2D.u(:))), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
title(sprintf('Optimal 2D: (%.2f, %.2f), A = %.4f', x_opt(1), x_opt(2), A_opt), 'FontSize', 12);

% Optimal vägg
subplot(1,3,2);
surf(Sol_opt_wall.x, Sol_opt_wall.y, abs(Sol_opt_wall.u));
shading interp;
light;
lighting gouraud;
view(-60, 60);
colormap('jet');
material shiny;
hold on;
plot3(opt_xs, yTV_fixed, max(abs(Sol_opt_wall.u(:))), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
title(sprintf('Optimal på vägg: (%.2f, %.2f), A = %.4f', opt_xs, yTV_fixed, opt_A), 'FontSize', 12);

% Dålig placering
subplot(1,3,3);
surf(Sol_bad.x, Sol_bad.y, abs(Sol_bad.u));
shading interp;
light;
lighting gouraud;
view(-60, 60);
colormap('jet');
material shiny;
hold on;
plot3(worst_xs, yTV_fixed, max(abs(Sol_bad.u(:))), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
title(sprintf('Dålig placering: (%.2f, %.2f), A = %.4f', worst_xs, yTV_fixed, A_values(worst_idx)), 'FontSize', 12);

%% Visualisering 4: Källfunktionen för optimal 2D-placering
% Visualisera själva ljudkällan för att visa dess form

figure;
% Skapa en funktion som visar källan centerad vid optimal position
S_func = @(x,y) S0_handle(sqrt((x-x_opt(1)).^2 + (y-x_opt(2)).^2));
surf(Sol_opt_2D.x, Sol_opt_2D.y, S_func(Sol_opt_2D.x, Sol_opt_2D.y));
shading interp;
light;
lighting gouraud;
colormap('autumn');
material shiny;
view(-60, 60);
title('Källfunktion för optimal 2D-placering', 'FontSize', 14);

%% ------------------- Hjälpfunktioner ------------------------------------

function A_val = evaluateA_x(xs, y_fixed, omega, aTV, S0_handle, N, sov_x_max, sov_y_min, sov_y_max)
% Utvärderar relativa ljudstyrkan A för en källa placerad vid väggen (xs, y_fixed)
% A = max_{Ω_sov}|u| / max_{Ω}|u|, där u erhålls med hhsolver.
%
% Parametrar:
%   xs         - x-position för källan
%   y_fixed    - fast y-position (0.6 för väggen)
%   omega      - frekvens
%   aTV        - källstyrka
%   S0_handle  - funktion som beskriver källans form
%   N          - upplösning för beräkningarna
%   sov_x_max  - högsta x-värde för sovhörnan
%   sov_y_min  - lägsta y-värde för sovhörnan
%   sov_y_max  - högsta y-värde för sovhörnan
%
% Returnerar:
%   A_val      - relativa ljudstyrkan i sovhörnan

    % Definiera källfunktionen centrerad vid (xs, y_fixed)
    S = @(x,y) aTV * S0_handle(sqrt((x - xs).^2 + (y - y_fixed).^2));
    
    % Lös Helmholtz-ekvationen med denna källa
    [~, Sol] = hhsolver(omega, S, N);
    
    % Platta ut rutnätsdata till vektorer för enklare bearbetning
    X = Sol.x(:);  % Alla x-koordinater
    Y = Sol.y(:);  % Alla y-koordinater
    U = Sol.u(:);  % Lösning för alla punkter
    
    % Hitta alla punkter som ligger i sovhörnan
    w = find(X <= sov_x_max & Y >= sov_y_min & Y <= sov_y_max);
    if isempty(w)
        error('Sovhörnans region hittades ej.');
    end
    
    % Beräkna maximala amplituden i sovhörnan och i hela rummet
    A_sov = max(abs(U(w)));  % Max i sovhörnan
    A_tot = max(abs(U));     % Max i hela rummet
    
    % Relativa ljudstyrkan = max i sovhörnan / max i hela rummet
    A_val = A_sov / A_tot;
end

function A_val = evaluateA_xy(x, omega, aTV, S0_handle, N, sov_x_max, sov_y_min, sov_y_max)
    % Robust utvärdering av relativa ljudstyrkan A för en källa placerad i 2D
    % Ger strafftillägg för positioner utanför rummet
    %
    % Parametrar:
    %   x           - vektor med [xs, ys] koordinater för källan
    %   övriga parametrar - samma som för evaluateA_x
    %
    % Returnerar:
    %   A_val       - relativa ljudstyrkan eller straffvärde vid ogiltig position
    
    xs = x(1);  % x-koordinat
    ys = x(2);  % y-koordinat
    
    % Kontrollera om positionen är utanför rummet
    outsideRoom = false;
    violation = 0;
    
    % 1. Kontrollera rumsgränser
    if xs < 0 || xs > 1 || ys < 0 || ys > 0.75
        outsideRoom = true;
        violation = violation + max(0, -xs) + max(0, xs-1) + max(0, -ys) + max(0, ys-0.75);
    end
    
    % 2. Kontrollera L-utskärningen
    if xs > 0.25 && xs < 0.6 && ys > 0.4 && ys < 0.6
        outsideRoom = true;
        violation = violation + 0.1;
    end
    
    % 3. Kontrollera väggen - VIKTIG: ingen position ovanför y=0.6 tillåts
    if ys >= 0.6
        outsideRoom = true;
        violation = violation + 0.2 + (ys - 0.6) * 5;  % Högre straff ju längre över väggen
    end
    
    % Om utanför rummet, returnera straffvärde
    if outsideRoom
        A_val = 1.0 + violation;
        return;
    end
    
    % Om innanför rummet, beräkna faktiska ljudstyrkan
    S = @(x,y) aTV * S0_handle(sqrt((x - xs).^2 + (y - ys).^2));
    [~, Sol] = hhsolver(omega, S, N);
    
    % Bearbeta lösningen
    X = Sol.x(:);
    Y = Sol.y(:);
    U = Sol.u(:);
    
    % Hitta punkter i sovhörnan
    w = find(X <= sov_x_max & Y >= sov_y_min & Y <= sov_y_max);
    if isempty(w)
        error('Sovhörnans region hittades ej.');
    end
    
    % Beräkna relativa ljudstyrkan
    A_sov = max(abs(U(w)));
    A_tot = max(abs(U));
    A_val = A_sov / A_tot;
end

function [opt_x, fval] = goldenSectionSearch(f, a, b, tol)
    % Utför gyllene snitt-sökning för att hitta minimum av en funktion i ett intervall
    %
    % Parametrar:
    %   f       - Funktionshandtag till funktionen som ska minimeras
    %   a, b    - Sökintervall [a, b]
    %   tol     - Tolerans (hur nära vi vill komma optimum)
    %
    % Returnerar:
    %   opt_x   - x-värde där funktionen har sitt minimum
    %   fval    - Funktionsvärdet vid minimum
    
    % Gyllene snittet - gyllene snittvärdet används för att få en effektiv sökning
    phi = (sqrt(5)-1)/2;  % ≈ 0.618
    
    % Initiera de inre punkterna i intervallet
    x1 = b - phi*(b - a);  % Vänstra inre punkten
    x2 = a + phi*(b - a);  % Högra inre punkten
    
    % Beräkna funktionsvärden vid de inre punkterna
    f1 = f(x1);
    f2 = f(x2);
    
    % Skriv ut rubrik för sökprocessen
    iter = 0;
    fprintf('Gyllene snitt-sökning:\n');
    fprintf('%-5s %-10s %-10s %-10s %-10s %-10s\n', 'Iter', 'a', 'b', 'x1', 'x2', 'f(x)');
    
    % Huvudloop: Fortsätt tills intervallet är mindre än toleransen
    while (b - a) > tol
        iter = iter + 1;
        
        % Kontrollera vilken av de inre punkterna som ger lägst funktionsvärde
        if f1 < f2
            % x1 är bättre, så välj nya intervallet [a, x2]
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi*(b - a);  % Beräkna ny vänstra punkt
            f1 = f(x1);  % Beräkna nytt funktionsvärde
            fprintf('%-5d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f\n', iter, a, b, x1, x2, f1);
        else
            % x2 är bättre, så välj nya intervallet [x1, b]
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi*(b - a);  % Beräkna ny högra punkt
            f2 = f(x2);  % Beräkna nytt funktionsvärde
            fprintf('%-5d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f\n', iter, a, b, x1, x2, f2);
        end
    end
    
    % När intervallet är tillräckligt litet, ta mittpunkten som optimum
    opt_x = (a + b) / 2;
    fval = f(opt_x);
    fprintf('Slutligt resultat: x = %.6f, f(x) = %.6f\n', opt_x, fval);
end