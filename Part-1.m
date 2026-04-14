%% Del 3 - Inversa problemet med isotrop källa
clear; close all; clc;

%% Parametrar
omega       = 19;                         % Frekvensen omega
trueParams  = [1, 0.4, 0.3];              % Sanna parametervärden [a_true, x0_true, y0_true]
noiseLevels = logspace(-2,0,10);          % Brusnivåer från 1% till 100% (logaritmiskt fördelade)
R           = 4 / sqrt(omega);            % Radien för integrationsdomänen (kvadratisk [-R,R]x[-R,R])
numPts      = 1000;                       % Antal punkter per axel för integration (högre = noggrannare)

%% Uppgift 1: Beräkna eta med egen trapetsregel i 2D
% eta = integral2d S0(r)*cos(omega * x) - konstanten som behövs för att lösa inversa problemet
eta = calculateEta2D(@sourceProfile, omega, R, numPts);  % Anropa integrationen
fprintf('Uppgift 1: eta = %.6f\n', eta);  % Skriv ut resultatet med 6 decimaler

%% Uppgift 2: Lös Helmholtz-ekvationen och visualisera u(x,y) och g(s)
[bd, sol] = solveHelmholtz(omega, trueParams(1), trueParams(2), trueParams(3));  % Lös med kända parametrar
visualizeSolution(sol, bd);  % Plotta lösningen u(x,y) och normalderivatan g(s)

%% Uppgift 3: Inversa problemet med Gauss–Newton + Levenberg–Marquardt
aVals = [0.5,1,1.5];                     % Möjliga startvärden för amplitud a
xVals = [0.2,0.4,0.6,0.8];               % Möjliga startvärden för x-koordinat
yVals = [0.2,0.3,0.4,0.5];               % Möjliga startvärden för y-koordinat
% --- 3.1 LM-dämpad GN (Levenberg–Marquardt) ---
[p_LM, err_LM] = invertSourceParameters(bd, eta, omega, aVals, xVals, yVals);
fprintf('Uppgift 3 (LM): a=%.4f, x0=%.4f, y0=%.4f, kv.fel=%.2e\n', p_LM, err_LM);

% --- 3.2 Ren GN (utan dämpning) ---
[p_GN, err_GN] = invertSourceParametersGN(bd, eta, omega, aVals, xVals, yVals);
fprintf('Uppgift 3 (GN): a=%.4f, x0=%.4f, y0=%.4f, kv.fel=%.2e\n\n', p_GN, err_GN);

% --- 3.3 Plotta sida vid sida för att jämföra ---
theta   = linspace(0, 2*pi, 50)';
Ic_data = computeIc(bd, omega, theta);
model   = @(p,t) p(1) * eta .* cos(omega * (p(2)*cos(t) + p(3)*sin(t)));

figure('Name','Del 3: GN vs LM');
subplot(1,2,1);
plot(theta, Ic_data, 'ko', 'MarkerSize', 4, 'DisplayName','Data'); hold on;
plot(theta, model(p_LM, theta), 'r-', 'LineWidth',1.5, 'DisplayName','LM-passning');
plot(theta, model(trueParams,theta), 'g--','LineWidth',1.5, 'DisplayName','Sann'); 
title('LM (Levenberg–Marquardt)'); xlabel('\theta'); ylabel('I_c(\theta)');
legend('Location','best'); grid on;

subplot(1,2,2);
plot(theta, Ic_data, 'ko', 'MarkerSize', 4, 'DisplayName','Data'); hold on;
plot(theta, model(p_GN, theta), 'b-', 'LineWidth',1.5, 'DisplayName','GN-passning');
plot(theta, model(trueParams,theta), 'g--','LineWidth',1.5, 'DisplayName','Sann');
title('GN (utan dämpning)'); xlabel('\theta'); ylabel('I_c(\theta)');
legend('Location','best'); grid on;

sgtitle('Jämförelse: ren GN vs LM');

%% Uppgift 4: Robusthetsanalys mot brus på g(s)
relErr_original = performNoiseRobustness(bd, eta, omega, trueParams, noiseLevels);  % Testa med olika brusnivåer

%% Uppgift 5: Analytiska startgissningar och jämförelse
initialParams = estimateInitialGuess(bd, eta, omega);  % Beräkna analytisk startgissning
[p_analytic, err_analytic] = invertSourceParametersSingle(bd, eta, omega, initialParams);  % Hitta parametrar med analytisk start
plotCompareResults(bd, eta, omega, p_analytic, p_LM, trueParams);  % Jämför resultat från olika startgissningar
relErr_analytic = performNoiseRobustnessWithAnalyticStart(bd, eta, omega, noiseLevels, initialParams, trueParams);  % Robusthet med analytisk start
compareRobustness(noiseLevels, relErr_original, relErr_analytic);  % Jämför robusthet för olika metoder
maxTolOriginal = findMaxTolerance(noiseLevels, relErr_original, 0.2);  % Hitta max tillåten brusnivå för original-metoden
maxTolAnalytic = findMaxTolerance(noiseLevels, relErr_analytic, 0.2);  % Hitta max tillåten brusnivå för analytisk metod
fprintf('Max brusnivå (original): %.1f%%\n', maxTolOriginal*100);  % Skriv ut resultat i procent
fprintf('Max brusnivå (analytisk): %.1f%%\n', maxTolAnalytic*100);  % Skriv ut resultat i procent

%% Uppgift 6: Analysera okända källor från .mat-filer
omegas      = [18, 22, 24, 26, 38]; 
numPtsForEta = 500;  % antal punkter för η-beräkning
results = analyzeUnknownSources('sources', omegas, numPtsForEta);  % Analysera källfiler i mappen 'sources'

%% Hjälpfunktioner

function val = sourceProfile(r)
    % Definierar källans profil som funktion av avstånd r från källans position
    val = cos(24*r) .* exp(-900*r.^2);
end

function eta = calculateEta2D(S0fun, omega, R, numPts)
    % Beräknar eta = integral S0(r)*cos(ωx) dx dy över hela planet med trapetsregeln
    % S0fun: funktionshandtag till källprofilen
    % omega: vinkelfrekvens
    % R: integrationsdomänens radie
    % numPts: antal punkter per axel
    
    x = linspace(-R,R,numPts);           % Skapa x-koordinater från -R till R
    dx = x(2)-x(1);                      % Steglängd i x-led
    [X,Y] = meshgrid(x,x);               % Skapa 2D rutnät av koordinater
    F = S0fun(hypot(X,Y)) .* cos(omega .* X);  % Beräkna integrandfunktionen (hypot = sqrt(x²+y²))
    w = ones(1,numPts);                  % Vikter för trapetsregeln
    w([1,end])=0.5;                      % Halv vikt vid ändpunkterna
    W = w'*w;                            % 2D-matris med vikter (ytterkanterna får halv eller kvarts vikt)
    eta = sum(F(:).*W(:)) * dx^2;        % Utför 2D-trapetsregeln (summera viktade funktionsvärden × area)
end

function [bd, sol] = solveHelmholtz(omega, a, x0, y0)
    % Löser Helmholtz-ekvationen för givna källparametrar
    % omega: vinkelfrekvens
    % a: källans styrka
    % x0, y0: källans position
    
    src = @(x,y) a .* sourceProfile(hypot(x-x0,y-y0));  % Definiera källfunktionen med rätt position och styrka
    N=200;                                              % Antal punkter i beräkningsnätet
    [bd, sol] = hhsolver(omega, src, N);                % Anropa Helmholtz-lösaren
    % bd innehåller randdata, sol innehåller lösningen i rummet
end

function visualizeSolution(sol, bd)
    % Visualiserar lösningen till Helmholtz-ekvationen
    % sol: lösningen i rummet (u)
    % bd: randdata (g = ∂u/∂n)
    
    figure('Name','Uppgift 2');                      % Skapa nytt figurfönster med titel
    
    subplot(1,2,1);                                  % Första delfiguren (1 rad, 2 kolumner, position 1)
    surf(sol.x, sol.y, sol.u);                       % Skapa 3D-yta med lösningen
    view(0,90);                                      % Visa uppifrån (2D-vy)
    shading interp;                                  % Mjuk skuggning
    axis equal;                                      % Lika skalning på axlarna
    title('u(x,y)');                                 % Titel för delfiguren
    xlabel('x');                                     % Etikett för x-axeln
    ylabel('y');                                     % Etikett för y-axeln
    colorbar;                                        % Visa färgskala
    
    subplot(1,2,2);                                  % Andra delfiguren
    plot(bd.s, bd.un,'LineWidth',1.5);               % Plotta normalderivatan mot båglängd
    title('g(s) = ∂u/∂n');                           % Titel för delfiguren
    xlabel('s');                                     % Etikett för x-axeln (båglängd)
    ylabel('g(s)');                                  % Etikett för y-axeln
    grid on;                                         % Visa rutnät
end

function Ic = computeIc(bd, omega, theta)
    % Beräknar cosinus-integralen Ic(α) för alla givna vinklar θ
    % Ic(α) = ∫ g(s)·cos(ω(x·cos(α)+y·sin(α)))ds
    % bd: randdata
    % omega: vinkelfrekvens
    % theta: vinklar för vilka Ic ska beräknas
    
    s=bd.s(:);                               % Båglängd längs randen (kolonnvektor)
    x=bd.x(:);                               % x-koordinater på randen (kolonnvektor)
    y=bd.y(:);                               % y-koordinater på randen (kolonnvektor)
    g=bd.un(:);                              % Normalderivatan på randen (kolonnvektor)
    
    ds=[s(2)-s(1);diff(s)];                  % Steglängder längs randen (första steget uppskattas)
    w=ones(size(s));                         % Vikter för trapetsregeln
    w([1,end])=0.5;                          % Halv vikt vid ändpunkterna
    
    % Beräkna Ic för varje vinkel theta med arrayfun
    Ic = arrayfun(@(t) sum(g .* cos(omega*(x*cos(t)+y*sin(t))) .* ds .* w), theta);
    % För varje t: summera g·cos(ω(x·cos(t)+y·sin(t)))·ds·w över alla randpunkter
end

function Is = computeIs_bd(bd, omega, theta)
    % Beräknar sinus-integralen Is(α) för alla givna vinklar θ
    % Is(α) = ∫ g(s)·sin(ω(x·cos(α)+y·sin(α)))ds
    % bd: randdata
    % omega: vinkelfrekvens
    % theta: vinklar för vilka Is ska beräknas
    
    s=bd.s(:);                               % Båglängd längs randen (kolonnvektor)
    x=bd.x(:);                               % x-koordinater på randen (kolonnvektor)
    y=bd.y(:);                               % y-koordinater på randen (kolonnvektor)
    g=bd.un(:);                              % Normalderivatan på randen (kolonnvektor)
    
    ds=[s(2)-s(1);diff(s)];                  % Steglängder längs randen (första steget uppskattas)
    w=ones(size(s));                         % Vikter för trapetsregeln
    w([1,end])=0.5;                          % Halv vikt vid ändpunkterna
    
    % Beräkna Is för varje vinkel theta med arrayfun
    Is = arrayfun(@(t) sum(g .* sin(omega*(x*cos(t)+y*sin(t))) .* ds .* w), theta);
    % För varje t: summera g·sin(ω(x·cos(t)+y·sin(t)))·ds·w över alla randpunkter
end

function [p, err] = invertSourceParameters(bd, eta, omega, aVals, xVals, yVals)
    % Löser det inversa problemet för att hitta källans parametrar
    % bd: randdata
    % eta: konstanten beräknad i uppgift 1
    % omega: vinkelfrekvens
    % aVals, xVals, yVals: vektorer med möjliga startvärden för källparametrarna
    
    theta=linspace(0,2*pi,30)';            % 30 vinklar från 0 till 2π (kolonnvektor)
    Ic=computeIc(bd,omega,theta);           % Beräkna Ic för alla vinklar
    
    % Modellen som ska anpassas till data: a·η·cos(ω(x0·cos(θ)+y0·sin(θ)))
    model=@(p,t) p(1)*eta .* cos(omega*(p(2)*cos(t)+p(3)*sin(t)));
    
    bestErr=inf;                           % Initiera bästa felet till oändligheten
    p=[nan nan nan];                       % Initiera bästa parametrar till NaN
    
    % Testa alla kombinationer av startvärden
    for A=aVals                            % För varje möjligt a-värde
        for X=xVals                         % För varje möjligt x0-värde
            for Y=yVals                      % För varje möjligt y0-värde
                [pTry,eTry]=dampedGaussNewton([A,X,Y],Ic,theta,eta,omega,1e-3,model);  % Kör optimering med dessa startvärden
                if eTry<bestErr               % Om resultatet är bättre än tidigare bästa
                    bestErr=eTry;              % Spara det nya bästa felet
                    p=pTry;                    % Spara de nya bästa parametrarna
                end
            end
        end
    end
    err=bestErr;                          % Returnera det bästa felet
end

function [p, err] = invertSourceParametersGN(bd, eta, omega, aVals, xVals, yVals)
    theta   = linspace(0, 2*pi, 30)';
    Ic_data = computeIc(bd, omega, theta);
    model   = @(p,t) p(1)*eta .* cos(omega*(p(2)*cos(t) + p(3)*sin(t)));
    err_best = inf; p = [nan nan nan];
    for A = aVals
        for X = xVals
            for Y = yVals
                [pTry, eTry] = gaussNewton([A, X, Y], Ic_data, theta, eta, omega, model);
                if eTry < err_best
                    err_best = eTry; p = pTry;
                end
            end
        end
    end
    err = err_best;
end


function [p, err] = dampedGaussNewton(p0, Ic, theta, eta, omega, lambda, model)
    % Implementerar Gauss-Newton-metoden med dämpning (Levenberg-Marquardt)
    % p0: startgissning för parametrarna [a, x0, y0]
    % Ic: cosinus-integralen för olika vinklar
    % theta: vinklar för vilka Ic beräknats
    % eta: konstanten från uppgift 1
    % omega: vinkelfrekvens
    % lambda: dämpningsparameter
    % model: modellen som ska anpassas till data
    
    p=p0;                                    % Börja med startgissningen
    
    % Iterera tills konvergens eller max 100 iterationer
    for i=1:100
        r = Ic-model(p,theta);                % Beräkna residualen (skillnad mellan data och modell)
        J = jacobianGN(p,theta,eta,omega);    % Beräkna Jacobi-matrisen (derivata av modellen)
        
        % Beräkna normalekvationerna med Levenberg-Marquardt-dämpning
        H = J'*J + lambda*diag(diag(J'*J));   % H: Hessian-matrisen med dämpning på diagonalen
        dp=H\(J'*r);                          % Lös linjära systemet för parameteruppdateringen
        
        pNew=p+dp';                           % Beräkna nya parametrar
        rNew=Ic-model(pNew,theta);            % Beräkna ny residual
        
        % Uppdatering av dämpningsparametern
        if sum(rNew.^2)<sum(r.^2)             % Om de nya parametrarna ger mindre residual
            p=pNew;                            % Acceptera uppdateringen
            lambda=lambda/10;                  % Minska dämpningen (mer Gauss-Newton-lik)
        else                                  % Annars (uppdateringen försämrar resultatet)
            lambda=lambda*10;                  % Öka dämpningen (mer gradient descent-lik)
        end
        
        if norm(dp)<1e-10                    % Konvergenstest: om uppdateringen är väldigt liten
            break;                             % Avbryt iterationen
        end
    end
    
    err=sum((Ic-model(p,theta)).^2);         % Beräkna slutligt kvadratfel
end

function [p, err] = gaussNewton(p0, Ic, theta, eta, omega, model)
    p = p0;
    for i = 1:100
        r = Ic - model(p, theta);
        a  = p(1); x0 = p(2); y0 = p(3);
        T  = omega * (x0*cos(theta) + y0*sin(theta));
        J  = [eta*cos(T), ...
              -a*eta*omega*sin(T).*cos(theta), ...
              -a*eta*omega*sin(T).*sin(theta)];
        dp = (J' * J) \ (J' * r);
        pNew = p + dp';
        if norm(dp) < 1e-10
            p = pNew; break;
        end
        p = pNew;
    end
    err = sum((Ic - model(p, theta)).^2);
end


function J = jacobianGN(p, theta, eta, omega)
    % Beräknar Jacobi-matrisen för model = a·η·cos(ω(x0·cos(θ)+y0·sin(θ)))
    % Varje rad är derivatan med avseende på [a, x0, y0]
    % p: aktuella parametervärden [a, x0, y0]
    % theta: vinklar för vilka modellen evalueras
    % eta: konstanten från uppgift 1
    % omega: vinkelfrekvens
    
    a=p(1);                                      % Extrahera parametervärden
    x0=p(2);
    y0=p(3);
    
    T=omega*(x0*cos(theta)+y0*sin(theta));       % Hjälpvariabel för att förenkla uttrycken
    
    % Jacobi-matrisen med partiella derivator
    % Viktigt: Skapa en matris med rätt dimensioner (antal_observationer x antal_parametrar)
    % Varje kolumn är derivatan med avseende på en parameter
    J = zeros(length(theta), 3);                  % Initiera Jacobi-matrisen med rätt storlek
    
    % Fyll i Jacobi-matrisen kolumn för kolumn
    J(:,1) = eta*cos(T);                          % ∂model/∂a = η·cos(ω(x0·cos(θ)+y0·sin(θ)))
    J(:,2) = -a*eta*omega*sin(T).*cos(theta);     % ∂model/∂x0 = -a·η·ω·sin(ω(x0·cos(θ)+y0·sin(θ)))·cos(θ)
    J(:,3) = -a*eta*omega*sin(T).*sin(theta);     % ∂model/∂y0 = -a·η·ω·sin(ω(x0·cos(θ)+y0·sin(θ)))·sin(θ)
end

function init = estimateInitialGuess(bd, eta, omega)
    % Uppskattar startgissning för källans parametrar med analytiska formler
    % bd: randdata
    % eta: konstanten från uppgift 1
    % omega: vinkelfrekvens
    
    theta0=0;                                    % Vinkel 0 (i x-riktningen)
    theta90=pi/2;                                % Vinkel π/2 (i y-riktningen)
    h=1e-3;                                      % Litet steg för numerisk derivering
    
    % Beräkna Ic och Is vid specifika vinklar
    Ic0=computeIc(bd,omega,theta0);              % Ic vid θ=0
    Is0=computeIs_bd(bd,omega,theta0);           % Is vid θ=0
    Ic90=computeIc(bd,omega,theta90);            % Ic vid θ=π/2
    Is90=computeIs_bd(bd,omega,theta90);         % Is vid θ=π/2
    
    % Beräkna derivatan av Ic med central differens
    dIc0=(computeIc(bd,omega,theta0+h)-computeIc(bd,omega,theta0-h))/(2*h);  % dIc/dθ vid θ=0
    dIc90=(computeIc(bd,omega,theta90+h)-computeIc(bd,omega,theta90-h))/(2*h); % dIc/dθ vid θ=π/2
    
    % Beräkna x0 och y0 med analytiska formler från uppgift 5
    x0=dIc90/(omega*Is90);                       % x0 = (dIc/dθ vid θ=π/2)/(ω·Is(π/2))
    y0=-dIc0/(omega*Is0);                        % y0 = -(dIc/dθ vid θ=0)/(ω·Is(0))
    
    % Beräkna amplituden
    aVal=sqrt(Ic0.^2+Is0.^2)/eta;                % a = √(Ic(0)²+Is(0)²)/η
    
    init=[aVal,x0,y0];                           % Sammanställ startgissningen
end


function relErr = performNoiseRobustness(bd, eta, omega, trueP, noiseLevels)
    % bd: struct med bd.s, bd.x, bd.y, bd.un
    % eta, omega: som tidigare
    % trueP = [a_true, x0_true, y0_true]
    % noiseLevels: vektor med brusnivåer t.ex. logspace(-2,0,10)
    
    % Hämta data från bd
    s  = bd.s(:);
    x  = bd.x(:);
    y  = bd.y(:);
    g0 = bd.un(:);                  % Ursprungs-g(s)
    ds = [s(2)-s(1); diff(s)];      % Trapets-segment
    w  = ones(size(s)); w([1,end]) = 0.5;
    
    % Bestäm vinklar för Ic
    Mtheta = 20;
    thetas = linspace(0,2*pi,Mtheta)';
    
   
    relErr = nan(3, numel(noiseLevels));
    
    % Välj tre brusnivåer för illustration
    exampleLevels = [noiseLevels(1), noiseLevels(round(end/2)), noiseLevels(end)];
    
    % --- 1) Plotta g och g_noisy för några nivåer ---
    figure('Name','Brus­exempel på g(s)');
    for i = 1:3
        nl = exampleLevels(i);
        gN = g0 + max(abs(g0))*randn(size(g0))*nl;
        
        subplot(3,1,i);
        plot(s, g0, 'LineWidth',1);
        hold on;
        plot(s, gN, '--','LineWidth',1);
        ylabel('g(s)');
        title(sprintf('Brusnivå = %.2f%%', nl*100));
        if i==3, xlabel('s (båglängd)'); end
        legend('g','g\_noisy','Location','Best');
        grid on;
    end
    
    % --- 2) Robusthetsanalys över alla brusnivåer ---
    % Modellfunktion för I_c
    model = @(p,th) p(1)*eta .* cos( omega*(p(2)*cos(th)+p(3)*sin(th)) );
    
    for k = 1:numel(noiseLevels)
        nl = noiseLevels(k);
        gN = g0 + max(abs(g0))*randn(size(g0))*nl;
        
        % Beräkna Ic med brus
        IcN = arrayfun(@(t) sum( gN .* cos(omega*(x*cos(t)+y*sin(t))) .* ds .* w ), thetas);
        
        
        obj = @(p) sum( (IcN - model(p, thetas)).^2 );
        pEst = fminsearch(obj, trueP);
        
        % Spara relativt fel
        relErr(:,k) = abs(pEst - trueP)' ./ trueP';
    end
    
    % --- 3) Plotta relativt fel vs brusnivå ---
    figure('Name','Robusthetsanalys mot brus');
    paramNames = {'a','x_0','y_0'};
    for i = 1:3
        subplot(1,3,i);
        loglog(noiseLevels*100, relErr(i,:)*100, '-o', 'LineWidth',1.2);
        xlabel('Brusnivå (%)');
        ylabel('Relativt fel (%)');
        title(paramNames{i});
        grid on;
    end
    sgtitle('Relativt fel som funktion av brusnivå');
end


function relErr = performNoiseRobustnessWithAnalyticStart(bd, eta, omega, noiseL, initP, trueP)
    % Undersöker algoritmens robusthet mot brus med analytisk startgissning
    % bd: originalranddata utan brus
    % eta: konstanten från uppgift 1
    % omega: vinkelfrekvens
    % noiseL: vektor med brusnivåer att testa
    % initP: initiala parametervärden (används ej direkt)
    % trueP: sanna parametervärden [a_true, x0_true, y0_true]
    
    s=bd.s(:);                                     % Extrahera randdata
    x=bd.x(:);
    y=bd.y(:);
    g0=bd.un(:);                                   % Original g utan brus
    
    ds=[s(2)-s(1);diff(s)];                        % Steglängder längs randen
    w=ones(size(s));                               % Vikter för trapetsregeln
    w([1,end])=0.5;                                % Halv vikt vid ändpunkterna
    
    theta=linspace(0,2*pi,30)';                    % 30 vinklar från 0 till 2π
    relErr=zeros(3,numel(noiseL));                 % Matris för att lagra relativa fel
    
    % Testa för varje brusnivå
    for k=1:numel(noiseL)
        % Lägg till brus på g
        gN=g0+max(abs(g0))*randn(size(g0))*noiseL(k);
        
        bd2=bd;                                    % Skapa en kopia av randdata
        bd2.un=gN;                                 % Ersätt g med brusig version
        
        % Beräkna analytisk startgissning baserat på brusig data
        init2=estimateInitialGuess(bd2,eta,omega);
        
        % Beräkna Ic med brusig data
        IcN=arrayfun(@(t)sum(gN.*cos(omega*(x*cos(t)+y*sin(t))).*ds.*w), theta);
        
        % Hitta parametrar med dampedGaussNewton och analytisk startgissning
        [pEst,~]=dampedGaussNewton(init2,IcN,theta,eta,omega,1e-3,@(p,t)p(1)*eta.*cos(omega*(p(2)*cos(t)+p(3)*sin(t))));
        
        % Beräkna relativt fel för varje parameter
        relErr(:,k)=abs(pEst-trueP)./trueP;
    end
end

function plotInverseFit(bd, eta, omega, estP, trueP)
    % Plottar anpassningen av modellen till data
    % bd: randdata
    % eta: konstanten från uppgift 1
    % omega: vinkelfrekvens
    % estP: estimerade parametrar [a_est, x0_est, y0_est]
    % trueP: sanna parametrar [a_true, x0_true, y0_true]
    
    theta=linspace(0,2*pi,50)';                    % 50 vinklar från 0 till 2π för plottning
    Ic=computeIc(bd,omega,theta);                  % Beräkna Ic för dessa vinklar
    
    % Modellen som anpassats
    model=@(p,t)p(1)*eta.*cos(omega*(p(2)*cos(t)+p(3)*sin(t)));
    
    figure;                                        % Skapa ny figur
    plot(theta,Ic,'ko');                           % Plotta data som svarta cirklar
    hold on;                                       % Behåll figuren för fler plottar
    plot(theta,model(estP,theta),'r-');            % Plotta anpassad modell som röd linje
    plot(theta,model(trueP,theta),'g--');          % Plotta sann modell som grön streckad linje
    legend('Data','Estimerad','Sann');             % Förklaring till de olika linjerna
    xlabel('t');                                   % Etikett för x-axeln
    ylabel('I_c');                                 % Etikett för y-axeln
    grid on;                                       % Visa rutnät
end

function plotCompareResults(bd, eta, omega, p1, p2, trueP)
    % Jämför två olika anpassningar av modellen
    % bd: randdata
    % eta: konstanten från uppgift 1
    % omega: vinkelfrekvens
    % p1: parametrar från metod 1 (analytisk startgissning)
    % p2: parametrar från metod 2 (grid av startgissningar)
    % trueP: sanna parametrar [a_true, x0_true, y0_true]
    
    theta=linspace(0,2*pi,50)';                    % 50 vinklar från 0 till 2π för plottning
    Ic=computeIc(bd,omega,theta);                  % Beräkna Ic för dessa vinklar
    
    % Modellen som anpassats
    model=@(p,t)p(1)*eta.*cos(omega*(p(2)*cos(t)+p(3)*sin(t)));
    
    figure;                                        % Skapa ny figur
    plot(theta,Ic,'k.');                           % Plotta data som svarta punkter
    hold on;                                       % Behåll figuren för fler plottar
    plot(theta,model(p1,theta),'r-','LineWidth',1.5); % Plotta anpassning från metod 1 (röd, tjock)
    plot(theta,model(p2,theta),'b-','LineWidth',1.5); % Plotta anpassning från metod 2 (blå, tjock)
    plot(theta,model(trueP,theta),'g--');          % Plotta sann modell (grön, streckad)
    legend('Data','Analytisk','Grid','Sann');      % Förklaring till de olika linjerna
    xlabel('t');                                   % Etikett för x-axeln
    ylabel('I_c');                                 % Etikett för y-axeln
    grid on;                                       % Visa rutnät
end

function compareRobustness(noiseL, r1, r2)
    % Jämför robustheten hos de två metoderna mot brus
    % noiseL: vektor med brusnivåer
    % r1: relativa fel för metod 1 (grid av startgissningar)
    % r2: relativa fel för metod 2 (analytisk startgissning)
    
    figure;                                       % Skapa ny figur
    for i=1:3                                     % Loopa över de tre parametrarna
        subplot(1,3,i);                           % Skapa delfigur (1 rad, 3 kolumner)
        loglog(noiseL*100,r1(i,:)*100,'b-o');     % Plotta resultat från metod 1 (blå)
        hold on;                                  % Behåll figuren för fler plottar
        loglog(noiseL*100,r2(i,:)*100,'r-s');     % Plotta resultat från metod 2 (röd)
        
        % Sätt etikett för delfiguren beroende på parameter
        if i == 1
            title('Parameter a');                  % Titel för parameter a
        elseif i == 2
            title('Parameter x0');                 % Titel för parameter x0
        else
            title('Parameter y0');                 % Titel för parameter y0
        end
        
        xlabel('Brus%');                          % Etikett för x-axeln (brusnivå i procent)
        ylabel('Fel%');                           % Etikett för y-axeln (fel i procent)
        grid on;                                  % Visa rutnät
    end
    
    legend('Grid','Analytisk');                   % Förklaring till de olika linjerna
    sgtitle('Robusthet mot brus');                % Övergripande titel för hela figuren
end

function maxTol = findMaxTolerance(noiseL, relErr, tol)
    % Hittar den maximala brusnivån där det relativa felet är under en viss tolerans
    % noiseL: vektor med brusnivåer
    % relErr: matris med relativa fel (3 x antal brusnivåer)
    % tol: tolerans för relativt fel
    
    % Hitta maximal brusnivå för varje parameter där felet är under toleransen
    vals = arrayfun(@(i) max(noiseL(relErr(i,:)<tol)),1:3);
    
    % Ta minimum av dessa (begränsande parameter)
    maxTol=min(vals);
end

function results = analyzeUnknownSources(folder, omegas, numPts)
    % Analyserar okända källor från .mat-filer
    % folder: mapp där källfilerna finns
    % omega: vinkelfrekvens
    % eta: konstanten från uppgift 1
    % initP: initial gissning för parametrarna
    
   files = dir(fullfile(folder, '*.mat'));
    nFiles = numel(files);

    % Kontroll: antalet omegas måste stämma överens med antal filer
    if nFiles ~= numel(omegas)
        error('Antalet frekvenser (%d) matchar inte antalet filer (%d).', ...
              numel(omegas), nFiles);
    end

    results = cell(nFiles, 3);

    for k = 1:nFiles
        % --- a) Läs in filen och hitta strukturen 'bd' ---
        data = load(fullfile(files(k).folder, files(k).name));
        fn   = fieldnames(data);
        bd   = [];
        for f = fn'
            v = data.(f{1});
            if isstruct(v) && all(isfield(v, {'s','x','y','un'}))
                bd = v;
                break;
            end
        end
        if isempty(bd)
            error('Kunde inte hitta en struktur med fält {s, x, y, un} i %s.', files(k).name);
        end

        % --- b) Hämta omega_k från vektorn 'omegas' ---
        omega_k = omegas(k);

        % --- c) Beräkna trunkeringsradie R_k och η_k för detta ω ---
        R_k   = 4 / sqrt(omega_k);
        eta_k = calculateEta2D(@sourceProfile, omega_k, R_k, numPts);

        % --- d) Skatta [a, x0, y0] med LM (dampedGaussNewton) ---
        initP = estimateInitialGuess(bd, eta_k, omega_k);
        [pEst, err] = invertSourceParametersSingle(bd, eta_k, omega_k, initP);

        % --- e) Spara resultatet i cell-arrayen ---
        results{k,1} = files(k).name;   % Filnamn (t.ex. 'Källa1.mat')
        results{k,2} = pEst;            % Skattade parametrar [a_est, x0_est, y0_est]
        results{k,3} = err;             % Kvadratsummafel
    end

    % --- f) Skriv ut sammanställning i kommandofönstret ---
    disp('Del 6 – Okända källor: Fil | [a, x0, y0] | kvadratsummafel');
    for k = 1:nFiles
        fn = results{k,1};
        p = results{k,2};
        e = results{k,3};
        fprintf('%-12s   [%.4f, %.4f, %.4f]   %.2e\n', fn, p(1), p(2), p(3), e);
    end
end
   

function [p,err]=invertSourceParametersSingle(bd,eta,omega,startP)
    % Inverterar källans parametrar med en given startgissning
    % bd: randdata
    % eta: konstanten från uppgift 1
    % omega: vinkelfrekvens
    % startP: startgissning för parametrarna
    
    theta=linspace(0,2*pi,30)';             % 30 vinklar från 0 till 2π
    Ic=computeIc(bd,omega,theta);           % Beräkna Ic för alla vinklar
    
    % Modellen som ska anpassas till data
    model=@(p,t)p(1)*eta.*cos(omega*(p(2)*cos(t)+p(3)*sin(t)));
    
    % Använd dampedGaussNewton för att hitta parametrarna
    [p,err]=dampedGaussNewton(startP,Ic,theta,eta,omega,1e-3,model);
end