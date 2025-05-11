%% Minimering av funktionen F(x, R)
clear all
close all

%%
% Startkonfigurationen av x
s1 = [0.3 0.3 0.3 0.3 0.3 0.3 0.2 0.2 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3];
s2 = [-0.3 -0.3 -0.3 -0.3 -0.3 -0.3 0.2 0.2 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3];
s_vec = [s1, s2];
sl = 65 * ones(1, 21);
x_start = [s_vec, sl];

% Konstant vektor R
R = [20000 10000 9000 8000 7000 6000 5000 4500 4000 3500 3000 2500 2000 1900 1800 1700 1600 1500 1400 1300 1200];

%%
%% Staggers and Span lengths

% Definiera funktionen F som ska minimeras
F = @(x) myF(x, R); % myF är en wrapper för din funktion F(x, R)

% %  - - - - - - -  - -- - - - GA- - - - - - - - - - - - - - - - -
% options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter');
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% % - - - - - - - - Pattern Search - - - - - - - - - - - - - - - - - - --
options = optimoptions('patternsearch', 'MaxIterations', 200, 'Display', 'iter', 'MeshTolerance', 1e-6);
% % -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- 

% Definiera eventuella begränsningar (bounds och linjära begränsningar)
N = numel(x_start);
LB = [-ones(1, 42)*0.4, 65 * ones(1, 21)]; % Nedre gräns för x
UB = [ones(1, 42)*0.4, Inf * ones(1, 21)];  % Övre gräns för x % Ändra till inf

%- - - -  Kör optimeringen --  -- - - - - - - 
% - - - - - - - - - - -Pattern-Search - - - - - -  - - - - - - - - 
[x_opt, fval] = patternsearch(F, x_start, [], [], [], [], LB, UB, options);
% - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - 

% Visa resultat
% fprintf('Optimalt värde av F: %f\n', fval);
% fprintf('Exitflag: %d\n', exitflag);
% disp('Optimal x-konfiguration:');
% disp(x_opt);

HIST_Generator(x_opt, R)

%% Only STAGGER, constant Sl
% Definiera funktionen F som ska minimeras (endast för de första 42 variablerna)
F = @(x) myF([x, 65 * ones(1, 21)], R); % Kombinerar x_opt med de fasta 65-värdena

% Skapa optimeringsalternativ
% 
%  %  - - - - - - - - - - GA- - - - - - - - - - - - - - - - - - - - - -
% options = optimoptions('ga', 'MaxGenerations', 200, 'Display', 'iter');
% % -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 
% % - - - - - - - - Pattern Search - - - - - - - - - - - - - - - - - - --
options = optimoptions('patternsearch', 'MaxIterations', 200, 'Display', 'iter', 'MeshTolerance', 1e-6);
% % -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- 

% - - - - - - - - - - - - PSO - - - - - - - - - - - -- - - - - - - - - - - 
% options = optimoptions('particleswarm', 'SwarmSize', 50, 'MaxIterations', 100, 'Display', 'iter');
% - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Definiera eventuella begränsningar (bounds och linjära begränsningar)
N_opt = 42; % Antal optimeringsvariabler
LB = -ones(1, N_opt)*0.4; % Nedre gräns för de första 42 variablerna
UB = ones(1, N_opt)*0.4;  % Övre gräns för de första 42 variablerna

% Startvärden endast för de första 42 variablerna
x_start_opt = x_start(1:N_opt); % Exempel: x_start definierar startvärdena för alla variabler

% Globala variabler för att lagra konvergenshistorik
global history_fval history_x;
history_fval = [];
history_x = [];

% Kör optimeringen
% [x_opt, fval, exitflag, output] = fmincon(F, x_start_opt, [], [], [], [], LB, UB, [], options);
% 
% % - - - - - - -  - - - - - GA - - - - - -  - - - - - - - - - - - 
% [x_opt, fval] = ga(F, N_opt, [], [], [], [], LB, UB, [], options);
% % - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% % - - - - - - - - - - -Pattern-Search - - - - - -  - - - - - - - - 
[x_opt, fval] = patternsearch(F, x_start_opt, [], [], [], [], LB, UB, options);
% % - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - 

%  - - - - - - - - - - PSO - - - - - - - - - - - - - - - - - -  - - 
% [x_opt, fval] = particleswarm(F, N_opt, LB, UB, options);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Återskapa den fullständiga lösningen (inkluderar de fasta sl värdena)
x_full = [x_opt, 65 * ones(1, 21)];
HIST_Generator(x_full, R)

% Visa resultat
fprintf('Optimalt värde av F: %f\n', fval);
fprintf('Exitflag: %d\n', exitflag);
disp('Optimal x-konfiguration (första 42 optimerade, resterande fasta):');
disp(x_full);

%% Wrapper-funktion för F(x, R)
function f = myF(x, R)
    % Här anropar vi din funktion F(x, R) som definierad i frågan
    %f = L2(x, R);
    term1 = L2(x, R);
    term2 = F_wind_left(x, R);
    term3 = F_wind_right(x, R);
    term4 = Long_spanlengths(x);
    term5 = DRAG_INDICATOR(x,R);
    K_windleft = 1000; %
    K_windright = 1000; %
    K_sl = 0; %
    K_drag = 500; %
    f = term1 + K_windleft * term2 + K_windright * term3 + K_sl * term4 + K_drag*term5; %alla
%     %f = term1 + K1 * term2 + K1 * term3 + K_drag*term5; % bara stagger. 
end

%% Utdatafunktion för att lagra konvergenshistorik
function stop = outfun(x, optimValues, state)
    global history_fval history_x;
    stop = false; % Fortsätt iterera
    switch state
        case 'iter'
            history_fval = [history_fval; optimValues.fval];
            history_x = [history_x; x'];
    end
end
