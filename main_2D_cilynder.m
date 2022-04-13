clc; clear; close all; more off
tic

Robin = 1; % Choose model with Robin or Neumann boundary conditions at cylinder base
snapshots = 1; % Print temperature profile every time

%---------------------------------------------------------------
% INPUTS
%---------------------------------------------------------------
t_horas = 168; % Número de horas para a simulação [ h ]
t_final = t_horas * 3600; % Tempo de estabilização [ s ]
Nr = 30; % Número de discretizações radiais
Nz = 30; % Número de discretizações axiais
T0 = 80 + 273.15; % Temperatura inicial da cortiça [ K ]
Tinf = 20 + 273.15; % Temperatura ambiente [ K ]
h = 10; % Coeficiente de transferÊncia de calor [ W / (m2.K) ]
%---------------------------------------------------------------

% Propriedades do cilindro de cortiça

R = 0.55; % Raio externo do cilindro [ m ]
L = 0.92; % Comprimento do cilindro [ m ]
k = 0.038; % Condutividade térmica da cortiça [ W / (m.k) ]
rho = 155; % Massa volúmica da cortiça [ kg / m3]
cp = 1900; % Capacidade térmica da cortiça [ J / (kg.K) ]
alpha = k / rho / cp; % Difusividade térmica da cortiça [ m2 / s ] (artigo = 1e-6)

% Parâmetros do modelo

Bir = h * R / k; % Número de Biot radial [ - ]
For = alpha * t_final / R^2; % Número de Fourier radial [ - ]
Biz = h * L / k; % Número de Biot axial [ - ]
Foz = alpha * t_final / L^2; % Número de Fourier axial [ - ]

par.Bir = Bir;
par.Biz = Biz;
par.For = For;
par.Foz = Foz;

% Geometrias do modelo

r = linspace(0, R, Nr) / R; % Número de discretizações radiais
dr = (r(2) - r(1)); % dr normalizado
dr2 = dr^2; % dr2 normalizado

z = linspace(0, L, Nz) / L; % Número de discretizações axiais
dz = (z(2) - z(1)); % dz normalizado
dz2 = dz^2; % dz2 normalizado

geo.r = r;
geo.dr = dr;
geo.dr2 = dr2;
geo.dz = dz;
geo.dz2 = dz2;


t_span = [0 t_final] / t_final; % Número de discretizações temporais

% Condições iniciais normalizadas
Nv = Nr * Nz;
T0_norm = ones(Nv, 1); % (T - Tinf) / (T0 - Tinf);

% Resolução do sistema de ODEs
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

if Robin == 1
    % Modelo com transferência de calor em todas as superfícies
    rhs = @(t, x) fun_model_Robin(t, x, par, geo, Nz, Nr);
else
    % Modelo sem transferência de calor na base
    rhs = @(t, x) fun_model_Neumann(t, x, par, geo, Nz, Nr);
end

[res_t, res_T] = ode15s(rhs, t_span, T0_norm, options);
toc

% Reformular o vetor solução obtido
kmax = length(res_t);
T_sol = zeros(Nr, Nz, kmax);
% Temperaturas para todos os tempos
% Uma snapshot do perfil de Temperaturas ao longo do tempo
for k = 1:kmax
    T_sol(:, :, k) = reshape(res_T(k, :), [Nr, Nz]);
end


%-------------------------------------------------------
% Gráfico do último instante (duplicado)
%-------------------------------------------------------
T_last = T_sol(:, :, end) * (T0 - Tinf) + Tinf - 273.15;
T_all = [flip(T_last(2:end, :)); T_last];
r_neg = -r(2:end);
r_all = [flip(r_neg) r];

figure('color', [1 1 1]);
mesh(z * L, r_all * R, T_all)

if T0 < Tinf
    axis([0 L -R R 20 Tinf - 273.15])
else
    axis([0 L -R R Tinf - 273.15 T0 - 273.15])
end

xlabel('\it z \rm [m]');
ylabel('\it r \rm [m]');
zlabel('\it T(z,r) \rm [^oC]');
% mytitle = sprintf('%s%7.1f%s', 't =', res_t(end) * t_final / 3600, ' h');
% title(mytitle)
view(148, 55)
camlight right;
lighting gouraud;

%----------------------------------------------------------%
if snapshots == 1
    % Snapshot ao longo do tempo

    r_neg = -r(2:end);
    r_all = [flip(r_neg) r];

    for i = 1:length(res_t)

        T_last = T_sol(:, :, i) * (T0 - Tinf) + Tinf - 273.15;
        T_all = [flip(T_last(2:end, :)); T_last];
        figure('color', [1 1 1]);
        mesh(z * L, r_all * R, T_all)
        axis([0 L -R R 20 80])
        xlabel('\it z \rm [m]');
        ylabel('\it r \rm [m]');
        zlabel('\it T(z,r) \rm [^oC]');
        mytitle = sprintf('%s%7.2f%s', 't =', res_t(i) * t_final / 3600, ' h');
        title(mytitle)
        print(['images/2plt_', num2str(j), 'cyl_t=', num2str(res_t(i) * t_final / 60), '.png'], '-dpng');
        view(148, 55)
        camlight right;
        lighting gouraud;
        pause(0.02)
    end

else
end

