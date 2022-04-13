function [dxdt] = fun_model_Neumann(t, x, par, geo, Nz, Nr)

    %----------------------------------------------------------------------------------------------%
    % Parâmetros do modelo
    %----------------------------------------------------------------------------------------------%
    Bir = par.Bir;
    Biz = par.Biz;
    For = par.For;
    Foz = par.Foz;

    %----------------------------------------------------------------------------------------------%
    % Geometrias do modelo
    %----------------------------------------------------------------------------------------------%

    r = geo.r;
    dr = geo.dr;
    dr2 = geo.dr2;
    dz = geo.dz;
    dz2 = geo.dz2;

    %----------------------------------------------------------------------------------------------%
    % Variáveis de estado
    %----------------------------------------------------------------------------------------------%
    % Ordenar o vetor de Nr * Nz odes para um formato mais percetivel, uma grelha de pontos em r e z
    T = reshape(x, Nr, Nz);

    % derivadas (pré-alocar memória)
    dTdt = zeros(Nr, Nz);

    %----------------------------------------------------------------------------------------------%
    % Criação do sistema de ODEs
    %----------------------------------------------------------------------------------------------%

    j = 1; % indice do z, quando z = 0

    i = 1; % indice do r, quando r = 0

    if r(i) == 0
        % Regra de L’Hospital, devido à indeterminação

        % Condições fronteira (CF)

        Tr0 = T(i + 1, j); % CF para r = 0 - perfil simétrico
        %Tz0 = Biz * T(i, j) * 2 * dz + T(i, j + 1); % CF para z=0 - convecção para o exterior
        Tz0 = T(i, j + 1); % CF para z = 0 - sem transferência de calor para o exterior

        % Modelo

        d2T_dr2 = (T(i + 1, j) - 2.0 * T(i, j) + Tr0) / dr2;
        d2T_dz2 = (T(i, j + 1) - 2.0 * T(i, j) + Tz0) / dz2;
        dTdt(i, j) = 2 * For * d2T_dr2 + Foz * d2T_dz2;

    else
        fprintf('Error!!!! r(i)=/=0\n')
    end

    for i = 2:Nr - 1 % quando r = >0 e <R

        % Condições fronteira (CF)
        %Tz0 = Biz * T(i, j) * 2 * dz + T(i, j + 1); % CF para z=0 - convecção para o exterior
        Tz0 = T(i, j + 1); % CF para z = 0 - sem transferência de calor para o exterior

        % Modelo
        d2T_dr2 = (T(i + 1, j) - 2.0 * T(i, j) + T(i - 1, j)) / dr2;
        d2T_dz2 = (T(i, j + 1) - 2.0 * T(i, j) + Tz0) / dz2;
        dT_dr = (T(i + 1, j) - T(i - 1, j)) / (2 * dr);
        dTdt(i, j) = For * (d2T_dr2 + (1 / r(i)) * dT_dr) + Foz * d2T_dz2;
    end

    i = Nr; % quando r = R

    % Condições fronteira (CF)

    Tz0 = T(i, j + 1); % CF para z = 0 - sem transferência de calor para o exterior
    TR =- Bir * T(i, j) * 2 * dr + T(i - 1, j); % CF para r = R - convecção para o exterior

    % Modelo
    d2T_dr2 = (TR - 2.0 * T(i, j) + T(i - 1, j)) / dr2;
    d2T_dz2 = (T(i, j + 1) - 2.0 * T(i, j) + Tz0) / dz2;
    dT_dr = (TR - T(i - 1, j)) / (2 * dr);
    dTdt(i, j) = For * (d2T_dr2 + (1 / r(i)) * dT_dr) + Foz * d2T_dz2;

    % ---------------------------------------------------------------------------------------%

    for j = 2:Nz - 1 % quando z = [z0+ L-]

        i = 1; % indice do r

        if r(i) == 0
            % Regra de L’Hospital

            % Condições fronteira (CF)

            Tr0 = T(i + 1, j); % CF para r = 0

            % Modelo
            d2T_dr2 = (T(i + 1, j) - 2.0 * T(i, j) + Tr0) / dr2;
            d2T_dz2 = (T(i, j + 1) - 2.0 * T(i, j) + T(i, j - 1)) / dz2;
            dTdt(i, j) = 2 * For * d2T_dr2 + Foz * d2T_dz2;

        else
            fprintf('Error!!!! r(i)=/=0\n')
        end

        for i = 2:Nr - 1 % Meio do cilindro

            % Modelo
            d2T_dr2 = (T(i + 1, j) - 2.0 * T(i, j) + T(i - 1, j)) / dr2;
            d2T_dz2 = (T(i, j + 1) - 2.0 * T(i, j) + T(i, j - 1)) / dz2;
            dT_dr = (T(i + 1, j) - T(i - 1, j)) / (2 * dr);
            dTdt(i, j) = For * (d2T_dr2 + (1 / r(i)) * dT_dr) + Foz * d2T_dz2;
        end

        i = Nr; % r = R

        % Condições fronteira (CF)
        TR =- Bir * T(i, j) * 2 * dr + T(i - 1, j); % CF para r = R - convecção para o exterior

        % Modelo
        d2T_dr2 = (TR - 2.0 * T(i, j) + T(i - 1, j)) / dr2;
        d2T_dz2 = (T(i, j + 1) - 2.0 * T(i, j) + T(i, j - 1)) / dz2;
        dT_dr = (TR - T(i - 1, j)) / (2 * dr);
        dTdt(i, j) = For * (d2T_dr2 + (1 / r(i)) * dT_dr) + Foz * d2T_dz2;

        % ---------------------------------------------------------------------------------------%
    end

    j = Nz; %quando z = L

    i = 1; % indice do r, para r = 0

    if r(i) == 0
        % Regra de L’Hospital

        % Condições fronteira (CF)
        Tr0 = T(i + 1, j); % CF para r = 0 - perfil simétrico
        TL =- Biz * T(i, j) * 2 * dz + T(i, j - 1); % CF para z = L - convecção para o exterior

        % Modelo
        d2T_dr2 = (T(i + 1, j) - 2.0 * T(i, j) + Tr0) / dr2;
        d2T_dz2 = (TL - 2.0 * T(i, j) + T(i, j - 1)) / dz2;
        dTdt(i, j) = 2 * For * d2T_dr2 + Foz * d2T_dz2;

    else
        fprintf('Error!!! r(i)=/=0\n')
    end

    for i = 2:Nr - 1

        % Condições fronteira (CF)
        TL =- Biz * T(i, j) * 2 * dz + T(i, j - 1); % CF para z = L - convecção para o exterior

        % Modelo
        d2T_dr2 = (T(i + 1, j) - 2.0 * T(i, j) + T(i - 1, j)) / dr2;
        d2T_dz2 = (TL - 2.0 * T(i, j) + T(i, j - 1)) / dz2;
        dT_dr = (T(i + 1, j) - T(i - 1, j)) / (2 * dr);
        dTdt(i, j) = For * (d2T_dr2 + (1 / r(i)) * dT_dr) + Foz * d2T_dz2;
    end

    i = Nr; % r = R

    % Condições fronteira (CF)
    TL =- Biz * T(i, j) * 2 * dz + T(i, j - 1); % CF para z = L - convecção para o exterior
    TR =- Bir * T(i, j) * 2 * dr + T(i - 1, j); % CF para r = R - convecção para o exterior

    % Modelo
    d2T_dr2 = (TR - 2.0 * T(i, j) + T(i - 1, j)) / dr2;
    d2T_dz2 = (TL - 2.0 * T(i, j) + T(i, j - 1)) / dz2;
    dT_dr = (TR - T(i - 1, j)) / (2 * dr);
    dTdt(i, j) = For * (d2T_dr2 + (1 / r(i)) * dT_dr) + Foz * d2T_dz2;

    % ---------------------------------------------------------------------------------------%

    % reordenar as derivadas para o formato original
    dxdt = reshape(dTdt, [Nr * Nz, 1]);

end
