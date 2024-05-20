M_vec = [2 4 6 8];
syms eta
eta_vec = 0.1:0.1:7;
outdata = [eta_vec.'];

for M = M_vec
    V_dist = abs(PA_model(M));
    C_elec = compute_Ct(V_dist, M);
    fvec = zeros(size(eta_vec));
    for k = 1:numel(eta_vec)
        fvec(k) = subs(C_elec, eta, 10.^eta_vec(k));
    end
    outdata = [outdata fvec.'];
end

writematrix(outdata, "electrical_model.csv");

figure(1);
semilogx(10.^outdata(:,1), outdata(:,2), "--"); hold on; grid on;
semilogx(10.^outdata(:,1), outdata(:,3), "--");
semilogx(10.^outdata(:,1), outdata(:,4), "--");
semilogx(10.^outdata(:,1), outdata(:,5), "--");
xlabel("\eta"); ylabel("C_{elec}");
title("Electrical Crosstalk");
ylim([0 1]);

%% Helper Functions

function [C_elec] = compute_Ct(V_dist, M)
    %% Expression for voltage crosstalk
    
    Vl = [];
    for i = 1:M
        for j = 1:M
            if (i == M | j == M) & ~(i == M & j == M)
                Vl = [Vl; V_dist(i,j)];
            end
        end
    end
    
    C_elec = max(Vl)/V_dist(M,M);
end

function [V_dist] = PA_model(M)
    %% Variables of Interest
    
    syms w_e t_e t_h rho_e l_0 eta
    syms x [1 M] 
    syms y [1 M]
    syms V0 [M M]
    syms V1 [M M]
    rho_r = rho_e * eta;
    R1 = (rho_e * l_0)/(w_e*t_e);
    RL = (rho_r * t_h)/(w_e^2);
    Rc = (rho_e * (w_e + w_e))/(w_e*t_e);
    
    
    %% Setup System of Equations
    sys_eq = [];
    
    for i = 1:M 
        for j = 1:M
            [p_rs,p_cs] = find_neighbours(i,j,M);
            
            % Equation 0
            eq_0 = 0;
            eq_0 = eq_0 + (V1(i,j) - V0(i,j))/RL;
            for i_n = 1:size(p_rs,1)
                eq_0 = eq_0 + (V1(i,j) - V1(p_rs(i_n,1), p_rs(i_n,2)))/Rc;
            end
            if i == 1
                eq_0 = eq_0 + (V1(i,j) - x(j))/R1;
            end
    
            % Equation 1
            eq_1 = 0;
            eq_1 = eq_1 + (V0(i,j) - V1(i,j))/RL;
            for i_n = 1:size(p_cs,1)
                eq_1 = eq_1 + (V0(i,j) - V0(p_cs(i_n,1), p_cs(i_n,2)))/Rc;
            end
            if j == 1
                eq_1 = eq_1 + (V0(i,j) - y(i))/R1;
            end
    
            sys_eq = [sys_eq; eq_0 == 0; eq_1 == 0];
        end
    end
    
    %% Substitutions and solve system of equations
    
    Vx = [reshape(V0, M*M,1); reshape(V1, M*M,1)];
    [A,B] = equationsToMatrix(sys_eq, Vx);
    A = subs(A, [w_e t_e t_h rho_e l_0], [3.5e-3 0.2e-3 0.16e-3 1e-5 5e-3]);
    xin = [0.5*ones(1,M-1),1]; yin = [0.5*ones(1,M-1),0];
    B = subs(B, x, xin); B = subs(B, y, yin);
    B = subs(B, [w_e t_e t_h rho_e l_0], [3.5e-3 0.2e-3 0.16e-3 1e-5 5e-3]);
    V_res = inv(A)*B;
    V1_res = reshape(V_res(1:M*M), M,M);
    V0_res = reshape(V_res(M*M+1:end), M,M);
    V_dist = V1_res - V0_res;
end

function [ngbr_i, ngbr_j] = find_neighbours(i,j,M)
    ngbr_i = []; ngbr_j = [];

    if i ~= 1 & i ~= M
        ngbr_i = [ngbr_i; i+1 j; i-1 j;];
    elseif i == 1
        ngbr_i = [ngbr_i; i+1 j;];
    elseif i == M
        ngbr_i = [ngbr_i; i-1 j;];
    end

    if j ~= 1 & j ~= M
        ngbr_j = [ngbr_j; i j+1; i j-1];
    elseif j == 1
        ngbr_j = [ngbr_j; i j+1;];
    elseif j == M
        ngbr_j = [ngbr_j; i j-1;];
    end
end