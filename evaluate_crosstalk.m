function [l_T_crosstalk, V_in] = evaluate_crosstalk(param_list, restartflag)
%EVALUATE_CROSSTALK Summary of this function goes here
%   Detailed explanation goes here


    %% Setups
    if restartflag
        homedir = pwd;
        cd('C:\D\COMSOL60\Multiphysics\bin\win64');
        system('comsolmphserver.exe &');
        cd('C:\D\COMSOL60\Multiphysics\mli');
        mphstart(2036);
        cd(homedir)
        import com.comsol.model.*;
        import com.comsol.model.util.*;
        ModelUtil.tags;
        ModelUtil.remove('Model');
    end
    
    %% Imports and Model Variables
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    ModelUtil.showProgress(false);
    model = ModelUtil.create('Model');
    model.modelPath(pwd)
    model.hist.disable;
    
    %% Create Component and Set Relevant Physics
    comp1 = model.component.create('comp1', true);
    geom1 = comp1.geom.create('geom1', 3);
    
    ht = comp1.physics.create('ht', 'HeatTransfer', 'geom1');
    ec = comp1.physics.create('ec', 'ConductiveMedia', 'geom1');
    emh = comp1.multiphysics.create('emh1', 'ElectromagneticHeating', 3);
    
    ss_study = model.study.create('std1');
    ss_study.create('stat', 'Stationary');
    ss_study.feature('stat').setSolveFor('/physics/ht', true);
    ss_study.feature('stat').setSolveFor('/physics/ec', true);
    
    %% Parameter List
    
    M = param_list(1); N = M; % MxN grid of pixels
    
    Mi = M; Ni = M; % Activated Pixels
    
    e_w = param_list(2); % Maximum Width of electrode
    e_wm = param_list(3) * e_w; % Minimum Width of electrode
    spacing = param_list(4) * param_list(2); % Spacing between pixels
    offset = 5e-3; % Starting offset between first pixel and start of electrode
    e_t = param_list(5); % Electrode Thickness
    hr_t = param_list(6); % Thickness of High Resistive Layer
    
    h = param_list(7); % Convective Heat Flux Coefficient
    m = 0.85; % Radiative Emissivity Constant
    T_inf = 20 + 273.15; % degC to K, Room Temperature
    
    % V_in_upper = sqrt(h*e_w^2*70*(1/hr_cond * hr_t/e_w^2))
    %            = sqrt(h*70*hr_t/hr_cond)
    V_in = 10; % Voltage applied to heat pixel
    
    % --- e_cond = 4.032258e4; hr_cond = 1.929e-1; % Experimental Data
    e_cond = param_list(8); % Electrode Electrical Conductivity
    conductivity_ratio = param_list(9); % Order of Magnitude difference between conductivities
    hr_cond = e_cond * 10^(-conductivity_ratio); % High Resistive Layer Electrical Conductivity
    
    e_k = param_list(10); % Electrode Thermal Conductivity
    hr_k = 0.25; % High Resistive Layer Thermal Conductivity
    
    T_req = 70+273.15; % Required Temperature
    mesh_size = 5; % Normal Mesh
    
    %% Geometry Creation and Setup
    
    % Aligned with X direction Bottom
    hor_array00 = geom1.feature.create('blke1', 'Block');
    hor_array00.set('size', [offset e_w e_t]);
    hor_array00.set('pos', [0 offset 0]);
    hor_array01 = geom1.feature.create('blke2', 'Block');
    hor_array01.set('size', [e_w e_w e_t]);
    hor_array01.set('pos', [offset offset 0]);
    hor_array02 = geom1.feature.create('blke3', 'Block');
    hor_array02.set('size', [spacing e_wm e_t]);
    hor_array02.set('pos', [offset+e_w offset+(e_w-e_wm)/2 0]);
    hor_array0 = geom1.feature.create('arr00', 'Array');
    hor_array0.selection('input').set('blke1');
    hor_array0.set('fullsize', [2 1 1]);
    hor_array0.set('displ', [offset+e_w*M+spacing*(M-1) 0 0]);
    hor_array1 = geom1.feature.create('arr01', 'Array');
    hor_array1.selection('input').set('blke2');
    hor_array1.set('fullsize', [M 1 1]);
    hor_array1.set('displ', [e_w+spacing 0 0]);
    hor_array2 = geom1.feature.create('arr02', 'Array');
    hor_array2.selection('input').set('blke3');
    hor_array2.set('fullsize', [M-1 1 1]);
    hor_array2.set('displ', [e_w+spacing 0 0]);
    hor_array = geom1.feature.create('unih', 'Union');
    hor_array.set('intbnd', false);
    hor_array.selection('input').set({'arr00' 'arr01' 'arr02'});
    hor_set = geom1.feature.create('arrh', 'Array');
    hor_set.selection('input').set('unih');
    hor_set.set('fullsize', [1 N 1]);
    hor_set.set('displ', [0 e_w+spacing 0]);
    
    % Aligned with Y direction Top
    vert_array00 = geom1.feature.create('blke4', 'Block');
    vert_array00.set('size', [e_w offset e_t]);
    vert_array00.set('pos', [offset 0 e_t+hr_t]);
    vert_array01 = geom1.feature.create('blke5', 'Block');
    vert_array01.set('size', [e_w e_w e_t]);
    vert_array01.set('pos', [offset offset e_t+hr_t]);
    vert_array02 = geom1.feature.create('blke6', 'Block');
    vert_array02.set('size', [e_wm spacing e_t]);
    vert_array02.set('pos', [offset+(e_w-e_wm)/2 offset+e_w e_t+hr_t]);
    vert_array0 = geom1.feature.create('arr10', 'Array');
    vert_array0.selection('input').set('blke4');
    vert_array0.set('fullsize', [1 2 1]);
    vert_array0.set('displ', [0 offset+e_w*N+spacing*(N-1) 0]);
    vert_array1 = geom1.feature.create('arr11', 'Array');
    vert_array1.selection('input').set('blke5');
    vert_array1.set('fullsize', [1 N 1]);
    vert_array1.set('displ', [0 e_w+spacing 0]);
    vert_array2 = geom1.feature.create('arr12', 'Array');
    vert_array2.selection('input').set('blke6');
    vert_array2.set('fullsize', [1 N-1 1]);
    vert_array2.set('displ', [0 e_w+spacing 0]);
    vert_array = geom1.feature.create('univ', 'Union');
    vert_array.set('intbnd', false);
    vert_array.selection('input').set({'arr10' 'arr11' 'arr12'});
    vert_set = geom1.feature.create('arrv', 'Array');
    vert_set.selection('input').set('univ');
    vert_set.set('fullsize', [M 1 1]);
    vert_set.set('displ', [e_w+spacing 0 0]);
    
    % High Resistive Squares in the middle layer
    hr_sq0 = geom1.feature.create('blkhr1', 'Block');
    hr_sq0.set('size', [e_w e_w hr_t]);
    hr_sq0.set('pos', [offset, offset, e_t]);
    hr_set = geom1.feature.create('arrhr', 'Array');
    hr_set.selection('input').set('blkhr1');
    hr_set.set('fullsize', [M N 1]);
    hr_set.set('displ', [e_w+spacing e_w+spacing 0]);
    
    % Finally Build Geometry
    geom1.run;
    
    %% Material Creation
    e_mat = comp1.material.create('e_mat');
    e_mat.label('Electrode Material');
    e_def = e_mat.materialmodel('def');
    e_def.set('electricconductivity', e_cond);
    e_def.set('thermalconductivity', e_k);
    e_def.set('density', 10.5e3);
    e_def.set('relpermittivity', 12.96);
    e_def.set('heatcapacity', 0.236);
    
    hr_mat = comp1.material.create('hr_mat');
    hr_mat.label('High Resistive Material');
    hr_def = hr_mat.materialmodel('def');
    hr_def.set('electricconductivity', hr_cond);
    hr_def.set('thermalconductivity', hr_k);
    hr_def.set('density', 2100);
    hr_def.set('relpermittivity', 2.5);
    hr_def.set('heatcapacity', 523.35);
    
    %% Helper Definitions
    
    tol_r = 1e-6;
    
    %% Assigning Geometry to Materials
    
    % Assign everything as electrode material
    e_mat.selection.all;
    
    % Create set of coords for hr layer
    clist = [];
    for i = 1:M
    for j = 1:N
        c0 = [offset+(i-1)*(e_w+spacing)-tol_r; offset+(j-1)*(e_w+spacing)-tol_r; e_t-tol_r];
        c1 = [offset+(i-1)*(e_w+spacing)+e_w+tol_r; offset+(j-1)*(e_w+spacing)+e_w+tol_r; e_t+hr_t+tol_r];
        idx = mphselectbox(model, 'geom1', [c0,c1], 'domain');
        clist(end+1) = idx;
    end
    end
    hr_mat.selection.set(clist);
    
    %% Setting up physics and BCs
    
    % Setting up Convective Heat Flux from all boundaries (should be changed) 
    ht.prop('PhysicalModelProperty').set('Tref', T_inf);
    ht.feature('init1').set('Tinit', T_inf);
    ht_hf = ht.create('hf1', 'HeatFluxBoundary', 2);
    ht_hf.selection.all;
    ht_hf.set('HeatFluxType', 'ConvectiveHeatFlux');
    ht_hf.set('h', h);
    ht_hf.set('Text', T_inf); % Text -> T_ext
    
    % Setting up Voltage Inputs
    gnd = ec.create('gnd', 'Ground', 2);
    c0 = [-tol_r; offset+(Ni-1)*(e_w+spacing)-tol_r; -tol_r];
    c1 = [+tol_r; offset+(Ni-1)*(e_w+spacing)+e_w+tol_r; e_t+tol_r];
    gnd_idx = mphselectbox(model, 'geom1', [c0,c1], 'boundary');
    gnd.selection.set(gnd_idx);
    vin = ec.create('pot', 'ElectricPotential', 2);
    c0 = [offset+(Mi-1)*(e_w+spacing)-tol_r; -tol_r; e_t+hr_t-tol_r];
    c1 = [offset+(Mi-1)*(e_w+spacing)+e_w+tol_r; +tol_r; e_t*2+hr_t+tol_r];
    vin_idx = mphselectbox(model, 'geom1', [c0,c1], 'boundary');
    vin.selection.set(vin_idx);
    vin.set('V0', V_in);
    
    %% Setting up the mesh
    
    mesh1 = comp1.mesh.create('mesh1', 'geom1');
    mesh1.automatic(true);
    mesh1.autoMeshSize(mesh_size);
    mesh1.run;
    
    %save(model, string(pwd+"\x.mph"));
    
    V_seed = [10 50];
    coeff_v = [];
    
    for V_in = V_seed
    vin.set('V0', V_in);
    ss_study.run;
    %% Helper Definitions
    
    tol_r = 1e-6;
    
    %% Compute Crosstalk
    
    T_dat = mpheval(model, 'T');
    T_s = T_dat.d1;
    x_s = T_dat.p(1,:);
    y_s = T_dat.p(2,:);
    z_s = T_dat.p(3,:);
    top_layer_indices = find(z_s >= (e_t*2 + hr_t - tol_r));
    x_top = x_s(top_layer_indices);
    y_top = y_s(top_layer_indices);
    T_top = T_s(top_layer_indices);
    pixel_avg_T = zeros(M,N);
    for i = 1:M
        for j = 1:N
            xb = [offset + e_w * (i-1) + spacing * (i-1), offset + e_w * i + spacing * (i-1)];
            yb = [offset + e_w * (j-1) + spacing * (j-1), offset + e_w * j + spacing * (j-1)];
            pixel_indices = intersect(top_layer_indices, intersect( find(x_s >= xb(1) & x_s <= xb(2)), find(y_s >= yb(1) & y_s <= yb(2))));
            pixel_Tset = T_s(pixel_indices);
    
            pixel_avg_T(i,j) = max(pixel_Tset);
        end
    end
    
    activated_pixel_temp = pixel_avg_T(Mi, Ni);
    is_pixel_active = activated_pixel_temp == max(pixel_avg_T, [], "all");
    
    line_set_T = [];
    
    for i = 1:M
        for j = 1:N
            % Creating line set
            if (j == Ni | i == Mi) & ~(j == Ni & i == Mi)
                line_set_T(end+1) = pixel_avg_T(i,j);
            end
        end
    end
    
    l_T_crosstalk = (max(line_set_T) - T_inf) ./ (pixel_avg_T(Mi,Ni) - T_inf);
    coeff_v(end+1) = (activated_pixel_temp - T_inf)/V_in^2;
    end
    
    coeff = mean(coeff_v);
    V_in = sqrt((T_req - T_inf)/coeff);
    vin.set('V0', V_in);
    ss_study.run;
    %% Helper Definitions
    
    tol_r = 1e-6;
    
    %% Compute Crosstalk
    
    T_dat = mpheval(model, 'T');
    T_s = T_dat.d1;
    x_s = T_dat.p(1,:);
    y_s = T_dat.p(2,:);
    z_s = T_dat.p(3,:);
    top_layer_indices = find(z_s >= (e_t*2 + hr_t - tol_r));
    x_top = x_s(top_layer_indices);
    y_top = y_s(top_layer_indices);
    T_top = T_s(top_layer_indices);
    pixel_avg_T = zeros(M,N);
    for i = 1:M
        for j = 1:N
            xb = [offset + e_w * (i-1) + spacing * (i-1), offset + e_w * i + spacing * (i-1)];
            yb = [offset + e_w * (j-1) + spacing * (j-1), offset + e_w * j + spacing * (j-1)];
            pixel_indices = intersect(top_layer_indices, intersect( find(x_s >= xb(1) & x_s <= xb(2)), find(y_s >= yb(1) & y_s <= yb(2))));
            pixel_Tset = T_s(pixel_indices);
    
            pixel_avg_T(i,j) = max(pixel_Tset);
        end
    end
    
    activated_pixel_temp = pixel_avg_T(Mi, Ni);
    is_pixel_active = activated_pixel_temp == max(pixel_avg_T, [], "all");
    
    line_set_T = [];
    
    for i = 1:M
        for j = 1:N
            % Creating line set
            if (j == Ni | i == Mi) & ~(j == Ni & i == Mi)
                line_set_T(end+1) = pixel_avg_T(i,j);
            end
        end
    end
    
    l_T_crosstalk = (max(line_set_T) - T_inf) ./ (pixel_avg_T(Mi,Ni) - T_inf);

    save(model, string(pwd)+ "\test.mph");

    %% Cleanup

    clear model;

    if restartflag
        system('Taskkill/IM WindowsTerminal.exe /f')
    end
end
