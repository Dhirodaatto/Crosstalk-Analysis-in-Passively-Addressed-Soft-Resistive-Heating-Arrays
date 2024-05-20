tic;
%% Define lb,ub for objective function
M = 8;
e_w = 5e-3; ewa = e_w;
L_tot = 15e-2; A_f = M*e_w/L_tot;   
lb = [0.1; 0.1; 1]; ub = [1; (L_tot - M*e_w)/((M-1)*e_w) ;10];
optimization_weights = [0.1; 0.9];

plist = [M;ewa];

x0 = [0.75;0.5;5];
                
%% Define Objective Function

objective_function = @(x) objective_function_wrapper(optimization_weights, plist, x);

%% Set Optimizer Options

setoutfun();
opts = optimoptions("patternsearch", ...
                    "MaxFunctionEvaluations", 300, ...
                    "OutputFcn",@outfun, ...
                    "AccelerateMesh",true, ...
                    "StepTolerance", 1e-2, ...
                    "MeshTolerance", 1e-2);
[xopt, fval, exitflag, op] = patternsearch(objective_function, x0, [], [], [], [], lb,ub, [], opts);

toc;

%% Show Results

global outdata;
writematrix(outdata, "tempopt.csv");
outdata = outdata(1:end-1, :);

hfig = figure(1); fname = "opt";
plot(outdata(:,1), outdata(:,5), "bo-");    
xlabel("\# of Iterations"); ylabel("Evaluated Objective Function");
title("Optimization Results");

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',15) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-vector');


%% Function Wrapper using Optimizing Variables
function [J] = objective_function_wrapper(w, plist, x)
    % fprintf("%.3f, %.3f", x(1), x(2));
    Vmax = 60;
    ctb = 0; jb = 0;
     param_list = [
        plist(1); % Number of Elements in Square Array
        plist(2); % Size of Pixel
        x(1); % Minimum Width factors
        x(2); % Spacing Factor
        0.2e-3; % Electrode Thickness
        0.16e-3; % High Resistive Layer Thickness
        10; % Convective Heat Flux Coefficient
        4.0322e4; % Base Electrode Conductivity
        x(3); % log10(Conductivity Ratio)
        70; % Thermal Conductivity of Electrode
    ];
    
    [ct, V] = evaluate_crosstalk(param_list, 0);

    if ct > 1
        ct = 1; ctb = 1;
    end

    if V > Vmax
        V = Vmax; jb = 1;
    end

    V = V/Vmax;

    J = w(1)*ct + w(2)*V;

    if ctb | jb
        J = 1;
    end

    fprintf("%.2f %.2f %.2f : %.2f %.2f %.2f\n", x, ct, V, J);
end

%% Define Output function
function [stop,options,optchanged] = outfun(optimval, options, flag)
stop = false;
optchanged = false;
global outdata;
outdata = [outdata; optimval.iteration optimval.x.' optimval.fval];

fprintf("Iteration # %d check:-\n -- Method = %s TolX = %.2f TolFun = %.2f MeshSize = %.2f\n", optimval.iteration, optimval.method, optimval.TolX, optimval.TolFun, optimval.meshsize);
fprintf(" ---- \n");
end

function [] = setoutfun()
global outdata;
outdata = [];
end