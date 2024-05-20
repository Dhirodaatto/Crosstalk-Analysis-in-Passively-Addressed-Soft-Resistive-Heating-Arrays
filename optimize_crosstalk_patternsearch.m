tic;
%% Define lb,ub for objective function
M = 4;
e_w = 3.5e-3; Actuation_L = 20*e_w; ewa = e_w;
lb = [0.1; 0.1]; ub = [1; (Actuation_L - M*ewa)/((M-1)*ewa)];

plist = [M;ewa];

x0 = [1,1];

%% Define Objective Function

objective_function = @(x) objective_function_wrapper(plist, x);

%% Set Optimizer Options

setoutfun();
opts = optimoptions("patternsearch", ...
                    "MaxFunctionEvaluations", 60, ...
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
plot(outdata(:,1), outdata(:,4), "bo-");
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
function [J] = objective_function_wrapper(plist, x)
     param_list = [
        plist(1); % Number of Elements in Square Array
        plist(2); % Size of Pixel
        x(1); % Minimum Width factor
        x(2); % Spacing Factor
        0.2e-3; % Electrode Thickness
        0.16e-3; % High Resistive Layer Thickness
        10; % Convective Heat Flux Coefficient
        4.0322e4; % Base Electrode Conductivity
        5.32; % log10(Conductivity Ratio)
        70; % Thermal Conductivity of Electrode
    ];
    
    [ct, V] = evaluate_crosstalk(param_list, 0);

    if ct > 1
        J = 1;
    else
        J = ct;
    end
end

%% Define Output function
function [stop,options,optchanged] = outfun(optimval, options, flag)
stop = false;
optchanged = false;
global outdata;
outdata = [outdata; optimval.iteration optimval.x.' optimval.fval];

fprintf("Iteration # %d check \n", optimval.iteration);
end

function [] = setoutfun()
global outdata;
outdata = [];
end