% Clear command window and workspace
clear
clc

% Store current working directory (home directory path)
home = pwd;

% Create path to the 'Save Runs' folder within the home directory
save_runs_folder = fullfile(home, 'Save Runs');

% Check if the 'Save Runs' folder exists; if not, create it
if ~exist(save_runs_folder, 'dir')
    mkdir(save_runs_folder);
end

% Create path to the 'Plots' folder within the home directory
PathForPlot =  fullfile(home, 'Plots');

% Check if the 'Plots' folder exists; if not, create it
if ~exist(PathForPlot, 'dir')
    mkdir(PathForPlot);
end

% Define arrays of parameters and values for iteration
title_para = { 'DCP', 'PCP','LCP'}; % Paradigm names 
gamma_para = {'gamma_0_5'}; % Gamma parameter values: Change to: 'gamma_0_7'

% Generate parameter value vectors for iteration
phi_pi_vec=linspace(0,1.25,250);
phi_Y_vec=linspace(0,2,250);

% Create grids of parameter values
[phi_pi_mat,phi_Y_mat]=meshgrid(phi_pi_vec,phi_Y_vec);

% Define the total number of iterations
num_iterations_1 = length(phi_pi_vec);
num_iterations_2 = length(phi_Y_vec);
total_iterations = num_iterations_1 * num_iterations_2;

% Define the update interval for waitbar
update_interval = 1000; % Update the waitbar every 1000 iterations

% Initialize variables
info_mat = NaN(size(phi_pi_mat));
options_.qz_criterium = 1 + 1e-6;

% Loop through each paradigm
for p = 1:length(title_para)

    % Create a waitbar object to display progress
h = waitbar(0, 'Progress: 0%');
    paradigm_currrent = title_para{p}; % Current paradigm name

    % Run Dynare according to the current paradigm and reset params
    run_folder = sprintf('%s/%s_Same_Cali/%s/', home, paradigm_currrent, gamma_para{1});
    PathForPlot_current =  sprintf('%s/Plots/SameCali/%s/%s', home, paradigm_currrent, gamma_para{1});
    switch paradigm_currrent
        case 'DCP'
            cd(run_folder)
            dynare DCP_D
            options_.qz_criterium = 1 + 1e-6;
        case 'PCP'
            cd(run_folder)
            dynare PCP_D
            options_.qz_criterium = 1 + 1e-6;
        case 'LCP'
            cd(run_folder)
            dynare LCP_D
            options_.qz_criterium = 1 + 1e-6;
    end

    % Initialize iteration counters
    iter_det = 1;
    iter_indet = 1;

    % Initialize data structures for storing results
    data_structure = struct('phi_Y_iter', [], ...
        'phi_pi_iter', []);
    data_structure_indet = struct('phi_Y_iter', [], ...
        'phi_pi_iter', []);

    % Nested loops to iterate through parameter combinations
    for phi_Y_iter = 1:length(phi_Y_vec)
        % set_param_value('phi_Y', phi_Y_vec(phi_Y_iter));
        for phi_pi_iter = 1:length(phi_pi_vec)
            set_param_value('phi_pi',phi_pi_mat(phi_pi_iter,phi_Y_iter));
            set_param_value('phi_Y',phi_Y_mat(phi_pi_iter,phi_Y_iter));


            % Resolve model and store results based on information criterion
            [dr, info] = resol(0, M_, options_, oo_);
            info_mat(phi_pi_iter,phi_Y_iter)=info(1);

            % Store data in appropriate structure based on information criterion value
            % Store data in appropriate structure based on information criterion
            if info_mat(phi_pi_iter,phi_Y_iter) == 0
                data_structure(iter_det) = struct('phi_Y_iter', phi_Y_vec(phi_Y_iter), ...
                    'phi_pi_iter', phi_pi_vec(phi_pi_iter));
                iter_det = iter_det + 1;
            else
                data_structure_indet(iter_indet) = struct('phi_Y_iter', phi_Y_vec(phi_Y_iter), ...
                    'phi_pi_iter', phi_pi_vec(phi_pi_iter));
                iter_indet = iter_indet + 1;
            end

            % Check if it's time to update the waitbar and display progress
            current_iteration = (phi_Y_iter - 1) * num_iterations_2 + phi_pi_iter;
            if mod(current_iteration, update_interval) == 0 || current_iteration == total_iterations
                waitbar(current_iteration / total_iterations, h, sprintf('%s Progress: %d%%', paradigm_currrent, round(current_iteration / total_iterations * 100)));
            end
        end

    end
    % Plotting based on information criterion values
    Z_plot = zeros(size(info_mat));
    Z_plot(info_mat == 0) = 1;

    % Create 3D scatter plot with alpha shapes
    figure()
    [c,h] =  contourf(phi_pi_mat,phi_Y_mat,Z_plot,1);
    map = [0.8 0.8 0.8; 0 0.8 0.15];
    colormap(map);
    h.FaceAlpha = 0.6;
    grid on

    xlabel('$\phi_{\pi}$','interpreter','latex', 'FontSize', 38, 'FontWeight', 'bold')
    ylabel('$\phi_Y$','interpreter','latex', 'FontSize', 38, 'FontWeight', 'bold')
    title(['Determinacy in ',sprintf('%s: ',title_para{p}),'$\gamma= $' num2str(gamma)],'interpreter','latex', 'FontSize',38)
    ax = gca;
    ax.FontWeight = 'bold';
    % Save the figure with a specific filename and path
    savefig(gcf, fullfile(PathForPlot_current, strrep(sprintf('%s_gamma_%s', title_para{p}, ...
        num2str(gamma)), '.', '_')))
    % Close the figure
    close all

    % Construct filenames with desired format and path for data structures
    filename = sprintf('data_structure_%s_Det_gamma_%s_plot_Same_Cali', title_para{p}, num2str(gamma));
    filepath = fullfile(save_runs_folder, strrep(filename, '.', '_'));

    % Save data structures to specified files
    save(filepath, 'data_structure');

end



% Change back to home directory
cd(home);
