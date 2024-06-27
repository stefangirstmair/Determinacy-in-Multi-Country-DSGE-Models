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
gamma_para = {'gamma_0_7'}; % Gamma parameter values

% Generate parameter value vectors for iteration
phi_pi_U_vec = linspace(0, 1.25, 50);
phi_pi_G_vec = linspace(0, 1.25, 50);
phi_pi_R_vec = linspace(0, 1.25, 50);
phi_Y_G_vec = [0 0.5000 1.0000];
phi_Y_U_vec = [0 0.5000 1.0000];
phi_Y_R_vec = [0 0.5000 1.0000];

% Create grids of parameter values
[phi_pi_U_mat, phi_pi_G_mat, phi_pi_R_mat] = meshgrid(phi_pi_U_vec, phi_pi_G_vec, phi_pi_R_vec);

% Define the total number of iterations
num_iterations_1 = length(phi_Y_G_vec);
num_iterations_2 = length(phi_Y_U_vec);
num_iterations_3 = length(phi_Y_R_vec);
num_iterations_4 = length(phi_pi_G_vec);
num_iterations_5 = length(phi_pi_U_vec);
num_iterations_6 = length(phi_pi_R_vec);
total_iterations = num_iterations_1 * num_iterations_2 * num_iterations_3 * num_iterations_4 * num_iterations_5 * num_iterations_6;

% Define the update interval for waitbar
update_interval = 1000; % Update the waitbar every 1000 iterations

% Create a waitbar object to display progress
h = waitbar(0, 'Progress: 0%');

% Initialize variables
info_mat = NaN(size(phi_pi_U_mat));
options_.qz_criterium = 1 + 1e-6;

% Loop through each paradigm
for p = 1:length(title_para)
    paradigm_currrent = title_para{p}; % Current paradigm name

    % Run Dynare according to the current paradigm and reset params
    run_folder = sprintf('%s/%s Model/%s/', home, paradigm_currrent, gamma_para{1});
    PathForPlot_current =  sprintf('%s/Plots/%s/%s', home, paradigm_currrent, gamma_para{1});
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
    data_structure = struct('phi_Y_G_iter', [], 'phi_Y_U_iter', [], 'phi_Y_R_iter', [], ...
        'phi_pi_G_iter', [], 'phi_pi_U_iter', [], 'phi_pi_R_iter', []);
    data_structure_indet = struct('phi_Y_G_iter', [], 'phi_Y_U_iter', [], 'phi_Y_R_iter', [], ...
        'phi_pi_G_iter', [], 'phi_pi_U_iter', [], 'phi_pi_R_iter', []);

    % Nested loops to iterate through parameter combinations
    for phi_Y_G_iter = 1:length(phi_Y_G_vec)
        set_param_value('phi_Y_G', phi_Y_G_vec(phi_Y_G_iter));
        for phi_Y_U_iter = 1:length(phi_Y_U_vec)
            set_param_value('phi_Y_U', phi_Y_U_vec(phi_Y_U_iter));
            for phi_Y_R_iter = 1:length(phi_Y_R_vec)
                set_param_value('phi_Y_R', phi_Y_R_vec(phi_Y_R_iter));
                for phi_pi_U_iter = 1:length(phi_pi_U_vec)
                    for phi_pi_G_iter = 1:length(phi_pi_G_vec)
                        for phi_pi_R_iter = 1:length(phi_pi_R_vec)
                            set_param_value('phi_pi_U', phi_pi_U_mat(phi_pi_U_iter, phi_pi_G_iter, phi_pi_R_iter));
                            set_param_value('phi_pi_G', phi_pi_G_mat(phi_pi_U_iter, phi_pi_G_iter, phi_pi_R_iter));
                            set_param_value('phi_pi_R', phi_pi_R_mat(phi_pi_U_iter, phi_pi_G_iter, phi_pi_R_iter));

                            % Resolve model and store results based on information criterion
                            [dr, info] = resol(0, M_, options_, oo_);
                            info_mat(phi_pi_U_iter, phi_pi_G_iter, phi_pi_R_iter) = info(1);

                            % Store data in appropriate structure based on information criterion value
                            % Store data in appropriate structure based on information criterion
                            if info_mat(phi_pi_U_iter, phi_pi_G_iter, phi_pi_R_iter) == 0
                                data_structure(iter_det) = struct('phi_Y_G_iter', phi_Y_G_vec(phi_Y_G_iter), ...
                                    'phi_Y_U_iter', phi_Y_U_vec(phi_Y_U_iter), ...
                                    'phi_Y_R_iter', phi_Y_R_vec(phi_Y_R_iter), ...
                                    'phi_pi_G_iter', phi_pi_G_vec(phi_pi_U_iter), ...
                                    'phi_pi_U_iter', phi_pi_U_vec(phi_pi_G_iter), ...
                                    'phi_pi_R_iter', phi_pi_R_vec(phi_pi_R_iter));
                                iter_det = iter_det + 1;
                            else
                                data_structure_indet(iter_indet) = struct('phi_Y_G_iter', phi_Y_G_vec(phi_Y_G_iter), ...
                                    'phi_Y_U_iter', phi_Y_U_vec(phi_Y_U_iter), ...
                                    'phi_Y_R_iter', phi_Y_R_vec(phi_Y_R_iter), ...
                                    'phi_pi_G_iter', phi_pi_G_vec(phi_pi_U_iter), ...
                                    'phi_pi_U_iter', phi_pi_U_vec(phi_pi_G_iter), ...
                                    'phi_pi_R_iter', phi_pi_R_vec(phi_pi_R_iter));
                                iter_indet = iter_indet + 1;
                            end

                            % Check if it's time to update the waitbar and display progress
                            current_iteration = (phi_Y_G_iter - 1) * num_iterations_2 * num_iterations_3 * num_iterations_4 * num_iterations_5 * num_iterations_6 + ...
                                (phi_Y_U_iter - 1) * num_iterations_3 * num_iterations_4 * num_iterations_5 * num_iterations_6 + ...
                                (phi_Y_R_iter - 1) * num_iterations_4 * num_iterations_5 * num_iterations_6 + ...
                                (phi_pi_U_iter - 1) * num_iterations_5 * num_iterations_6 + ...
                                (phi_pi_G_iter - 1) * num_iterations_6 + ...
                                phi_pi_R_iter;
                            if mod(current_iteration, update_interval) == 0 || current_iteration == total_iterations
                                waitbar(current_iteration / total_iterations, h, sprintf('%s Progress: %d%%', paradigm_currrent, round(current_iteration / total_iterations * 100)));
                            end
                        end
                    end
                end

                % Plotting based on information criterion values
                Z_plot = zeros(size(info_mat));
                Z_plot(info_mat == 0) = 1;

                % Create 3D scatter plot with alpha shapes
                figure()
                x = [phi_pi_U_mat(:); phi_pi_U_mat(:); phi_pi_U_mat(:)];
                y = [phi_pi_G_mat(:); phi_pi_G_mat(:); phi_pi_G_mat(:)];
                z = [phi_pi_R_mat(:); phi_pi_R_mat(:); phi_pi_R_mat(:)];
                determ = [Z_plot(:); Z_plot(:); Z_plot(:)];

                % Filter data based on determinacy
                fulfilling_points = determ == 1;
                x_fulfilling = x(fulfilling_points);
                y_fulfilling = y(fulfilling_points);
                z_fulfilling = z(fulfilling_points);

                nfulfilling_points = determ == 0;
                x_nfulfilling = x(nfulfilling_points);
                y_nfulfilling = y(nfulfilling_points);
                z_nfulfilling = z(nfulfilling_points);

                % alpha_value = 0.0221;

                % Create alpha shapes for fulfilling and non-fulfilling points
                alpha_shape = alphaShape(x_fulfilling, y_fulfilling, z_fulfilling);
                alpha_value = alpha_shape.Alpha;
                alpha_shape.Alpha = 0.0221;%alpha_value;
                alpha_shape2 = alphaShape(x_nfulfilling, y_nfulfilling, z_nfulfilling);
                alpha_value2 = alpha_shape2.Alpha;
                alpha_shape2.Alpha = 0.0221;%alpha_value2;

                % Plot alpha shapes
                plot(alpha_shape, 'FaceColor', 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
                hold on
                plot(alpha_shape2, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

                % Set labels and title for the plot
                xlabel('$\phi_{\pi}^U$', 'interpreter', 'latex', 'FontSize', 18)
                ylabel('$\phi_{\pi}^G$', 'interpreter', 'latex', 'FontSize', 18)
                zlabel('$\phi_{\pi}^R$', 'interpreter', 'latex', 'FontSize', 18)
                title(['Determinacy in ', title_para{p}, ': $\gamma = $', num2str(gamma), newline, ...
                    ' $\phi_{y}^G = $', num2str(phi_Y_G), ', $\phi_{y}^U = $', num2str(phi_Y_U), ...
                    ', $\phi_{y}^R = $', num2str(phi_Y_R)], 'interpreter', 'latex', 'FontSize', 18)

                % Save the figure with a specific filename and path
                savefig(gcf, fullfile(PathForPlot_current, strrep(sprintf('%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s', title_para{p}, ...
                    num2str(gamma), num2str(phi_Y_G), num2str(phi_Y_U), num2str(phi_Y_R)), '.', '_')))
                % Close the figure
                close all
            end
        end
    end
    % Construct filenames with desired format and path for data structures
    filename = sprintf('data_structure_%s_Det_gamma_%s_plot', title_para{p}, num2str(gamma));
    filepath = fullfile(save_runs_folder, strrep(filename, '.', '_'));

    % Save data structures to specified files
    save(filepath, 'data_structure');

    filename_Indet = sprintf('data_structure_%s_Indet_gamma_%s_plot', title_para{p}, num2str(gamma));
    filepath_Indet = fullfile(save_runs_folder, strrep(filename_Indet, '.', '_'));

    % Save data structures to specified files
    save(filepath_Indet, 'data_structure_indet');
end

% Close the waitbar after loops finish
close(h);

% Change back to home directory
cd(home);
