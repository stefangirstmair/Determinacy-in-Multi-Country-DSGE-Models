%% phiY analysis 1

phiy_ana   = {'U','G','R'};
title_para = {'DCP','PCP','LCP'}; %

home = pwd;

gamma_now = {'0_7','0_5'};

PathForPlot_Combined =  sprintf('%s/Plots/Plots_Combined',home);

phiYG_all = {'0','0_5','1'};
phiYU_all = {'0','0_5','1'};
phiYR_all = {'0','0_5','1'};

for p = 1:length(title_para)

    paradigm_currrent = title_para(p);

    % Set nec paths
    switch paradigm_currrent{1}
        case 'DCP'
            PathOfFiles_gamma_05 = sprintf('%s/Plots/DCP/gamma_0_5',home);
            PathOfFiles_gamma_07 = sprintf('%s/Plots/DCP/gamma_0_7',home);
        case 'PCP'
            PathOfFiles_gamma_05 = sprintf('%s/Plots/PCP/gamma_0_5',home);
            PathOfFiles_gamma_07 = sprintf('%s/Plots/PCP/gamma_0_7',home);
        case'LCP'
            PathOfFiles_gamma_05 = sprintf('%s/Plots/LCP/gamma_0_5',home);
            PathOfFiles_gamma_07 = sprintf('%s/Plots/LCP/gamma_0_7',home);
    end

    for k = 1:length(phiy_ana)

        phi_Y_currrent = phiy_ana(k);


        % Set nec paths
        switch phi_Y_currrent{1}
            case 'G'
                % Define which figures to include in plot
                plot_paths = {sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{1},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{2},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{3},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{1},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{2},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{3},phiYU_all{1},phiYR_all{1})};
            case 'U'
                % Define which figures to include in plot
                plot_paths = {sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{1},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{1},phiYU_all{2},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{1},phiYU_all{3},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{1},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{1},phiYU_all{2},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{1},phiYU_all{3},phiYR_all{1})};
            case'R'
                % Define which figures to include in plot
                plot_paths = {sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{1},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{1},phiYU_all{1},phiYR_all{2}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_07,paradigm_currrent{1},gamma_now{1},phiYG_all{1},phiYU_all{1},phiYR_all{3}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{1},phiYU_all{1},phiYR_all{1}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{1},phiYU_all{1},phiYR_all{2}),...
                    sprintf('%s/%s_gamma_%s_phiyG_%s_phiyU_%s_phiyR_%s.fig',PathOfFiles_gamma_05,paradigm_currrent{1},gamma_now{2},phiYG_all{1},phiYU_all{1},phiYR_all{3})};
        end


        % Open full screen figure
        figure('units','normalized','outerposition',[0 0 1 1])

        % Set tiled layout with according size (rows, columns)
        tc = tiledlayout(2,3);

        % Adjust TileSpacing and Padding properties
        tc.TileSpacing = 'compact';  % Set TileSpacing to 'compact' for tighter spacing between tiles
        tc.Padding = 'compact';      % Set Padding to 'compact' for tighter padding around the layout

        % Preallocate for plotting (does not have to be adjusted)
        N=min(numel(plot_paths), prod(tc.GridSize));
        ax=gobjects(N,1);

        % Plot
        for i=1:N
            ax(i) = loadTile(tc,plot_paths{i},i);
        end

        % Adjust the outer position of the figure
        figure_handle = gcf;
        % figure_handle.OuterPosition(4) = 0.6; % Adjust the height as needed

        % Save the figure as a PNG
        filename = sprintf('%s_phiy%s_analysis_1_07_05.png',paradigm_currrent{1},phi_Y_currrent{1});
        % filename = 'DCP_phiyR_ana.png';
        % saveas(figure_handle, filename);
        saveas(figure_handle, fullfile(PathForPlot_Combined,filename)); % strrep replaces dot from decimal with underscore

        % Close the figure
        close(figure_handle);
    end
end

%%

%% phiY analysis 2

clear
clc

home = pwd;

title_para = {'DCP','PCP','LCP'}; %
gamma_now = {'0_7','0_5'};

PathForPlot_Combined =  sprintf('%s/Plots/Plots_Combined',home);


for p = 1:length(title_para)

    paradigm_currrent = title_para(p);

    PathOfFiles = sprintf('%s/Plots/%s',home,paradigm_currrent{1});

    % Define which figures to include in plot
    plot_paths = {sprintf('%s/gamma_0_7/%s_gamma_0_7_phiyG_0_5_phiyU_0_phiyR_0_5.fig',PathOfFiles,paradigm_currrent{1}),...
        sprintf('%s/gamma_0_7/%s_gamma_0_7_phiyG_0_5_phiyU_0_5_phiyR_0_5.fig',PathOfFiles,paradigm_currrent{1}),...
        sprintf('%s/gamma_0_7/%s_gamma_0_7_phiyG_0_5_phiyU_1_phiyR_0_5.fig',PathOfFiles,paradigm_currrent{1}),...
        sprintf('%s/gamma_0_5/%s_gamma_0_5_phiyG_0_5_phiyU_0_phiyR_0_5.fig',PathOfFiles,paradigm_currrent{1}),...
        sprintf('%s/gamma_0_5/%s_gamma_0_5_phiyG_0_5_phiyU_0_5_phiyR_0_5.fig',PathOfFiles,paradigm_currrent{1}),...
        sprintf('%s/gamma_0_5/%s_gamma_0_5_phiyG_0_5_phiyU_1_phiyR_0_5',PathOfFiles,paradigm_currrent{1})};


    % Open full screen figure
    figure('units','normalized','outerposition',[0 0 1 1])

    % Set tiled layout with according size (rows, columns)
    tc = tiledlayout(2,3);

    % Adjust TileSpacing and Padding properties
    tc.TileSpacing = 'compact';  % Set TileSpacing to 'compact' for tighter spacing between tiles
    tc.Padding = 'compact';      % Set Padding to 'compact' for tighter padding around the layout

    % Preallocate for plotting (does not have to be adjusted)
    N=min(numel(plot_paths), prod(tc.GridSize));
    ax=gobjects(N,1);

    % Plot
    for i=1:N
        ax(i) = loadTile(tc,plot_paths{i},i);
    end

    % Adjust the outer position of the figure
    figure_handle = gcf;
    % figure_handle.OuterPosition(4) = 0.6; % Adjust the height as needed

    % Save the figure as a PNG
    filename = sprintf('%s_phiyU_analysis_2_07_05.png',paradigm_currrent{1});
    % filename = 'DCP_all_phiy_ana2.png';
    % saveas(figure_handle, filename);
    saveas(figure_handle, fullfile(PathForPlot_Combined,filename)); % strrep replaces dot from decimal with underscore

    % Close the figure
    close(figure_handle);


end


%% Plot Correlation

clear
clc

home = pwd;
Plots_Correlation =  sprintf('%s/Plots/Correlations', home);


title_para = {'DCP','PCP','LCP'};

make_corr_fig = 1;

if make_corr_fig == 1
    for p = 1:length(title_para)

        paradigm_current = title_para{p};

        for j = 1:2

            if j == 1
                gamma = 0.5;
            elseif j == 2
                gamma = 0.7;
            end

            % Load data
            % load(strrep(sprintf('%s/Save Runs/data_structure_%s_Det_gamma_%s_long',home,title_para{p},num2str(gamma)),'.','_'))
            load(strrep(sprintf('%s/Save Runs/data_structure_%s_Det_gamma_%s_plot',home,title_para{p},num2str(gamma)),'.','_'))
            % load(strrep(sprintf('%s/Save Runs/data_structure_%s_Indet_gamma_%s_long',home,title_para{p},num2str(gamma)),'.','_'))

            % Create data matrices
            data_matrix = [data_structure.phi_Y_G_iter; data_structure.phi_Y_U_iter; data_structure.phi_Y_R_iter; data_structure.phi_pi_G_iter; data_structure.phi_pi_U_iter; data_structure.phi_pi_R_iter]';
            % data_matrix_indet = [data_structure_indet.phi_Y_G_iter; data_structure_indet.phi_Y_U_iter; data_structure_indet.phi_Y_R_iter; data_structure_indet.phi_pi_G_iter; data_structure_indet.phi_pi_U_iter; data_structure_indet.phi_pi_R_iter]';

            parameter_names = {'\phi_{\pi}^G', '\phi_{\pi}^U', '\phi_{\pi}^R','\phi_{y}^G', '\phi_{y}^U', '\phi_{y}^R'};

            % Compute correlation matrices and p-values
            [correlation_matrix, p_values] = corrcoef(data_matrix);
            % [correlation_matrix_indet, p_values_indet] = corrcoef(data_matrix_indet);

            % Create custom colormap
            num_colors = 128/2; % Number of colors in the colormap
            blue_shades = winter(num_colors); % Shades of blue for negative correlations
            yellow_red_shades = flipud(autumn(num_colors)); % Shades of yellow to red for positive correlations
            light_gray_color = [0.8, 0.8, 0.8]; % Light gray color for zero correlation
            custom_colormap = [blue_shades; light_gray_color; yellow_red_shades]; % Combine blue, light gray, and yellow to red shades

            % Visualize the correlation matrix with the custom colormap
            figure('units','normalized','outerposition',[0 0 1 1])
            imagesc(correlation_matrix, [-1, 1]);
            colorbar;
            title(['Determinacy in ',sprintf('%s: ', title_para{p}), '$\gamma= $' num2str(gamma), newline, 'Correlation Matrix' ], 'interpreter','latex', 'FontSize',38);
            xticklabels(parameter_names);
            yticklabels(parameter_names);
            colormap(custom_colormap);
            set(gca, 'FontSize', 20);

            % Add text annotations with correlation coefficients and significance stars (lower triangular part)
            [rows, cols] = size(correlation_matrix);
            for r = 1:rows
                for c = 1:r % Iterate only over the lower triangular part and main diagonal
                    % Display correlation coefficient
                    correlation_value = correlation_matrix(r, c);
                    text_str = sprintf('%.3f', correlation_value);

                    % Determine number of stars based on p-value
                    if p_values(r, c) < 0.01
                        stars = '***'; % Three stars for p < 0.01
                    elseif p_values(r, c) < 0.05
                        stars = '**'; % Two stars for p < 0.05
                    elseif p_values(r, c) < 0.1
                        stars = '*'; % One star for p < 0.1
                    else
                        stars = ''; % No star if not significant
                    end

                    % Append stars to the correlation coefficient text
                    text_str = [text_str stars];

                    % Display the text
                    text(c, r, text_str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontSize', 12);
                end
            end

            % Save figure
            cd(Plots_Correlation)
            savefig(strrep(sprintf('Correlation_Det_%s_gamma_%s',title_para{p},num2str(gamma)),'.','_')) % strrep replaces dot from decimal with underscore
            cd(home)

            % Close all figures
            close all

        end
    end
end



% PathForPlot_Corr =  '/home/stefan/Dropbox/Stefan/UNI/PhD/DCP Determinacy/Code/Det_Analysis/Plots/Plots_Correlation';
%
% Define which figures to include in plot
plot_paths = {sprintf('%s/Correlation_Det_DCP_gamma_0_7.fig',Plots_Correlation),...
    sprintf('%s/Correlation_Det_PCP_gamma_0_7.fig',Plots_Correlation),...
    sprintf('%s/Correlation_Det_LCP_gamma_0_7.fig',Plots_Correlation),...
    sprintf('%s/Correlation_Det_DCP_gamma_0_5.fig',Plots_Correlation),...
    sprintf('%s/Correlation_Det_PCP_gamma_0_5.fig',Plots_Correlation),...'
    sprintf('%s/Correlation_Det_LCP_gamma_0_5.fig',Plots_Correlation)};

% Open full screen figure
figure('units','normalized','outerposition',[0 0 1 1])

% Set tiled layout with according size (rows, columns)
tc = tiledlayout(2,3);

% Apply custom colormap
colormap(custom_colormap);


% Adjust TileSpacing and Padding properties
tc.TileSpacing = 'compact';  % Set TileSpacing to 'compact' for tighter spacing between tiles
tc.Padding = 'compact';      % Set Padding to 'compact' for tighter padding around the layout

% Preallocate for plotting (does not have to be adjusted)
N=min(numel(plot_paths), prod(tc.GridSize));
ax=gobjects(N,1);

% Plot
for i=1:N
    ax(i) = loadTile(tc,plot_paths{i},i);
end

% Adjust the outer position of the figure
figure_handle = gcf;
% figure_handle.OuterPosition(4) = 0.6; % Adjust the height as needed

% Save the figure as a PNG
% filename = 'Correlation_Plot.png';
filename = 'Correlation_Plot_long_07_05.png';
cd Plots/Correlations
saveas(figure_handle, filename);
cd ..
cd ..
% saveas(figure_handle, fullfile(PathForPlot_Corr,filename));
% Close the figure
close(figure_handle);



%% Same Cali

clear
clc

home = pwd;


% for p = 1:length(title_para)

paradigm_currrent = {'DCP'};%title_para(p);

PathOfFiles = sprintf('%s/Plots/SameCali/%s',home,paradigm_currrent{1});

% Define which figures to include in plot
plot_paths = {sprintf('%s/gamma_0_5/%s_gamma_0_5.fig',PathOfFiles,paradigm_currrent{1}),...
    sprintf('%s/gamma_0_7/%s_gamma_0_7.fig',PathOfFiles,paradigm_currrent{1})};


% Open full screen figure
figure('units','normalized','outerposition',[0 0 1 1])

% Set tiled layout with according size (rows, columns)
tc = tiledlayout(1,2);

% Adjust TileSpacing and Padding properties
tc.TileSpacing = 'compact';  % Set TileSpacing to 'compact' for tighter spacing between tiles
tc.Padding = 'compact';      % Set Padding to 'compact' for tighter padding around the layout

% Preallocate for plotting (does not have to be adjusted)
N=min(numel(plot_paths), prod(tc.GridSize));
ax=gobjects(N,1);

% Plot
for i=1:N
    ax(i) = loadTile(tc,plot_paths{i},i);
end

map = [0.8 0.8 0.8; 0 0.8 0.15];
colormap(map);
h.FaceAlpha = 0.6;

% Adjust the outer position of the figure
figure_handle = gcf;
% figure_handle.OuterPosition(4) = 0.6; % Adjust the height as needed

% Save the figure as a PNG
filename = sprintf('Same_Cali_%s.png',paradigm_currrent{1});
% filename = 'DCP_all_phiy_ana2.png';
% saveas(figure_handle, filename);
saveas(figure_handle, fullfile(PathForPlot_Combined,filename)); % strrep replaces dot from decimal with underscore

% Close the figure
close(figure_handle);





PathOfFiles = sprintf('%s/Plots/SameCali/%s',home);

% Define which figures to include in plot
plot_paths = {sprintf('%s/PCP/gamma_0_5/PCP_gamma_0_5.fig',PathOfFiles),...
    sprintf('%s/PCP/gamma_0_7/PCP_gamma_0_7.fig',PathOfFiles),...
    sprintf('%s/LCP/gamma_0_5/LCP_gamma_0_5.fig',PathOfFiles),...
    sprintf('%s/LCP/gamma_0_7/LCP_gamma_0_7.fig',PathOfFiles)};


% Open full screen figure
figure('units','normalized','outerposition',[0 0 1 1])

% Set tiled layout with according size (rows, columns)
tc = tiledlayout(2,2);

% Adjust TileSpacing and Padding properties
tc.TileSpacing = 'compact';  % Set TileSpacing to 'compact' for tighter spacing between tiles
tc.Padding = 'compact';      % Set Padding to 'compact' for tighter padding around the layout

% Preallocate for plotting (does not have to be adjusted)
N=min(numel(plot_paths), prod(tc.GridSize));
ax=gobjects(N,1);

% Plot
for i=1:N
    ax(i) = loadTile(tc,plot_paths{i},i);
end

map = [0.8 0.8 0.8; 0 0.8 0.15];
colormap(map);
h.FaceAlpha = 0.6;

% Adjust the outer position of the figure
figure_handle = gcf;
% figure_handle.OuterPosition(4) = 0.6; % Adjust the height as needed

% Save the figure as a PNG
filename = 'Same_Cali_PCP_and_LCP.png';
% filename = 'DCP_all_phiy_ana2.png';
% saveas(figure_handle, filename);
saveas(figure_handle, fullfile(PathForPlot_Combined,filename)); % strrep replaces dot from decimal with underscore

% Close the figure
close(figure_handle);

%% Same Cali all in one plot

clear
clc

home = pwd;


PathOfFiles = sprintf('%s/Plots/SameCali/%s',home);

% Define which figures to include in plot
plot_paths = {sprintf('%s/DCP/gamma_0_5/DCP_gamma_0_5.fig',PathOfFiles),...
    sprintf('%s/DCP/gamma_0_7/DCP_gamma_0_7.fig',PathOfFiles),...
    sprintf('%s/PCP/gamma_0_5/PCP_gamma_0_5.fig',PathOfFiles),...
    sprintf('%s/PCP/gamma_0_7/PCP_gamma_0_7.fig',PathOfFiles),...
    sprintf('%s/LCP/gamma_0_5/LCP_gamma_0_5.fig',PathOfFiles),...
    sprintf('%s/LCP/gamma_0_7/LCP_gamma_0_7.fig',PathOfFiles)};


% Open full screen figure
figure('units','normalized','outerposition',[0 0 1 1])

% Set tiled layout with according size (rows, columns)
tc = tiledlayout(3,2);

% Adjust TileSpacing and Padding properties
tc.TileSpacing = 'compact';  % Set TileSpacing to 'compact' for tighter spacing between tiles
tc.Padding = 'compact';      % Set Padding to 'compact' for tighter padding around the layout

% Preallocate for plotting (does not have to be adjusted)
N=min(numel(plot_paths), prod(tc.GridSize));
ax=gobjects(N,1);

% Plot
for i=1:N
    ax(i) = loadTile(tc,plot_paths{i},i);
end

map = [0.8 0.8 0.8; 0 0.8 0.15];
colormap(map);
h.FaceAlpha = 0.6;

% Adjust the outer position of the figure
figure_handle = gcf;
% figure_handle.OuterPosition(4) = 0.6; % Adjust the height as needed

% Save the figure as a PNG
filename = 'Same_Cali_DCP_PCP_and_LCP.png';
% filename = 'DCP_all_phiy_ana2.png';
% saveas(figure_handle, filename);
saveas(figure_handle, fullfile(PathForPlot_Combined,filename)); % strrep replaces dot from decimal with underscore

% Close the figure
close(figure_handle);

%% Same Cali all in one plot v2

clear
clc

home = pwd;


PathOfFiles = sprintf('%s/Plots/SameCali/%s',home);

% Define which figures to include in plot
plot_paths = {sprintf('%s/DCP/gamma_0_5/DCP_gamma_0_5_diff_calibration.fig',PathOfFiles),...
    sprintf('%s/DCP/gamma_0_7/DCP_gamma_0_7_diff_calibration.fig',PathOfFiles),...
    sprintf('%s/PCP/gamma_0_5/PCP_gamma_0_5_diff_calibration.fig',PathOfFiles),...
    sprintf('%s/PCP/gamma_0_7/PCP_gamma_0_7_diff_calibration.fig',PathOfFiles),...
    sprintf('%s/LCP/gamma_0_5/LCP_gamma_0_5_diff_calibration.fig',PathOfFiles),...
    sprintf('%s/LCP/gamma_0_7/LCP_gamma_0_7_diff_calibration.fig',PathOfFiles)};


% Open full screen figure
figure('units','normalized','outerposition',[0 0 1 1])

% Set tiled layout with according size (rows, columns)
tc = tiledlayout(3,2);

% Adjust TileSpacing and Padding properties
tc.TileSpacing = 'compact';  % Set TileSpacing to 'compact' for tighter spacing between tiles
tc.Padding = 'compact';      % Set Padding to 'compact' for tighter padding around the layout

% Preallocate for plotting (does not have to be adjusted)
N=min(numel(plot_paths), prod(tc.GridSize));
ax=gobjects(N,1);

% Plot
for i=1:N
    ax(i) = loadTile(tc,plot_paths{i},i);
end

map = [0.8 0.8 0.8; 0 0.8 0.15];
colormap(map);
h.FaceAlpha = 0.6;

% Adjust the outer position of the figure
figure_handle = gcf;
% figure_handle.OuterPosition(4) = 0.6; % Adjust the height as needed

% Save the figure as a PNG
filename = 'Same_Cali_DCP_PCP_and_LCP_v2.png';
% filename = 'DCP_all_phiy_ana2.png';
% saveas(figure_handle, filename);
saveas(figure_handle, fullfile(PathForPlot_Combined,filename)); % strrep replaces dot from decimal with underscore

% Close the figure
close(figure_handle);

