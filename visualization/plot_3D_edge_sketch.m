%> Automatically plot every 3D edge file in the folder with a different color
%
%> (c) LEMS, Brown University
%> chiang-heng chien

% Define the folder path containing the 3D edge output files
data_folder_path = "/gpfs/data/bkimia/zqiwu/3D_Edge_Sketch/outputs/";

% Specify the common pattern in the file names
file_pattern = "Gamma1s_*.txt";

% Get all files matching the pattern
edge_files = dir(fullfile(data_folder_path, file_pattern));

% Define a set of colors to be used for different files
colors = lines(length(edge_files));  % Generate a set of distinct colors using the 'lines' colormap

% Create a figure for plotting
figure;
hold on;

% Loop through each file and plot its edges in 3D
for i = 1:length(edge_files)
    % Read the current file
    current_file_path = fullfile(data_folder_path, edge_files(i).name);
    edges_file_read = fopen(current_file_path, 'r');
    ldata = textscan(edges_file_read, '%f\t%f\t%f', 'CollectOutput', true);
    edges_3d = double(ldata{1,1});
    fclose(edges_file_read);

    % Plot the edges using a different color for each file
    plot3(edges_3d(:,1), edges_3d(:,2), -edges_3d(:,3), 'Color', colors(i, :), 'Marker', '.', 'LineStyle', 'none');
end

% Set the plot settings
[az, el] = view(gca);
axis equal;
axis off;
set(gcf, 'color', 'w');

% Add a legend for each file
legend({edge_files.name}, 'Interpreter', 'none');
hold off;
