%> Automatically plot every 3D edge file in the folder with a different color
%
%> (c) LEMS, Brown University
%> chiang-heng chien
clear;
clc;
close all;

%> Define the folder path containing the 3D edge output files
data_folder_name = 'outputs';
data_folder_path = fullfile(fileparts(mfilename('fullpath')), '..', data_folder_name);

%> Specify the common pattern in the file names
edge_location_file_pattern = "3D_edges_*.txt";
edge_orientation_file_pattern = "3D_tangents_*.txt";

%> Get all files matching the pattern
edge_location_files = dir(fullfile(data_folder_path, edge_location_file_pattern));
edge_tangents_files = dir(fullfile(data_folder_path, edge_orientation_file_pattern));

%> Define a set of colors to be used for different files through 'lines' colormap
colors = lines(length(edge_location_files)); 

%> Create a figure for plotting
figure;
ax = axes();
hold on;

quiver_mag = 0.005;

%> Loop through each file and plot its edges in 3D
hLegs = [];
for i = 1:length(edge_location_files)
    %> Read the current edge location file
    current_location_file_path = fullfile(data_folder_path, edge_location_files(i).name);
    edge_location_file_read = fopen(current_location_file_path, 'r');
    disp(current_location_file_path);

    %> Read the current edge orientation file
    current_tangent_file_path = fullfile(data_folder_path, edge_tangents_files(i).name);
    edge_tangent_file_read = fopen(current_tangent_file_path, 'r');
    disp(current_tangent_file_path);

    %> parse 3D edge location data
    ldata = textscan(edge_location_file_read, '%f\t%f\t%f', 'CollectOutput', true);
    edges_location_3d = double(ldata{1,1});
    fclose(edge_location_file_read);

    %> parse 3D edge orientation data
    ldata = textscan(edge_tangent_file_read, '%f\t%f\t%f', 'CollectOutput', true);
    edges_tangent_3d = double(ldata{1,1});
    fclose(edge_tangent_file_read);

    %> Get the legend
    hypothesis_view1_index = extractBetween(edge_location_files(i).name, 'hypo1_', '_hypo2');
    hypothesis_view2_index = extractBetween(edge_location_files(i).name, 'hypo2_', '_t');
    show_legend = strcat("3D edges from hypothesis views ", hypothesis_view1_index, " and ", hypothesis_view2_index);

    %> (Optional) remove 3D edges that are too "isolated" in the 3D space
    % invalid_edge_index = [];
    % for ei = 1:size(edges_location_3d, 1)
    %     dist_norms = vecnorm(edges_location_3d(ei,:) - edges_location_3d, 2, 2);
    %     yy = find(dist_norms < 0.04);
    %     if length(yy) < 4
    %         invalid_edge_index = [invalid_edge_index; ei];
    %     end
    % end
    % edges_3d(invalid_edge_index, :) = [];

    %> Plot the edges using a different color for each file
    h = plot3(edges_location_3d(:,1), edges_location_3d(:,2), edges_location_3d(:,3), ...
             'Color', colors(i, :), 'Marker', '.', 'MarkerSize', 3.5, 'LineStyle', 'none', ...
             'DisplayName', show_legend);
    quiver3(edges_location_3d(:,1), edges_location_3d(:,2), edges_location_3d(:,3), ...
            quiver_mag*edges_tangent_3d(:,1), quiver_mag*edges_tangent_3d(:,2), quiver_mag*edges_tangent_3d(:,3), 0, ...
            'Color', colors(i, :), 'LineWidth', 0.5, 'MaxHeadSize', 0.5);

    %> Make marker size larger in the legend. This would make visualization clearer especially when multiple passes of hypothesis views are used for 3D edge sketch.
    hLeg = copyobj(h, ax);
    set(hLeg, 'XData', NaN, 'YData', NaN, 'ZData', NaN, 'MarkerSize', 20);
    hLegs = [hLegs, hLeg];
end

%> Set the plot settings
view(3);
axis equal;
axis off;
set(gcf, 'color', 'w');

%> Avoid 3D edge points vanish when zomming in
ax = gca;
ax.Clipping = "off";

%> Add a legend for each file
legend(hLegs);
hold off;
