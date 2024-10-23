%> Read from the output .txt 3D edges file and plot 3D edges as 3D points
%
%> TODO: add 3D tangents
%
%> (c) LEMS, Brown University
%> chiang-heng chien

% Define the path from which the output 3D edges file comes
data_folder_path = "/gpfs/data/bkimia/zqiwu/3D_Edge_Sketch/outputs/";

% Some settings which define the output file name
dataset_name = "ABC-NEF";
object_name  = "00000006";
delta        = "03";
theta        = "15";
N            = "4";

% Define the first set of indices (6 and 8)
H1_index_1 = "6";
H2_index_1 = "8";

output_file_name_1 = strcat("Gamma1s_", dataset_name, "_", object_name, "_", H1_index_1, "n", H2_index_1, ...
                          "_t32to0_delta", delta, "_theta", theta, "_N", N, ".txt");

% Define the second set of indices (28 and 47)
H1_index_2 = "28";
H2_index_2 = "47";

output_file_name_2 = strcat("Gamma1s_", dataset_name, "_", object_name, "_", H1_index_2, "n", H2_index_2, ...
                          "_t32to0_delta", delta, "_theta", theta, "_N", N, ".txt");

% Read the first file
edges_file_read_1 = fopen(fullfile(data_folder_path, output_file_name_1), 'r');
ldata_1 = textscan(edges_file_read_1, '%f\t%f\t%f', 'CollectOutput', true);
edges_3d_1 = double(ldata_1{1,1});
fclose(edges_file_read_1);

% Read the second file
edges_file_read_2 = fopen(fullfile(data_folder_path, output_file_name_2), 'r');
ldata_2 = textscan(edges_file_read_2, '%f\t%f\t%f', 'CollectOutput', true);
edges_3d_2 = double(ldata_2{1,1});
fclose(edges_file_read_2);

% Plot the first set of edges
figure(2);
hold on;
plot3(edges_3d_1(:,1), edges_3d_1(:,2), -edges_3d_1(:,3), 'Color', 'r', 'Marker', '.', 'LineStyle', 'none');

% Plot the second set of edges
plot3(edges_3d_2(:,1), edges_3d_2(:,2), -edges_3d_2(:,3), 'Color', 'b', 'Marker', '.', 'LineStyle', 'none');

% Set the plot settings
[az, el] = view(gca);
axis equal;
axis off;
set(gcf, 'color', 'w');
legend('H1=6, H2=8', 'H1=28, H2=47');
