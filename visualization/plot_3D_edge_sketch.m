%> Read from the output .txt 3D edges file and plot 3D edges as 3D points
%
%> TODO: add 3D tangents
%
%> (c) LEMS, Brown University
%> chiang-heng chien

%> Define the path from which the output 3D egdes file comes
data_folder_path = "/gpfs/data/bkimia/cchien3/3D_Edge_Sketch/outputs/";

%> Some settings which define the output file name
dataset_name = "ABC-NEF";
object_name  = "00000006";
delta        = "03";
theta        = "15";
N            = "4";
H1_index     = "6";
H2_index     = "8";

output_file_name = strcat("Gamma1s_", dataset_name, "_", object_name, "_", H1_index, "n", H2_index, ...
                          "_t32to0_delta", delta, "_theta", theta, "_N", N, ".txt");

edges_3D_file_name = fullfile(data_folder_path, output_file_name);
edges_file_read = fopen(edges_3D_file_name, 'r');
ldata = textscan(edges_file_read, '%f\t%f\t%f', 'CollectOutput', true);
edges_3d = double(ldata{1,1});

figure(2);
plot3(edges_3d(:,1), edges_3d(:,2), -edges_3d(:,3), 'Color', 'r', 'Marker', '.', 'LineStyle','none');

[az, el] = view(gca);
axis equal;
% xlabel('x');
% ylabel('y');
% zlabel('z');
axis off;
% set(gca,'FontSize',15);
set(gcf,'color','w');
