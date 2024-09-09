close all;
clear;

import yaml.loadFile

%> Set 1 to Write_to_Projmatrix_Files is writing projection matrices as a 
%  series of .projmatrix files, if necessary. This is mainly used for multiview curve sketch
Write_to_Projmatrix_Files = 1;
Write_rotations_and_translations_in_files = 1;

%> curve points sampled from ABC dataset's parametrized representation
% input_curves = load("curves.mat").curve_points;

%> yml file containing matrices of all the views
media_storage = "/media/chchien/843557f5-9293-49aa-8fb8-c1fc6c72f7ea/";
dataset_name = "ABC-NEF/";
object_tag_name = "00000162";
mfiledir = strcat(media_storage, "/home/chchien/datasets/", dataset_name);
ymlPath = fullfile(mfiledir, object_tag_name, "transforms_train.json");

data = yaml.loadFile(ymlPath);
view_matrix = data.frames;
K_view = {};
RT_view = {};
for i = 1:size(view_matrix, 2)
    %> get K
    tmp = view_matrix{i}.camera_intrinsics;
    K(1, :) = cell2mat(tmp{1});
    K(2, :) = cell2mat(tmp{2});
    K(3, :) = cell2mat(tmp{3});
    %> get RT
    tmp = view_matrix{i}.transform_matrix;
    RT(1, :) = cell2mat(tmp{1});
    RT(2, :) = cell2mat(tmp{2});
    RT(3, :) = cell2mat(tmp{3});
    RT(4, :) = cell2mat(tmp{4});
    
    K_view{end + 1} = K;
    RT_view{end + 1} = RT;
end

viewCnt = size(K_view, 2);
projection_matrix_by_view = {};
rotation_matrix_by_view = zeros(3*viewCnt, 3);
translation_vector_by_view = zeros(3*viewCnt, 1);

for i = 1:viewCnt
    %> image of the object from the view i
    % imgPath = fullfile(mfiledir, object_tag_name, "/train_img/", sprintf("%d_colors.png", i - 1));
    %> the dataset's matrix is from camera to world. Use inverse to get the
    %matrix from world to camera
    trans = inv(RT_view{i});
    trans = trans(1:3, :);

    K = K_view{i};
    
    %> Projection matrix. Before multiply K, let x be -x
    projMat = K * [-1 0 0; 0 1 0; 0 0 1] * trans;

    %> Rotation matrix.
    rotation_matrix_by_view(3*(i-1)+1:3*(i-1)+3, :) = trans(1:3, 1:3);

    %> Translation vector
    translation_vector_by_view(3*(i-1)+1:3*(i-1)+3, 1) = trans(1:3, 4);
    
    %> visualization. Curve points projection should overlap the edge of
    %the object 
    % imshow(imread(imgPath));
    % hold on;
    % for j = 1:size(input_curves, 2)
    %     c = input_curves{j};
    %     c = [c'; ones(1, size(c, 1))];
    %     proj = projMat * c;
    %     proj = proj(1:2, :) ./ proj(3, :);
    %     plot(proj(1, :), proj(2, :));
    %     hold on
    % end
    % hold off;

    projection_matrix_by_view{end + 1} = projMat;

    fprintf("-");
end
fprintf("\n");

%> save the matrix converting result
% save_path = fullfile(mfiledir, "projection.mat");
% save(save_path, "projection_matrix_by_view");

%> Write projection matrix for the use of curve sketch
if Write_to_Projmatrix_Files == 1
    for i = 1:length(projection_matrix_by_view)
        proj_matrix = projection_matrix_by_view{i};
        path = fullfile(mfiledir, object_tag_name, "train_img", sprintf("%d_colors.projmatrix", i - 1));
        writematrix(proj_matrix, path, "FileType", "text", "Delimiter", " ");
    end
end

%> Write rotations and translations for the use of edge sketch
if Write_rotations_and_translations_in_files == 1
    rot_path = fullfile(mfiledir, object_tag_name, "R_matrix.txt");
    transl_path = fullfile(mfiledir, object_tag_name, "T_matrix.txt");
    writematrix(rotation_matrix_by_view, rot_path, "FileType", "text", "Delimiter", "\t");
    writematrix(translation_vector_by_view, transl_path, "FileType", "text", "Delimiter", "\t");
end
