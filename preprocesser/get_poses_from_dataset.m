close all;
clear;

import yaml.loadFile

%> Set 1 to Write_to_Projmatrix_Files is writing projection matrices as a 
%  series of .projmatrix files, if necessary. This is mainly used for multiview curve sketch
Write_to_Projmatrix_Files = 0;
Write_rotations_and_translations_in_files = 1;
Check_Camera_Rot_Transl = 0;

%> Define the path of the dataset
media_storage = "/gpfs/data/bkimia/Datasets/";
dataset_name = "ABC-NEF/";  %> ABC-NEF or DTU
object_tag_name = "00008605/";
mfiledir = strcat(media_storage, dataset_name);
if strcmp(dataset_name, "ABC-NEF/")
    json_file_name = "transforms_train.json";
    extra_folder_name = "train_img";
elseif strcmp(dataset_name, "DTU/")
    json_file_name = "meta_data.json";
    extra_folder_name = "color";
elseif strcmp(dataset_name, "Replica/")
    json_file_name = "meta_data.json";
    extra_folder_name = "";
else
    error("Invalid dataset name!")
end

%> yml file containing matrices of all the views
ymlPath = fullfile(mfiledir, object_tag_name, json_file_name);

data = yaml.loadFile(ymlPath);
view_matrix = data.frames;
K_view = {};
RT_view = {};
for i = 1:size(view_matrix, 2)
    %> get intrinsic matrix K
    if strcmp(dataset_name, "ABC-NEF/")
        tmp = view_matrix{i}.camera_intrinsics;
    elseif strcmp(dataset_name, "DTU/")
        tmp = view_matrix{i}.intrinsics;
    elseif strcmp(dataset_name, "Replica/")
        tmp = view_matrix{i}.intrinsics;
    end
    K_(1, :) = cell2mat(tmp{1});
    K_(2, :) = cell2mat(tmp{2});
    K_(3, :) = cell2mat(tmp{3});
    if size(K_,2) > 3
        K = K_(1:3,1:3);
    else
        K = K_;
    end
    
    %> get camera to world ground-truth pose (R, T)
    if strcmp(dataset_name, "ABC-NEF/")
        tmp = view_matrix{i}.transform_matrix;
    elseif strcmp(dataset_name, "DTU/")
        tmp = view_matrix{i}.camtoworld;
    elseif strcmp(dataset_name, "Replica/")
        tmp = view_matrix{i}.camtoworld;
    end
    RT(1, :) = cell2mat(tmp{1});
    RT(2, :) = cell2mat(tmp{2});
    RT(3, :) = cell2mat(tmp{3});
    RT(4, :) = cell2mat(tmp{4});
    
    K_view{end + 1} = K;

    if strcmp(dataset_name, "DTU/")
        RT(:,2) = RT(:,2).*(-1);
        RT(:,3) = RT(:,3).*(-1);
    end
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
    
    %> Projection matrix. Before multiplying K, let x be -x
    % if strcmp(dataset_name, "ABC-NEF/")
    trans = [-1 0 0; 0 1 0; 0 0 1] * trans;
    % end
    projMat = K * trans;

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
        if strcmp(dataset_name, "ABC-NEF/") || strcmp(dataset_name, "DTU/")
            file_name = sprintf("%d_colors.projmatrix", i - 1);
        elseif strcmp(dataset_name, "Replica/")
            file_name = pad(string(i-1), 6, 'left', '0');
            file_name = strcat(file_name, "_rgb.projmatrix");
        end
        path = fullfile(mfiledir, object_tag_name, extra_folder_name, file_name);
        writematrix(proj_matrix, path, "FileType", "text", "Delimiter", " ");
    end
end

%> Write rotations and translations for the use of edge sketch
if Write_rotations_and_translations_in_files == 1
    [status, msg, ~] = mkdir(strcat(media_storage, dataset_name, object_tag_name, "RnT"));
    disp(msg);
    rot_path = fullfile(mfiledir, object_tag_name, "RnT/R_matrix.txt");
    transl_path = fullfile(mfiledir, object_tag_name, "RnT/T_matrix.txt");
    writematrix(rotation_matrix_by_view, rot_path, "FileType", "text", "Delimiter", "\t");
    writematrix(translation_vector_by_view, transl_path, "FileType", "text", "Delimiter", "\t");
end

if Check_Camera_Rot_Transl == 1
    R1 = rotation_matrix_by_view(1:3,:);
    T1 = translation_vector_by_view(1:3,:);

    R2 = rotation_matrix_by_view(4:6,:);
    T2 = translation_vector_by_view(4:6,:);

    img0_file_name = "0_colors";
    img1_file_name = "1_colors";

    img0 = imread(fullfile(mfiledir, object_tag_name, extra_folder_name, strcat(img0_file_name, ".png")));
    img1 = imread(fullfile(mfiledir, object_tag_name, extra_folder_name, strcat(img1_file_name, ".png")));

    Rel_R = R2 * R1';
    Rel_T = -R2 * R1' * T1 + T2;
    skew_T = @(T)[0, -T(3,1), T(2,1); T(3,1), 0, -T(1,1); -T(2,1), T(1,1), 0];
    E = skew_T(Rel_T) * Rel_R;
    F = inv(K)' * E * inv(K);

    cols = size(img1, 2);
    test_pt = [250; 150; 1];
    a = F(1,:)*test_pt;
    b = F(2,:)*test_pt;
    c = F(3,:)*test_pt;
    yMin = -c/b;
    yMax = (-c - a*cols) / b;

    figure(1);
    imshow(img0); hold on;
    plot(test_pt(1), test_pt(2), 'cs');
    hold off;
    pause(0.5);
    figure(2);
    imshow(img1); hold on;
    line([1, cols], [yMin, yMax], 'Color', 'c', 'LineWidth', 2);

end


