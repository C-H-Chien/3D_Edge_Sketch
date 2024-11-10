close all;
clear;

import yaml.loadFile

%> Set 1 to Write_to_Projmatrix_Files is writing projection matrices as a 
%  series of .projmatrix files, if necessary. This is mainly used for multiview curve sketch
Write_to_Projmatrix_Files = 0;
Write_rotations_and_translations_in_files = 1;

%> curve points sampled from ABC dataset's parametrized representation
% input_curves = load("curves.mat").curve_points;

%> yml file containing matrices of all the views
media_storage   = "/oscar/data/bkimia/Datasets/";
dataset_name    = "ABC-NEF/";
% object_tag_name = "00000568";
RnT_folder_name = "RnT/";
mfiledir = strcat(media_storage, dataset_name);

All_Scene_Name_Folders = ["00000077/", "00000858/", "00001519/", "00001911/", "00002532/", "00003138/", "00003619/", ...
                        "00004605/", "00005109/", "00005894/", "00006645/", "00007551/", "00008609/", "00009579/", "00000146/", ...
                        "00000871/", "00001559/", "00001970/", "00002633/", "00003173/", "00003676/", "00004698/", "00005139/", ...
                        "00005896/", "00006730/", "00007617/", "00008750/", "00009581/", "00000168/", "00000952/", "00001567/", ...
                        "00002030/", "00002718/", "00003178/", "00003717/", "00004733/", "00005222/", "00005898/", "00006845/", ...
                        "00007644/", "00009016/", "00009683/", "00000325/", "00001164/", "00001614/", "00002063/", "00002750/", ...
                        "00003280/", "00003823/", "00004857/", "00005225/", "00006007/", "00007025/", "00007809/", "00009017/", ...
                        "00009685/", "00000369/", "00001204/", "00001766/", "00002211/", "00002852/", "00003337/", "00003884/", ...
                        "00004888/", "00005787/", "00006053/", "00007042/", "00007871/", "00009045/", "00009741/", "00000699/", ...
                        "00001225/", "00001776/", "00002345/", "00002906/", "00003525/", "00004054/", "00004926/", "00005788/", ...
                        "00006248/", "00007128/", "00008100/", "00009127/", "00009779/", "00000710/", "00001313/", "00001777/", ...
                        "00002412/", "00002954/", "00003558/", "00004383/", "00005084/", "00005793/", "00006394/", "00007192/", ...
                        "00008231/", "00009239/", "00009886/", "00000843/", "00001314/", "00001855/", "00002524/", "00003014/", ...
                        "00003598/", "00004451/", "00005087/", "00005893/", "00006464/", "00007202/", "00008605/", "00009270/"];

for fi = 1:length(All_Scene_Name_Folders)

    object_tag_name = All_Scene_Name_Folders(fi);
    ymlPath = fullfile(mfiledir, object_tag_name, "transforms_train.json");
    fprintf("Converting camera absolute poses for scene %s\n", object_tag_name);

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
    
        %> The transformation matrix also has to make x become -x
        trans = [-1 0 0; 0 1 0; 0 0 1] * trans;
    
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
        rot_path = fullfile(mfiledir, object_tag_name, RnT_folder_name, "R_matrix.txt");
        transl_path = fullfile(mfiledir, object_tag_name, RnT_folder_name, "T_matrix.txt");
        writematrix(rotation_matrix_by_view, rot_path, "FileType", "text", "Delimiter", "\t");
        writematrix(translation_vector_by_view, transl_path, "FileType", "text", "Delimiter", "\t");
    end
end
