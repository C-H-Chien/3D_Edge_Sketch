
%> Third-Order Edge Detector Credit:
%> Paper: Kimia, Benjamin B., Xiaoyan Li, Yuliang Guo, and Amir Tamrakar. 
%         "Differential geometry in edge detection: accurate estimation of 
%         position, orientation and curvature." IEEE transactions on 
%         pattern analysis and machine intelligence 41, no. 7 (2018): 
%         1573-1586.
%> Implementations: (1) https://github.com/yuliangguo/Differential_Geometry_in_Edge_Detection
%                   (2) https://github.com/C-H-Chien/Third-Order-Edge-Detector

% mfiledir = fileparts(mfilename('fullpath'));
%> Path to where the dataset is
Dataset_Path        = '/oscar/data/bkimia/Datasets/';
Dataset_Name        = 'ABC-NEF/';
Image_Folder_Name   = 'train_img/';
postfix             = '.png';
Edges_Folder_Name   = 'Edges/';
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

%> Whether to save .edg files for curve sketch
save_edg_files = 0;

%> Settings for the Third-Order Edge Detector
All_Threshs = [1, 2, 4, 8, 16, 32];
% thresh = 1;
sigma = 1;
n = 1;
format long;

%> Loop over all folders under the ABC-NEF dataset
for fi = 1:length(All_Scene_Name_Folders)
    Scene_Name = All_Scene_Name_Folders(fi);
    All_Images = dir(strcat(Dataset_Path, Dataset_Name, Scene_Name, Image_Folder_Name, '*', postfix));
    Full_Accessible_Path = [Dataset_Path, Dataset_Name, Scene_Name, Image_Folder_Name];
    fprintf("Extracting TO edges for scene %s\n", Scene_Name);

    for i = 1:size(All_Images, 1)
        src_Data_Path = strcat(Dataset_Path, Dataset_Name, Scene_Name, Image_Folder_Name, All_Images(i).name);
        img_ = imread(src_Data_Path);
        img_ = double(rgb2gray(img_));

        for j = 1:length(All_Threshs)
            thresh = All_Threshs(j);
            [TO_edges, ~, ~, ~] = third_order_edge_detector(img_, sigma, n, thresh, 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % outputs of *third_order_edge_detector*
            % TO_edges = [Subpixel_X Subpixel_Y Orientation Confidence]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            TO_Edges_Name = extractBefore(All_Images(i).name, postfix);
            yy = find(TO_edges(:,1) < 10 | TO_edges(:,1) > size(img_,1)-10 | TO_edges(:,2) < 10 | TO_edges(:,2) > size(img_,2)-10);
            TO_edges(yy,:) = [];
        
            %> Save as .edg file
            if save_edg_files == 1
                save_edg([Full_Accessible_Path, TO_Edges_Name, '.edg'], TO_edges, [size(img_, 1), size(img_, 2)]);
            end
        
            %> Save as .txt file
            img_index               = extractBefore(string(All_Images(i).name), "_");
            output_edges_file_txt   = strcat("Edge_", img_index, "_t", string(thresh), ".txt");
            output_file_path        = fullfile(Dataset_Path, Dataset_Name, Scene_Name, Edges_Folder_Name, output_edges_file_txt);
            writematrix(TO_edges, output_file_path, 'Delimiter', 'tab');
            
            %> Monitor the progress
            % fprintf(". ");
        
            % figure(1);
            % src_Data_Path = strcat(Full_Accessible_Path, All_Images(i).name);   %> 01.jpg
            % img_ = imread(src_Data_Path);
            % img_ = double(rgb2gray(img_));
            % imshow(uint8(img_)); hold on;
            % plot(TO_edges(:,1), TO_edges(:,2), 'c.');
        end
        fprintf(". ");
    end
    fprintf("\n");
end

fprintf("\n");

