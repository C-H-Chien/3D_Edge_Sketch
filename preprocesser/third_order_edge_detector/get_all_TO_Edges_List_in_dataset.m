%> This code generates *all* third-order edges for the entire dataset
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
Dataset_Path        = '/gpfs/data/bkimia/Datasets/';
Dataset_Name        = 'ABC-NEF/';      %> ABC-NEF or DTU
Image_Folder_Name   = 'train_img/';    %> train_img for ABC-NEF, color for DTU
postfix             = '.png';          %> .png for ABC-NEF and DTU
Edges_Folder_Name   = 'Edges/';


All_object_names_file = "ABC-NEF-object-names.txt";
All_object_names_path = fullfile(Dataset_Path, Dataset_Name, All_object_names_file);
list_of_objects = importdata(All_object_names_path);

%> Settings for the Third-Order Edge Detector
All_TOED_threshs = [1, 2, 4, 8, 16, 32];
sigma = 1;
n = 1;

%> Whether to save .edg files for curve sketch
save_edg_files = 0;

for i = 1:length(list_of_objects)
    object_name = pad(string(list_of_objects(i)), 8, 'left', '0');
    fprintf(strcat("Object ", object_name, "\n"));

    [status, msg, ~] = mkdir(strcat(Dataset_Path, Dataset_Name, object_name, Edges_Folder_Name));
    disp(msg);

    for j = 1:length(All_TOED_threshs)
        
        toed_thresh = All_TOED_threshs(j);
        fprintf(strcat("TOED Threshold ", string(toed_thresh), " "));

        All_Images = dir(strcat(Dataset_Path, Dataset_Name, object_name, "/", Image_Folder_Name, '*', postfix));
        Full_Accessible_Path = [Dataset_Path, Dataset_Name, object_name, "/", Image_Folder_Name];

        for k = 1:size(All_Images, 1)
            src_Data_Path = strcat(Dataset_Path, Dataset_Name, object_name, "/", Image_Folder_Name, All_Images(k).name);
            img_ = imread(src_Data_Path);
            img_ = double(rgb2gray(img_));
            [TO_edges, ~, ~, ~] = third_order_edge_detector(img_, sigma, n, toed_thresh, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % outputs of *third_order_edge_detector*
            % TO_edges = [Subpixel_X Subpixel_Y Orientation Confidence]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            TO_Edges_Name = extractBefore(All_Images(k).name, postfix);
            yy = find(TO_edges(:,1) < 10 | TO_edges(:,1) > size(img_,1)-10 | TO_edges(:,2) < 10 | TO_edges(:,2) > size(img_,2)-10);
            TO_edges(yy,:) = [];
        
            %> Save as .edg file
            if save_edg_files == 1
                save_edg([Full_Accessible_Path, TO_Edges_Name, '.edg'], TO_edges, [size(img_, 1), size(img_, 2)]);
            end
        
            %> Save as .txt file
            img_index               = extractBefore(string(All_Images(k).name), "_");
            output_edges_file_txt   = strcat("Edge_", img_index, "_t", string(toed_thresh), ".txt");
            output_file_path        = fullfile(Dataset_Path, Dataset_Name, object_name, Edges_Folder_Name, output_edges_file_txt);
            writematrix(TO_edges, output_file_path, 'Delimiter', 'tab');
            
            %> Monitor the progress
            fprintf(". ");
        
            % figure(1);
            % src_Data_Path = strcat(Full_Accessible_Path, All_Images(i).name);   %> 01.jpg
            % img_ = imread(src_Data_Path);
            % img_ = double(rgb2gray(img_));
            % imshow(uint8(img_)); hold on;
            % plot(TO_edges(:,1), TO_edges(:,2), 'c.');
        end
        fprintf("\n");
    end
    fprintf("\n");
end
