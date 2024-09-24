
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
Dataset_Path = '/oscar/data/bkimia/zqiwu/3D_Edge_Sketch/datasets/';
Dataset_Name = 'ABC-NEF/';      %> ABC-NEF or Replica
Scene_Name = '00000006/';          %> 00009779
Image_Folder_Name = 'train_img/';   %> train_img for ABC-NEF, color for Replica
postfix = '.png';               %> .png for ABC-NEF, .jpg for Replica
All_Images = dir(strcat(Dataset_Path, Dataset_Name, Scene_Name, Image_Folder_Name, '*', postfix));
Full_Accessible_Path = [Dataset_Path, Dataset_Name, Scene_Name, Image_Folder_Name];

%> Settings for the Third-Order Edge Detector
thresh = 2;
sigma = 1;
n = 1;
format long;

for i = 1:size(All_Images, 1)
    src_Data_Path = strcat(Dataset_Path, Dataset_Name, Scene_Name, Image_Folder_Name, All_Images(i).name);
    img_ = imread(src_Data_Path);
    img_ = double(rgb2gray(img_));
    [TO_edges, ~, ~, ~] = third_order_edge_detector(img_, sigma, n, thresh, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % outputs of *third_order_edge_detector*
    % TO_edges = [Subpixel_X Subpixel_Y Orientation Confidence]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TO_Edges_Name = extractBefore(All_Images(i).name, postfix);
    yy = find(TO_edges(:,1) < 10 | TO_edges(:,1) > size(img_,1)-10 | TO_edges(:,2) < 10 | TO_edges(:,2) > size(img_,2)-10);
    TO_edges(yy,:) = [];

    %> Save as .edg file
    save_edg([Full_Accessible_Path, TO_Edges_Name, '.edg'], TO_edges, [size(img_, 1), size(img_, 2)]);

    %> Save as .txt file
    img_index               = extractBefore(string(All_Images(i).name), "_");
    output_edges_file_txt   = strcat("Edge_", img_index, "_t", string(thresh), ".txt");
    output_file_path        = fullfile(Dataset_Path, Dataset_Name, Scene_Name, Image_Folder_Name, output_edges_file_txt);
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

%% 
%> An Example of super-imposing third-order edges on an image
figure;
src_Data_Path = strcat(Full_Accessible_Path, All_Images(1).name);   %> 01.jpg
img_ = imread(src_Data_Path);
img_ = double(rgb2gray(img_));
[TO_edges, ~, ~, ~] = third_order_edge_detector(img_, sigma, n, thresh, 1);
imshow(uint8(img_)); hold on;
plot(TO_edges(:,1), TO_edges(:,2), 'c.');
