close all;
clear all;

%> Dataset and sequence paths
datasetFolder = 'F:\ICL-NUIM\';
sequenceName  = 'traj1_frei_png.tar\';

%> Data paths and file names
FileName_GT_Poses = 'traj1.gt.freiburg';
FilePath_RGB_Imgs = strcat(datasetFolder, sequenceName, 'rgb\');
FilePath_Depths   = strcat(datasetFolder, sequenceName, 'depth\');

%> Get Ground Truth Poses (Assume there are N poses.)
%  GT_Transl: Ground Truth Translation   (3xN)
%  GT_Rot   : Ground Truth Rotation      (3x3xN)
%  GT_Center: Ground Truth Camera Center (3xN)
%  K        : Camera Calibration Matrix  (3x3)
GT_Poses_dir  = strcat(datasetFolder, sequenceName, FileName_GT_Poses);
GT_Poses_file = fopen(GT_Poses_dir, 'r');
ldata         = textscan( GT_Poses_file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f', 'CollectOutput', true );
GT_C          = ldata{1,1}(:,2:4);
GT_Q          = ldata{1,1}(:,5:8);
GT_Transl     = zeros(3, size(GT_C, 1));
GT_Rot        = zeros(3, 3, size(GT_C, 1));
GT_Center     = GT_C';
for i = 1:size(GT_C, 1)
    quat = [GT_Q(i,4), GT_Q(i,1), GT_Q(i,2), GT_Q(i,3)];
    GT_Rot(:,:,i) = quat2rotm(quat);
    GT_Transl(:,i) = -GT_Rot(:,:,i)*GT_Center(:,i);
end
K = [481.20, 0.0, 319.50; 0.0, -480.00, 239.50; 0.0, 0.0, 1.0];
invK = inv(K);

imageArray = {;,;,;};
R_matrix   = zeros(3,3,50);
T_matrix   = zeros(3,1,50);
% get imageArray and other required data
figure(1)
for(i = 1:50)
    % 
    KeyFrame_dir    = strcat(FilePath_RGB_Imgs, string(i*19), '.png');
    imageArray{i}   = imread(KeyFrame_dir);
    R_matrix(:,:,i) = GT_Rot(:,:,i*19);
    T_matrix(:,:,i) = GT_Transl(:,i*19);
    imshow(imageArray{i});
end

R_matrix1   = zeros(3,3,25);
T_matrix1   = zeros(3,1,25);
for(i = 1:25)
    R_matrix1(:,:,i) = R_matrix(:,:,i+25);
    T_matrix1(:,:,i) = T_matrix(:,:,i+25);
end