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
    GT_Rot(:,:,i) = quat2rotm(quat)';
    GT_Transl(:,i) = -GT_Rot(:,:,i)*GT_Center(:,i);
end
K = [481.20, 0.0, 319.50; 0.0, -480.00, 239.50; 0.0, 0.0, 1.0];
invK = inv(K);


%> Below is an example of fetching two images, compute their relative pose,
%  and also pick one point on the first image, find the corresponding point
%  on the second image using the ground truth depth

%> Two example images
HYPO1_IMG_INDX = 6*19;
HYPO2_IMG_INDX = 3*19;
Example_Point_on_HYPO1_IMG = [236, 236];

%> Read Images
HYPO1_IMG_dir = strcat(FilePath_RGB_Imgs, string(HYPO1_IMG_INDX), '.png');
HYPO2_IMG_dir = strcat(FilePath_RGB_Imgs, string(HYPO2_IMG_INDX), '.png');
HYPO1_IMG = imread(HYPO1_IMG_dir);
HYPO2_IMG = imread(HYPO2_IMG_dir);

%> Get Relative Pose
abs_R1 = GT_Rot(:,:,HYPO1_IMG_INDX);
abs_R2 = GT_Rot(:,:,HYPO2_IMG_INDX);
abs_C1 = GT_Center(:,HYPO1_IMG_INDX);
abs_C2 = GT_Center(:,HYPO2_IMG_INDX);
R21    = abs_R2 * abs_R1';
T21    = abs_R2 * (abs_C1 - abs_C2);
T_x = @(T)[0,      -T(3,1),  T(2,1); ...
           T(3,1),  0,      -T(1,1); ...
          -T(2,1),  T(1,1),  0];
E21 = T_x(T21) * R21;
F21 = invK' * E21 * invK;

%> Compute the coefficients of the epipolar line
coeffs = F21 * [Example_Point_on_HYPO1_IMG'; 1];
Epipolar_Coeffs.A = coeffs(1,1);
Epipolar_Coeffs.B = coeffs(2,1);
Epipolar_Coeffs.C = coeffs(3,1);

%> Read Depths
HYPO1_Depth_dir = strcat(FilePath_Depths, string(HYPO1_IMG_INDX), '.png');
HYPO2_Depth_dir = strcat(FilePath_Depths, string(HYPO2_IMG_INDX), '.png');
HYPO1_Depth = double(imread(HYPO1_Depth_dir));
HYPO2_Depth = double(imread(HYPO2_Depth_dir));
HYPO1_Depth = HYPO1_Depth ./ 5000;
HYPO2_Depth = HYPO2_Depth ./ 5000;

%> Get ground truth correspondence point gamma2 on the second image from
%  gamma1 on the first image
e3 = [0;0;1];
gamma1 = invK * [Example_Point_on_HYPO1_IMG'; 1];
rho1   = HYPO1_Depth(Example_Point_on_HYPO1_IMG(2), Example_Point_on_HYPO1_IMG(1));
gamma2 = (rho1*R21*gamma1 + T21)/((e3'*R21*gamma1)*rho1 + e3'*T21);
pt2    = K*gamma2;

%> Display
figure;
imshow(HYPO1_IMG); hold on;
plot(Example_Point_on_HYPO1_IMG(1), Example_Point_on_HYPO1_IMG(2), 'ms');

pause(0.5);
figure;
imshow(HYPO2_IMG); hold on;
cols = size(HYPO2_IMG, 2);
yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
line([1, cols], [yMin, yMax], 'Color', 'y', 'LineWidth', 2); hold on;
plot(pt2(1), pt2(2), 'bs'); hold on;
set(gcf,'color','w');



