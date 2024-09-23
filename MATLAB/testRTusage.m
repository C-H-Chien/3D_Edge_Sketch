clear all
close all
load('MATData/TO_Edges_10.mat');
load('MATData/imageArray_10.mat')
load('MATData/R_matrix_10.mat')
load('MATData/T_matrix_10.mat')
load('MATData/K_10.mat')


% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 1;
params.HYPO2_VIEW_INDX          = 6;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.99;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.001*5;
params.ICL_DATA                 = 0;
params.multiK                   = 1;

%> Convert from orientation to tangent vector for each edgel
TO_Edges_Tangents = cell(params.NUM_OF_IMGS, 1);
for i = 1:params.NUM_OF_IMGS
    TO_Edgels_Orientations = collect_third_order_edges{i,1}(:,3);
    
    Tgts = zeros(size(TO_Edgels_Orientations, 1), 2);
    for j = 1:size(TO_Edgels_Orientations, 1)
        Tgts(j,:) = [cos(TO_Edgels_Orientations(j,1)), ...
                     sin(TO_Edgels_Orientations(j,1))];
    end
    TO_Edges_Tangents{i,1} = Tgts;
end

%> Array of view indices
view_idx = [];
for i = 1:params.NUM_OF_IMGS; view_idx = [view_idx; i]; end
validation_view_indices = view_idx;
validation_view_indices([params.HYPO1_VIEW_INDX; params.HYPO2_VIEW_INDX],:) = [];

%> Fetch edgels of the two hypothesis views
TO_Edges_HYPO1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,1};
TO_Edges_HYPO2 = collect_third_order_edges{params.HYPO2_VIEW_INDX,1};
Gamma1s        = [];

hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
rows      = size(hypo_img2, 1);
cols      = size(hypo_img2, 2);

edgel_HYPO1 = TO_Edges_HYPO1(100, :);

K1 = K(:,:,params.HYPO1_VIEW_INDX);
K2 = K(:,:,params.HYPO2_VIEW_INDX);
invK1 = inv(K1);
invK2 = inv(K2);
abs_R1 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_T1 =  T_matrix(:,params.HYPO1_VIEW_INDX);
abs_C1 = -abs_R1' * abs_T1;
abs_R2 =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_T2 =  T_matrix(:,params.HYPO2_VIEW_INDX);
abs_C2 = -abs_R2' * abs_T2;

%> Revised by CH
abs_R1 = abs_R1';
abs_R2 = abs_R2';
abs_C1 = -abs_C1;
abs_C2 = -abs_C2;

R21 = abs_R2' * abs_R1;
T21 = abs_R2' * (abs_C1 - abs_C2);

T_x = @(T)[0, -T(3,1), T(2,1); T(3,1), 0, -T(1,1); -T(2,1), T(1,1), 0];
E  = T_x(T21) * R21;
F  = inv(K2)' * E * inv(K1);

e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% Calculation for epipole
epipole_met_view1 = [(e1'*R21'*T21) / (e3'*R21'*T21); (e2'*R21'*T21) / (e3'*R21'*T21)];
epipole_met_view2 = [(e1'*T21) / (e3'*T21); (e2'*T21) / (e3'*T21)];
epipole_pix_view1 = K1 * [epipole_met_view1; 1];
epipole_pix_view2 = K2 * [epipole_met_view2; 1];

coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
Apixel = coeffs(1,1);
Bpixel = coeffs(2,1);
Cpixel = coeffs(3,1);
Epipolar_Coeffs.A = Apixel;
Epipolar_Coeffs.B = Bpixel;
Epipolar_Coeffs.C = Cpixel;

x=[0, size(hypo_img2,2)];
% express the equation for x and y in pixel representation
y1=(-Apixel*x-Cpixel)/Bpixel;

figure;
cols = size(hypo_img2, 2);
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
subplot(1,2,1);
imshow(uint8(hypo_img1)); hold on;
title 'hypothesis view 1'
% y = bx + c
hypo1  = edgel_HYPO1(:,1:2)' - epipole_pix_view1(1:2,:);
b1     = hypo1(2,:)./hypo1(1,:);
c1     = edgel_HYPO1(:,2) - b1 * edgel_HYPO1(:,1);
y1_min = b1*1 + c1;
y1_max = b1*cols + c1;
line([1, cols], [y1_min, y1_max], 'Color', 'y', 'LineWidth', 0.5);
hold on;
plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',0.25);
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mo','MarkerSize',5);
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mx','MarkerSize',5);
set(gcf,'color','w');

subplot(1,2,2);
imshow(uint8(hypo_img2)); hold on;
title 'hypothesis view 2'
hold on;
plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',0.25);
hold on;
plot(x,y1,'r', 'LineWidth', 0.05);