load('TO_Edges_ICL-NUIM_ofkt1_26to50.mat');
load('imageArray_ICL-NUIM_ofkt1_26to50.mat')
load('R_matrix_ICL-NUIM_ofkt1_26to50.mat')
load('T_matrix_ICL-NUIM_ofkt1_26to50.mat')
load('K_ICL-NUIM_ofkt1.mat')
icl =1;
% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 44-25;
params.HYPO2_VIEW_INDX          = 31-25;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.9;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.001*5;
invK = inv(K);

TO_Edges_HYPO1 = collect_third_order_edges{params.HYPO1_VIEW_INDX,1};
TO_Edges_HYPO2 = collect_third_order_edges{params.HYPO2_VIEW_INDX,1};
% choose the edge pair
edgel_HYPO1  = TO_Edges_HYPO1(12368,:);
edgels_HYPO2 = TO_Edges_HYPO2(6796,:);


% choose the view of validation
VALID_INDX = 6;
VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
% get edges in the validation view
TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};


abs_R1 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C1 = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' * T_matrix(:,params.HYPO1_VIEW_INDX);
abs_R2 =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C2 = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' * T_matrix(:,params.HYPO2_VIEW_INDX);
abs_R3 =  R_matrix(:,:,VALID_INDX);
abs_C3 = -R_matrix(:,:,VALID_INDX)' * T_matrix(:,VALID_INDX);

R21 = abs_R2' * abs_R1;
T21 = abs_R2' * (abs_C1 - abs_C2);
R31 = abs_R3' * abs_R1;
T31 = abs_R3' * (abs_C1 - abs_C3);

e1 = [1;0;0];
e3 = [0;0;1];

point1            = [edgel_HYPO1(:,1:2)'; 1];
gamma1            = invK*[edgel_HYPO1(:,1:2)'; 1];
tangent1          = [cos(edgel_HYPO1(:,3)); ...
    sin(edgel_HYPO1(:,3))];
pt1_tgt_to_pixels = point1(1:2,1) + tangent1;
pt1_tgt_to_meters = invK * [pt1_tgt_to_pixels; 1];
tgt1_meters       = pt1_tgt_to_meters - gamma1;
edge_pos_gamma3   = zeros(size(edgels_HYPO2,1), 2);
edge_tgt_gamma3   = zeros(size(edgels_HYPO2,1), 2);


%> fetch hypothesis pair in the second hypothesis view
point2 = [edgels_HYPO2(1, 1:2)'; 1];
gamma2 = invK*[edgels_HYPO2(1, 1:2)'; 1];

%> Get the position of the reprojected point
rho1   = (e1'*T21 - (e3' * T21) * (e1'*gamma2)) / ...
    ((e3'*R21*gamma1) * (e1'*gamma2) - (e1'*R21*gamma1));
rho3   = rho1 * (e3'*R31*gamma1) + e3'*T31;
gamma3 = 1 / rho3 * (R31*rho1*gamma1 + T31);
point3 = K * gamma3;

edge_pos_gamma3(1, :) = point3(1:2, 1)';

%> Get the orientation of the reprojected point
tangent2          = [cos(edgels_HYPO2(1,3)); ...
    sin(edgels_HYPO2(1,3))];
pt2_tgt_to_pixels = point2(1:2,1) + tangent2;
pt2_tgt_to_meters = invK * [pt2_tgt_to_pixels; 1];
tgt2_meters       = pt2_tgt_to_meters - gamma2;

P1   = gamma1;
P2   = gamma2;
t1   = tgt1_meters;
t2   = tgt2_meters;
n1   = cross(t1, P1);
n2   = R21' * cross(t2, P2);
T_v1 = cross(n1,n2) ./ norm(cross(n1,n2));

T_v3 = R31 * T_v1;
t_v3 = T_v3 - (e3'*T_v3)*gamma3;
t_v3 = t_v3 ./ norm(t_v3);

if(icl == 1)
    t_v3(2,1) = -1*t_v3(2,1);
end

tangent3              = t_v3(1:2,1);
edge_tgt_gamma3(1,:) = tangent3';

% visualize result
len = 2;
figure;
subplot(1,2,1);
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
imshow(uint8(hypo_img1)); 
hold on;
plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'mx', 'MarkerSize',7, 'LineWidth', 2);
x1 = edgel_HYPO1(1,1)+len*tangent1(1,1);
x2 = edgel_HYPO1(1,1)-len*tangent1(1,1);
y1 = edgel_HYPO1(1,2)+len*tangent1(2,1);
y2 = edgel_HYPO1(1,2)-len*tangent1(2,1);
hold on;
plot([x1 x2], [y1 y2], 'm', 'LineWidth', 1);
title 'hypothesis view 1'
subplot(1,2,2);
imshow(uint8(hypo_img2)); 
hold on;
plot(TO_Edges_HYPO2(:, 1), TO_Edges_HYPO2(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
plot(edgels_HYPO2(1, 1), edgels_HYPO2(1, 2), 'mx', 'MarkerSize',7, 'LineWidth', 2);
x1 = edgels_HYPO2(1,1)+len*tangent2(1,1);
x2 = edgels_HYPO2(1,1)-len*tangent2(1,1);
y1 = edgels_HYPO2(1,2)+len*tangent2(2,1);
y2 = edgels_HYPO2(1,2)-len*tangent2(2,1);
hold on;
plot([x1 x2], [y1 y2], 'm', 'LineWidth', 1);
title 'hypothesis view 2'

figure;
imshow(uint8(VALID_IMG)); 
hold on;
plot(TO_Edges_VALID(:, 1), TO_Edges_VALID(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
plot(edge_pos_gamma3(1, 1), edge_pos_gamma3(1, 2), 'gx', 'MarkerSize',7, 'LineWidth', 2);
x1 = edge_pos_gamma3(1,1)+len*tangent3(1,1);
x2 = edge_pos_gamma3(1,1)-len*tangent3(1,1);
y1 = edge_pos_gamma3(1,2)+len*tangent3(2,1);
y2 = edge_pos_gamma3(1,2)-len*tangent3(2,1);
hold on;
plot([x1 x2], [y1 y2], 'g', 'LineWidth', 1);
title 'validation view'