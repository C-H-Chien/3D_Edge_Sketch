% Epipolar wedges method is used in this main code.
clear all; 
close all;
Nview = 4;
deltastr = '03';
thetaangle = 15;
fileID1 = fopen('GenData/pairededge_ABC0006_6n8_t32to0excludehypo1n2_delta03_theta15_N4.txt','r');
% fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
formatSpec = '%f';
A1 = fscanf(fileID1,formatSpec);
B1 = reshape(A1,50,size(A1,1)/50);
B1=B1+1;
paired_edge1 = B1';

load('MATData/TO_Edges_ABC00000006.mat');
load('MATData/imageArray_ABC00000006.mat')
load('MATData/R_matrix_ABC00000006.mat')
load('MATData/T_matrix_ABC00000006.mat')
load('MATData/K_ABC00000006.mat')

% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 7;
params.HYPO2_VIEW_INDX          = 9;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 15;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.ICL_DATA                 = 0;
params.multiK                   = 1;
params.syn                      = 0;
params.cols                     = size(double(rgb2gray(imageArray{1})), 2);
params.rows                     = size(double(rgb2gray(imageArray{1})), 1);
params.delta                    = 0.05;
params.circle                   = 55;
params.parallelangle            = 15;
ref_idx = params.HYPO1_VIEW_INDX;
if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
else
    invK = inv(K);
end

%> Fetch edgels of the two hypothesis views
% load('F:\T-less\edges50_1.mat');
TO_Edges_HYPO1  = collect_third_order_edges{params.HYPO1_VIEW_INDX,6};
TO_Edges_HYPO2  = collect_third_order_edges{params.HYPO2_VIEW_INDX,6};
Gamma1s_1         = [];
GammaTangent_3D_1 = [];

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

for(i = 1 : size(paired_edge1,1))
    if((size(find(paired_edge1(i,:)>0),2)-2)<Nview)
        continue;
    end
    edgel_HYPO1 = TO_Edges_HYPO1(paired_edge1(i,1), :);
    edgel_HYPO2 = TO_Edges_HYPO2(paired_edge1(i,2), :);
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    [edgels_HYPO2_corrected, ...
     edgels_HYPO1_corrected, e_i] = edgelsHYPO2correct(Epipolar_Coeffs, ...
                                                       edgel_HYPO2, ...
                                                       edgel_HYPO1, ...
                                                       R_matrix, ...
                                                       T_matrix, ...
                                                       params, K);
    edgels_HYPO2_final = edgels_HYPO2_corrected;
    edgels_HYPO1_final = edgels_HYPO1_corrected;
    if(params.multiK == 1)
                    pt1 = invK1 * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK2 * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                else
                    pt1 = invK * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                end
    pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
    Gamma = linearTriangulation(2, pts_meters, R21, T21);
    GammaTgt_3D     = get3DOrientation(edgels_HYPO1_final,edgels_HYPO2_final, K, params, R21);
    Gamma1s_1 = [Gamma1s_1, Gamma];
    GammaTangent_3D_1 = [GammaTangent_3D_1, GammaTgt_3D];
%     Gamma = Gamma*1000;
end


I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
                 T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_orin_1 = -R21afterBA'*(Gamma1s_1-T21afterBA);
Gamma1s_ore_orin1_1 = -R21afterBA'*(GammaTangent_3D_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID2 = fopen('GenData/pairededge_ABC0006_22n26_t32to0excludehypo1n2_delta03_theta15_N4.txt','r');
% fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
formatSpec = '%f';
A2 = fscanf(fileID2,formatSpec);
B2 = reshape(A2,50,size(A2,1)/50);
B2=B2+1;
paired_edge2 = B2';

params.HYPO1_VIEW_INDX          = 23;
params.HYPO2_VIEW_INDX          = 27;

if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
else
    invK = inv(K);
end

%> Fetch edgels of the two hypothesis views
% load('F:\T-less\edges50_1.mat');
TO_Edges_HYPO1  = collect_third_order_edges{params.HYPO1_VIEW_INDX,6};
TO_Edges_HYPO2  = collect_third_order_edges{params.HYPO2_VIEW_INDX,6};
Gamma1s_2         = [];
GammaTangent_3D_2 = [];

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

for(i = 1 : size(paired_edge2,1))
    if((size(find(paired_edge1(i,:)>0),2)-2)<Nview)
        continue;
    end
    edgel_HYPO1 = TO_Edges_HYPO1(paired_edge2(i,1), :);
    edgel_HYPO2 = TO_Edges_HYPO2(paired_edge2(i,2), :);
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    [edgels_HYPO2_corrected, ...
     edgels_HYPO1_corrected, e_i] = edgelsHYPO2correct(Epipolar_Coeffs, ...
                                                       edgel_HYPO2, ...
                                                       edgel_HYPO1, ...
                                                       R_matrix, ...
                                                       T_matrix, ...
                                                       params, K);
    edgels_HYPO2_final = edgels_HYPO2_corrected;
    edgels_HYPO1_final = edgels_HYPO1_corrected;
    if(params.multiK == 1)
                    pt1 = invK1 * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK2 * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                else
                    pt1 = invK * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                end
    pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
    Gamma = linearTriangulation(2, pts_meters, R21, T21);
    GammaTgt_3D     = get3DOrientation(edgels_HYPO1_final,edgels_HYPO2_final, K, params, R21);
    Gamma1s_2 = [Gamma1s_2, Gamma];
    GammaTangent_3D_2 = [GammaTangent_3D_2, GammaTgt_3D];
end


abs_R1afterBA =  R_matrix(:,:,ref_idx);
abs_C1afterBA = -R_matrix(:,:,ref_idx)' *...
    T_matrix(:,ref_idx);
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_2trans = R21afterBA'*(Gamma1s_2-T21afterBA);
Gamma1s_ore_2trans = R21afterBA'*(GammaTangent_3D_2);

I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,ref_idx);
abs_C2afterBA = -R_matrix(:,:,ref_idx)' *...
                 T_matrix(:,ref_idx);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_orin_2 = -R21afterBA'*(Gamma1s_2trans-T21afterBA);
Gamma1s_ore_orin1_2 = -R21afterBA'*(Gamma1s_ore_2trans);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID3 = fopen('GenData/pairededge_ABC0006_42n47_t32to0excludehypo1n2_delta03_theta15_N4.txt','r');
% fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
formatSpec = '%f';
A3 = fscanf(fileID3,formatSpec);
B3 = reshape(A3,50,size(A3,1)/50);
B3=B3+1;
paired_edge3 = B3';

params.HYPO1_VIEW_INDX          = 43;
params.HYPO2_VIEW_INDX          = 48;

if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
else
    invK = inv(K);
end

%> Fetch edgels of the two hypothesis views
% load('F:\T-less\edges50_1.mat');
TO_Edges_HYPO1  = collect_third_order_edges{params.HYPO1_VIEW_INDX,6};
TO_Edges_HYPO2  = collect_third_order_edges{params.HYPO2_VIEW_INDX,6};
Gamma1s_3         = [];
GammaTangent_3D_3 = [];

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

for(i = 1 : size(paired_edge3,1))
    if((size(find(paired_edge1(i,:)>0),2)-2)<Nview)
        continue;
    end
    edgel_HYPO1 = TO_Edges_HYPO1(paired_edge3(i,1), :);
    edgel_HYPO2 = TO_Edges_HYPO2(paired_edge3(i,2), :);
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    [edgels_HYPO2_corrected, ...
     edgels_HYPO1_corrected, e_i] = edgelsHYPO2correct(Epipolar_Coeffs, ...
                                                       edgel_HYPO2, ...
                                                       edgel_HYPO1, ...
                                                       R_matrix, ...
                                                       T_matrix, ...
                                                       params, K);
    edgels_HYPO2_final = edgels_HYPO2_corrected;
    edgels_HYPO1_final = edgels_HYPO1_corrected;
    if(params.multiK == 1)
                    pt1 = invK1 * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK2 * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                else
                    pt1 = invK * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                end
    pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
    Gamma = linearTriangulation(2, pts_meters, R21, T21);
    GammaTgt_3D     = get3DOrientation(edgels_HYPO1_final,edgels_HYPO2_final, K, params, R21);
    Gamma1s_3 = [Gamma1s_3, Gamma];
    GammaTangent_3D_3 = [GammaTangent_3D_3, GammaTgt_3D];
end


abs_R1afterBA =  R_matrix(:,:,ref_idx);
abs_C1afterBA = -R_matrix(:,:,ref_idx)' *...
    T_matrix(:,ref_idx);
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_3trans = R21afterBA'*(Gamma1s_3-T21afterBA);
Gamma1s_ore_3trans = R21afterBA'*(GammaTangent_3D_3);

I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,ref_idx);
abs_C2afterBA = -R_matrix(:,:,ref_idx)' *...
                 T_matrix(:,ref_idx);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_orin_3 = -R21afterBA'*(Gamma1s_3trans-T21afterBA);
Gamma1s_ore_orin1_3 = -R21afterBA'*(Gamma1s_ore_3trans);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID4 = fopen('GenData/pairededge_ABC0006_27n46_t32to0excludehypo1n2_delta03_theta15_N4.txt','r');
% fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
formatSpec = '%f';
A4 = fscanf(fileID4,formatSpec);
B4 = reshape(A4,50,size(A4,1)/50);
B4=B4+1;
paired_edge4 = B4';

params.HYPO1_VIEW_INDX          = 28;
params.HYPO2_VIEW_INDX          = 47;

if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
else
    invK = inv(K);
end

%> Fetch edgels of the two hypothesis views
% load('F:\T-less\edges50_1.mat');
TO_Edges_HYPO1  = collect_third_order_edges{params.HYPO1_VIEW_INDX,6};
TO_Edges_HYPO2  = collect_third_order_edges{params.HYPO2_VIEW_INDX,6};
Gamma1s_4         = [];
GammaTangent_3D_4 = [];

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

for(i = 1 : size(paired_edge4,1))
    if((size(find(paired_edge1(i,:)>0),2)-2)<Nview)
        continue;
    end
    edgel_HYPO1 = TO_Edges_HYPO1(paired_edge4(i,1), :);
    edgel_HYPO2 = TO_Edges_HYPO2(paired_edge4(i,2), :);
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    [edgels_HYPO2_corrected, ...
     edgels_HYPO1_corrected, e_i] = edgelsHYPO2correct(Epipolar_Coeffs, ...
                                                       edgel_HYPO2, ...
                                                       edgel_HYPO1, ...
                                                       R_matrix, ...
                                                       T_matrix, ...
                                                       params, K);
    edgels_HYPO2_final = edgels_HYPO2_corrected;
    edgels_HYPO1_final = edgels_HYPO1_corrected;
    if(params.multiK == 1)
                    pt1 = invK1 * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK2 * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                else
                    pt1 = invK * [edgels_HYPO1_final(1, 1:2)'; 1]; pt1 = pt1(1:2, 1);
                    pt2 = invK * [edgels_HYPO2_final(1, 1:2)'; 1]; pt2 = pt2(1:2, 1);
                end
    pts_meters = [pt1(1), pt2(1); pt1(2), pt2(2)];
    Gamma = linearTriangulation(2, pts_meters, R21, T21);
    GammaTgt_3D     = get3DOrientation(edgels_HYPO1_final,edgels_HYPO2_final, K, params, R21);
    Gamma1s_4 = [Gamma1s_4, Gamma];
    GammaTangent_3D_4 = [GammaTangent_3D_4, GammaTgt_3D];
end


abs_R1afterBA =  R_matrix(:,:,ref_idx);
abs_C1afterBA = -R_matrix(:,:,ref_idx)' *...
    T_matrix(:,ref_idx);
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_4trans = R21afterBA'*(Gamma1s_4-T21afterBA);
Gamma1s_ore_4trans = R21afterBA'*(GammaTangent_3D_4);

I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,ref_idx);
abs_C2afterBA = -R_matrix(:,:,ref_idx)' *...
                 T_matrix(:,ref_idx);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_orin_4 = -R21afterBA'*(Gamma1s_4trans-T21afterBA);
Gamma1s_ore_orin1_4 = -R21afterBA'*(Gamma1s_ore_4trans);
%}

Gamma1_all = [Gamma1s_orin_1,Gamma1s_orin_2,Gamma1s_orin_3, Gamma1s_orin_4];
figure
plot3(Gamma1_all(1,:), Gamma1_all(2,:), Gamma1_all(3,:), 'b.');
namefile = ['ABC/ABC0006_delta',deltastr,'theta',num2str(thetaangle),'N',num2str(Nview),'.mat'];
save(namefile, "Gamma1_all");