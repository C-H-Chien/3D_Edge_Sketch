fileIDtest = fopen('GenData/Gamma1s_ABC0006_46n20_t32to0excludehypo1n2_delta03_theta15_N4.txt','r');
% fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
formatSpec = '%f';
Atest = fscanf(fileIDtest,formatSpec);
Btest = reshape(Atest,3,size(Atest,1)/3);
Gamma1test = Btest;

load('MATData/TO_Edges_ABC00000006.mat');
load('MATData/imageArray_ABC00000006.mat')
load('MATData/R_matrix_ABC00000006.mat')
load('MATData/T_matrix_ABC00000006.mat')
load('MATData/K_ABC00000006.mat')

% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 47;
params.HYPO2_VIEW_INDX          = 21;
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
params.delta                    = 0.3;
params.circle                   = 55;
params.parallelangle            = 15;
% params.SUPPORT_OREN_THRESH_DEG  = 15;

hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));

if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
else
    invK = inv(K);
end

%> Convert from orientation to tangent vector for each edgel
TO_Edges_Tangents = cell(params.NUM_OF_IMGS, 1);
if(params.syn == 0)
    for i = 1:params.NUM_OF_IMGS
        TO_Edgels_Orientations = collect_third_order_edges{i,6}(:,3);

        Tgts = zeros(size(TO_Edgels_Orientations, 1), 2);
        for j = 1:size(TO_Edgels_Orientations, 1)
            Tgts(j,:) = [cos(TO_Edgels_Orientations(j,1)), ...
                sin(TO_Edgels_Orientations(j,1))];
        end
        TO_Edges_Tangents{i,1} = Tgts;
    end

else
    for i = 1:params.NUM_OF_IMGS
        TO_Edges_Tangents{i,1} = collect_third_order_edges{i,6}(:,3:4);
    end
end

% tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
% tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};

%> Array of view indices
view_idx = [];
% figure;
for i = 1:params.NUM_OF_IMGS; view_idx = [view_idx; i]; end
validation_view_indices = view_idx;
validation_view_indices([params.HYPO1_VIEW_INDX; ...
                         params.HYPO2_VIEW_INDX],:) = [];

%> Fetch edgels of the two hypothesis views
% load('F:\T-less\edges50_1.mat');
TO_Edges_HYPO1  = collect_third_order_edges{params.HYPO1_VIEW_INDX,6};
TO_Edges_HYPO2  = collect_third_order_edges{params.HYPO2_VIEW_INDX,6};
Gamma1s         = [];
GammaTangent_3D = [];

hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
rows      = size(hypo_img2, 1);
cols      = size(hypo_img2, 2);

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

finalEdgePair    = [];
% paired_edge      = [];

% get the orientation list between hypo1 and hypo2
[ore_list1bar_all, ...
 ore_list2_sorted, ...
 epipole_pix_view1, ...
 epipole_pix_view2]= getOreList_New(TO_Edges_HYPO1, TO_Edges_HYPO2, ...
                                     R_matrix, T_matrix, ...
                                     params.HYPO1_VIEW_INDX, ...
                                     params.HYPO2_VIEW_INDX, params, K);

I  = [1,1,1];
abs_R1afterBA =  diag(I);
abs_C1afterBA = -diag(I)' *...
                 [0; 0; 0];
abs_R2afterBA =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2afterBA = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
                 T_matrix(:,params.HYPO1_VIEW_INDX);
R21afterBA    = abs_R2afterBA * abs_R1afterBA';
T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
Gamma1s_orintest = -R21afterBA'*(Gamma1test-T21afterBA);
% Gamma1s_ore_orin1 = -R21afterBA'*(GammaTangent_3D);

figure;
plot3(Gamma1s_orintest(1,:)*100, Gamma1s_orintest(2,:)*100, Gamma1s_orintest(3,:)*100, 'b.');
axis equal;
