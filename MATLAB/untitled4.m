load('MATData/TO_Edges_10.mat');
load('MATData/imageArray_10.mat')
R_matrixbeforeBA = load('MATData/R_matrix_10.mat').R_matrix;
T_matrixbeforeBA = load('MATData/T_matrix_10.mat').T_matrix;
R_matrixafterBA  = load('MATData/R_matrix.mat').R_matrix;
T_matrixafterBA  = load('MATData/T_matrix.mat').T_matrix;
load('MATData/K_10.mat')
% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 6;
params.HYPO2_VIEW_INDX          = 16;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.99;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.005*5;
params.ANGLE_FOR_EPIPOLE        = 0.03;
params.ICL_DATA                 = 0;
params.multiK                   = 1;
params.cols                     = size(double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX})), 2);
%> Convert from orientation to tangent vector for each edgel
len = 2;
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

load('collective\correcthnh03_BA.mat')
recons_coor2 = Gamma1s';
load('collective\correcthnh03_tgt_BA.mat')
GammaTangent_3D2 = GammaTangent_3D;
fignum = 1;
subnum = 1;
e3 = [0;0;1];
abs_R1afterBA  =  R_matrixafterBA(:,:,params.HYPO1_VIEW_INDX);
abs_C1afterBA  = -R_matrixafterBA(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrixafterBA(:,params.HYPO1_VIEW_INDX);

reproject_error1_all = [];
reproject_error2_all = [];
validateidex = [6,16,7];
rangeplot = [203   204   208   214   215   216   221   222   229   230   238   239   247   255];
reproject_coor2  = [];
reproject_tgt2   = [];
reproject_error2 = [];

for(vidx = 1:3)
    vi = validateidex(1,vidx);
    VALID_INDX = vi;%validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    Kmatrix    = K(:,:,VALID_INDX);
    invK       = inv(Kmatrix);
    abs_R2afterBA =  R_matrixafterBA(:,:,VALID_INDX);
    abs_C2afterBA = -R_matrixafterBA(:,:,VALID_INDX)' *...
        T_matrixafterBA(:,VALID_INDX);
    R21afterBA    = abs_R2afterBA * abs_R1afterBA';
    T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);

    extrafterBA           = zeros(3,4);
    extrafterBA(1:3, 1:3) = R21afterBA;
    extrafterBA(1:3, 4)   = T21afterBA;
    PmatrixafterBA        = Kmatrix*extrafterBA;


    i = rangeplot(1,1);
    reprojection = PmatrixafterBA*[recons_coor2(i,1);recons_coor2(i,2);recons_coor2(i,3);1];
    x_coor       = reprojection(1,1)/reprojection(3,1);
    y_coor       = reprojection(2,1)/reprojection(3,1);
    reproject_coor2 = [reproject_coor2; [x_coor, y_coor]];

    T_v3 = R21afterBA * GammaTangent_3D2(:,i);
    t_v3 = T_v3 - (e3'*T_v3)*invK*[reproject_coor2(1,:),1]';
    t_v3_reproj = t_v3 ./ norm(t_v3);
    reproject_tgt2  = [reproject_tgt2; t_v3_reproj(1,1), t_v3_reproj(2,1)];

    %{}
    figure(1);
    subplot(1,3,subnum)
    imshow(imageArray{VALID_INDX});
    hold on;
    plot(TO_Edges_VALID(:,1), TO_Edges_VALID(:,2), 'y.', 'MarkerSize',2, 'LineWidth', 1);
    plot(reproject_coor2(vidx,1), reproject_coor2(vidx,2), 'gx', 'MarkerSize',10, 'LineWidth', 1);
    edgel_tgt2  = reproject_tgt2(vidx, 1:2);
    hold on;
    x1 = reproject_coor2(vidx, 1)+len*edgel_tgt2(1,1);
    x2 = reproject_coor2(vidx, 1)-len*edgel_tgt2(1,1);
    y1 = reproject_coor2(vidx, 2)+len*edgel_tgt2(1,2);
    y2 = reproject_coor2(vidx, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.25);
    title(['view ', num2str(vi)])
    subnum = subnum+1;
end
subnum = 1;
refidx = [3, 7, 10];
newidx = [7, 10, 3; 10, 3, 7];
subidx = [2, 3, 1; 3, 1, 2]; 
colorlin3 = [0 1 1;
             1 0 1;
             0.4660 0.6740 0.1880];
R_relative = [];
T_relative = [];
for(vidx = 1:3)
    idx_ref = refidx(1,vidx);
    idx_new = newidx(1,vidx);
    abs_R1 =  R_matrixafterBA(:,:,idx_ref);
    abs_C1 = -R_matrixafterBA(:,:,idx_ref)' *...
        T_matrixafterBA(:,idx_ref);
    abs_R2 =  R_matrixafterBA(:,:,idx_new);
    abs_C2 = -R_matrixafterBA(:,:,idx_new)' *...
        T_matrixafterBA(:,idx_new);

    % CPP modification
    R21    = abs_R2 * abs_R1';
    % CPP modification
    T21    = abs_R2 * (abs_C1 - abs_C2);
    R_relative{vidx*2-1} = R21;
    T_relative{vidx*2-1} = T21;
    T_x = @(T)[0,      -T(3,1),  T(2,1); ...
        T(3,1),  0,      -T(1,1); ...
        -T(2,1),  T(1,1),  0];
    E   = T_x(T21) * R21;

    % Calculate fundamental matrix
    if(params.multiK == 1)
        K1 = K(:,:,idx_ref);
        K2 = K(:,:,idx_new);
        invK1 = inv(K1);
        invK2 = inv(K2);
        F   = invK2'*E*invK1;
    else
        invK = inv(K);
        F   = invK'*E*invK;
    end
    coeffs = F * [reproject_coor2(vidx,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*params.cols) ./ Epipolar_Coeffs.B;
    subplot(1,3,subidx(1,vidx))
    hold on;
    line([0, params.cols], [yMin, yMax], 'Color', colorlin3(vidx,:), 'LineWidth', 0.5);

    %%%%%
    idx_new = newidx(2,vidx);
    abs_R1 =  R_matrixafterBA(:,:,idx_ref);
    abs_C1 = -R_matrixafterBA(:,:,idx_ref)' *...
        T_matrixafterBA(:,idx_ref);
    abs_R2 =  R_matrixafterBA(:,:,idx_new);
    abs_C2 = -R_matrixafterBA(:,:,idx_new)' *...
        T_matrixafterBA(:,idx_new);

    % CPP modification
    R21    = abs_R2 * abs_R1';
    % CPP modification
    T21    = abs_R2 * (abs_C1 - abs_C2);
    R_relative{vidx*2} = R21;
    T_relative{vidx*2} = T21;
    T_x = @(T)[0,      -T(3,1),  T(2,1); ...
        T(3,1),  0,      -T(1,1); ...
        -T(2,1),  T(1,1),  0];
    E   = T_x(T21) * R21;

    % Calculate fundamental matrix
    if(params.multiK == 1)
        K1 = K(:,:,idx_ref);
        K2 = K(:,:,idx_new);
        invK1 = inv(K1);
        invK2 = inv(K2);
        F   = invK2'*E*invK1;
    else
        invK = inv(K);
        F   = invK'*E*invK;
    end
    coeffs = F * [reproject_coor2(vidx,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*params.cols) ./ Epipolar_Coeffs.B;
    subplot(1,3,subidx(2,vidx))
    hold on;
    line([0, params.cols], [yMin, yMax], 'Color', colorlin3(vidx,:), 'LineWidth', 0.5);
end

R_rel_pret = R_relative;
T_rel_pret = T_relative;

% % % for(perturbidx = 1:6)
% % %     eulangleZYX = rotm2eul(R_relative{perturbidx});
% % %     R_rel_pret{perturbidx} = eul2rotm((0.99 + (0.01)*rand(1))*eulangleZYX);
% % %     T_rel_pret{perturbidx} = (0.99 + (0.01)*rand(1))*T_relative{perturbidx};
% % % end

% % % for(vidx = 1:3)
% % %     figure(1);
% % %     subplot(1,3,subnum)
% % % %     imshow(imageArray{VALID_INDX});
% % %     hold on;
% % %     plot(reproject_coor2(vidx,1), reproject_coor2(vidx,2), 'gx', 'MarkerSize',10, 'LineWidth', 1);
% % % %     title(['view ', num2str(vi)])
% % %     subnum = subnum+1;
% % % end
reproject_coor2_pret  = [];
reproject_tgt2_pret   = [];
for(vidx = 1:3)
    vi = validateidex(1,vidx);
    VALID_INDX = vi;%validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    Kmatrix    = K(:,:,VALID_INDX);
    invK       = inv(Kmatrix);
    abs_R2afterBA =  R_matrixafterBA(:,:,VALID_INDX);
    abs_C2afterBA = -R_matrixafterBA(:,:,VALID_INDX)' *...
        T_matrixafterBA(:,VALID_INDX);
    R21afterBA    = abs_R2afterBA * abs_R1afterBA';
    T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);

    % perturb here
    eulangleZYX = rotm2eul(R21afterBA);
    R21afterBA_pret = eul2rotm((0.99 + (0.01)*rand(1))*eulangleZYX);
    T21afterBA_pret = (0.99 + (0.01)*rand(1))*T21afterBA;

    extrafterBA           = zeros(3,4);
    extrafterBA(1:3, 1:3) = R21afterBA_pret;
    extrafterBA(1:3, 4)   = T21afterBA_pret;
    PmatrixafterBA        = Kmatrix*extrafterBA;


    i = rangeplot(1,1);
    reprojection = PmatrixafterBA*[recons_coor2(i,1);recons_coor2(i,2);recons_coor2(i,3);1];
    x_coor       = reprojection(1,1)/reprojection(3,1);
    y_coor       = reprojection(2,1)/reprojection(3,1);
    reproject_coor2_pret = [reproject_coor2_pret; [x_coor, y_coor]];

    T_v3 = R21afterBA_pret * GammaTangent_3D2(:,i);
    t_v3 = T_v3 - (e3'*T_v3)*invK*[reproject_coor2(1,:),1]';
    t_v3_reproj = t_v3 ./ norm(t_v3);
    reproject_tgt2_pret  = [reproject_tgt2_pret; t_v3_reproj(1,1), t_v3_reproj(2,1)];

    %{}
    figure(2);
    subplot(1,3,subnum)
    imshow(imageArray{VALID_INDX});
    hold on;
    plot(TO_Edges_VALID(:,1), TO_Edges_VALID(:,2), 'y.', 'MarkerSize',2, 'LineWidth', 1);
    plot(reproject_coor2_pret(vidx,1), reproject_coor2_pret(vidx,2), 'gx', 'MarkerSize',10, 'LineWidth', 1);
    edgel_tgt2  = reproject_tgt2_pret(vidx, 1:2);
    hold on;
    x1 = reproject_coor2_pret(vidx, 1)+len*edgel_tgt2(1,1);
    x2 = reproject_coor2_pret(vidx, 1)-len*edgel_tgt2(1,1);
    y1 = reproject_coor2_pret(vidx, 2)+len*edgel_tgt2(1,2);
    y2 = reproject_coor2_pret(vidx, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.25);
    title(['view ', num2str(vi)])
    subnum = subnum+1;
end
subnum = 1;
refidx = [3, 7, 10];
newidx = [7, 10, 3; 10, 3, 7];
subidx = [2, 3, 1; 3, 1, 2]; 
colorlin3 = [0 1 1;
             1 0 1;
             0.4660 0.6740 0.1880];
R_relative = [];
T_relative = [];
for(vidx = 1:3)
    idx_ref = refidx(1,vidx);
    idx_new = newidx(1,vidx);
    abs_R1 =  R_matrixafterBA(:,:,idx_ref);
    abs_C1 = -R_matrixafterBA(:,:,idx_ref)' *...
        T_matrixafterBA(:,idx_ref);
    abs_R2 =  R_matrixafterBA(:,:,idx_new);
    abs_C2 = -R_matrixafterBA(:,:,idx_new)' *...
        T_matrixafterBA(:,idx_new);

    % CPP modification
    R21    = abs_R2 * abs_R1';
    % CPP modification
    T21    = abs_R2 * (abs_C1 - abs_C2);
    %%%
    eulangleZYX = rotm2eul(R21);
    R21 = eul2rotm((0.99 + (0.01)*rand(1))*eulangleZYX);
    T21 = (0.99 + (0.01)*rand(1))*T21;
%     R_relative{vidx*2-1} = R21;
%     T_relative{vidx*2-1} = T21;
    T_x = @(T)[0,      -T(3,1),  T(2,1); ...
        T(3,1),  0,      -T(1,1); ...
        -T(2,1),  T(1,1),  0];
    E   = T_x(T21) * R21;

    % Calculate fundamental matrix
    if(params.multiK == 1)
        K1 = K(:,:,idx_ref);
        K2 = K(:,:,idx_new);
        invK1 = inv(K1);
        invK2 = inv(K2);
        F   = invK2'*E*invK1;
    else
        invK = inv(K);
        F   = invK'*E*invK;
    end
    coeffs = F * [reproject_coor2_pret(vidx,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*params.cols) ./ Epipolar_Coeffs.B;
    subplot(1,3,subidx(1,vidx))
    hold on;
    line([0, params.cols], [yMin, yMax], 'Color', colorlin3(vidx,:), 'LineWidth', 0.5);

    %%%%%
    idx_new = newidx(2,vidx);
    abs_R1 =  R_matrixafterBA(:,:,idx_ref);
    abs_C1 = -R_matrixafterBA(:,:,idx_ref)' *...
        T_matrixafterBA(:,idx_ref);
    abs_R2 =  R_matrixafterBA(:,:,idx_new);
    abs_C2 = -R_matrixafterBA(:,:,idx_new)' *...
        T_matrixafterBA(:,idx_new);

    % CPP modification
    R21    = abs_R2 * abs_R1';
    % CPP modification
    T21    = abs_R2 * (abs_C1 - abs_C2);
    %%%
    eulangleZYX = rotm2eul(R21);
    R21 = eul2rotm((0.99 + (0.01)*rand(1))*eulangleZYX);
    T21 = (0.99 + (0.01)*rand(1))*T21;
%     R_relative{vidx*2} = R21;
%     T_relative{vidx*2} = T21;
    T_x = @(T)[0,      -T(3,1),  T(2,1); ...
        T(3,1),  0,      -T(1,1); ...
        -T(2,1),  T(1,1),  0];
    E   = T_x(T21) * R21;

    % Calculate fundamental matrix
    if(params.multiK == 1)
        K1 = K(:,:,idx_ref);
        K2 = K(:,:,idx_new);
        invK1 = inv(K1);
        invK2 = inv(K2);
        F   = invK2'*E*invK1;
    else
        invK = inv(K);
        F   = invK'*E*invK;
    end
    coeffs = F * [reproject_coor2_pret(vidx,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
    yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*params.cols) ./ Epipolar_Coeffs.B;
    subplot(1,3,subidx(2,vidx))
    hold on;
    line([0, params.cols], [yMin, yMax], 'Color', colorlin3(vidx,:), 'LineWidth', 0.5);
end