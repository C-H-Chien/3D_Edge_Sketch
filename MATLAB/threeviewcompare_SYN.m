% close all
clear all
rng(0)
load('MATData/TO_Edges_syn.mat');
load('MATData/imageArray_syn.mat')
R_g  = load('MATData/R_matrix_syn.mat').R_matrix;
T_g  = load('MATData/T_matrix_syn.mat').T_matrix;
load('MATData/K_syn.mat')

% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 3;
params.HYPO2_VIEW_INDX          = 5;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.94;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.005*5;
params.ANGLE_FOR_EPIPOLE        = 0.001;
params.ICL_DATA                 = 0;
params.multiK                   = 0;
params.syn                      = 1;
params.cols                     = size(double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX})), 2);
params.rows                     = size(double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX})), 1);
params.delta                    = 0.3;
%> Convert from orientation to tangent vector for each edgel
len = 2;
cols = params.cols;
TO_Edges_Tangents = cell(params.NUM_OF_IMGS, 1);
for i = 1:params.NUM_OF_IMGS
    TO_Edges_Tangents{i,1} = collect_third_order_edges{i,1}(:,3:4);
end

load('MATData\TO_Edges3D_syn.mat')
recons_coor2 = points3D;
load('MATData\TO_Edges3D_syn_tgt.mat')
GammaTangent_3D2 = tangent3D';
% fignum = 1;
subnum = 1;
e3 = [0;0;1];
abs_R1afterBA  =  R_g(:,:,params.HYPO1_VIEW_INDX);
abs_C1afterBA  = -R_g(:,:,params.HYPO1_VIEW_INDX)' *...
                  T_g(:,params.HYPO1_VIEW_INDX);

reproject_error1_all = [];
reproject_error2_all = [];

validateidex = [params.HYPO1_VIEW_INDX,params.HYPO2_VIEW_INDX,10];

VALID_INDX = 10;
rangeplot = [297   204   208   214   215   216   221   222   229   230   238   239   247   255];
reproject_coor2  = [];
reproject_tgt2   = [];
reproject_error2 = [];

% rangeplot = [203]%   204   208   214   215   216   221   222   229   230   238   239   247   255];
% % figure(5);
% % plot3(recons_coor2(:,1), -recons_coor2(:,2), -recons_coor2(:,3), 'b.', 'MarkerSize', 3, 'LineWidth', 1);
% % hold on;
% % plot3(recons_coor2(rangeplot(1,1),1), -recons_coor2(rangeplot(1,1),2), -recons_coor2(rangeplot(1,1),3), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
% % axis equal;
% % title (['Synthetic 3D edge with chosen 3D edge plotted in green'])
% % i = rangeplot(1,1);
% % for(vidx = 1:3)
% %     vi = validateidex(1,vidx);
% %     VALID_INDX = vi;%validation_view_indices(vi, 1);
% %     VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    TO_Edges_VALI = collect_third_order_edges{VALID_INDX,1};



    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    Kmatrix    = K;
    invK       = inv(Kmatrix);
close all;
TO_Edges_HYPO1          = collect_third_order_edges{params.HYPO1_VIEW_INDX,1};
TO_Edges_HYPO2          = collect_third_order_edges{params.HYPO2_VIEW_INDX,1};
TO_Edges_VALID          = collect_third_order_edges{VALID_INDX,1};
% [ore_list1bar, ore_list2_sorted, epipole_pix_view1, epipole_pix_view2]= ...
%     getOreList(TO_Edges_HYPO1, TO_Edges_HYPO2, R_g, T_g, params, K);
[ore_list1bar_all, ore_list2_sorted1, epipole_pix_view_1, epipole_pix_view_2]= ...
    getOreList_New(TO_Edges_HYPO1, TO_Edges_HYPO2, R_g, T_g, params.HYPO1_VIEW_INDX, params.HYPO2_VIEW_INDX, params, K);
R1                      = R_g(:,:,params.HYPO1_VIEW_INDX);
T1                      = T_g(:,:,params.HYPO1_VIEW_INDX);
R2                      = R_g(:,:,params.HYPO2_VIEW_INDX);
T2                      = T_g(:,:,params.HYPO2_VIEW_INDX);
% pix_perturb             = 0.1;
% params.ANGLE_FOR_EPIPOLE= 0.006;
[R21, T21, E, F]        = getRelativePose(R_g, T_g, params, K);
pix_pert        = [];
percentage_vali = [];
TO_Edges_VALID = TO_Edges_VALI;
for(pix_perturb = 0.1:0.05:params.delta*2)

%     TO_Edges_VALID = zeros(size(TO_Edges_VALI));
% for(valid_i = 1 : size(TO_Edges_VALI,1))
%     angle_perv = rand(1)*2*pi;
%     pertxv     = pix_perturb*cos(angle_perv);
%     pertyv     = pix_perturb*sin(angle_perv);
%     TO_Edges_VALID(valid_i, :) = [TO_Edges_VALI(valid_i, 1)+(0.75+0.25*rand(1))*pertxv, ...
%         TO_Edges_VALI(valid_i,2)+(0.75+0.25*rand(1))*pertyv, ...
%         TO_Edges_VALI(valid_i,3:4)];
% end

    supported_indices_stack = [];
    supported_link          = [];
    supported_simi          = [];
    reprojerror             = [];
    realvalidinside         = [];
    realvalidoutside        = [];
    emptyquad               = [];
for(idx=3087:size(recons_coor2,1))
    reproject_coor2 = [];
    reproject_tgt2  = [];
    reprojection = Kmatrix * (R1*[recons_coor2(idx,1);recons_coor2(idx,2);recons_coor2(idx,3)] + T1);
    x_coor       = reprojection(1,1)/reprojection(3,1);
    y_coor       = reprojection(2,1)/reprojection(3,1);
    reproject_coor2 = [reproject_coor2; [x_coor, y_coor]];
    T_v3 = R1 * GammaTangent_3D2(:,idx);
    t_v3 = T_v3 - (e3'*T_v3)*invK*[reproject_coor2(1,:),1]';
    t_v3_reproj = t_v3 ./ norm(t_v3);
    reproject_tgt2  = [reproject_tgt2; t_v3_reproj(1,1), t_v3_reproj(2,1)];
    reprojection = Kmatrix * (R2*[recons_coor2(idx,1);recons_coor2(idx,2);recons_coor2(idx,3)] + T2);
    x_coor       = reprojection(1,1)/reprojection(3,1);
    y_coor       = reprojection(2,1)/reprojection(3,1);
    reproject_coor2 = [reproject_coor2; [x_coor, y_coor]];
    T_v3 = R2 * GammaTangent_3D2(:,idx);
    t_v3 = T_v3 - (e3'*T_v3)*invK*[reproject_coor2(2,:),1]';
    t_v3_reproj = t_v3 ./ norm(t_v3);
    reproject_tgt2  = [reproject_tgt2; t_v3_reproj(1,1), t_v3_reproj(2,1)];

    %{}    
        angle_per1 = rand(1)*2*pi;
        angle_per2 = rand(1)*2*pi;
        pertx1     = pix_perturb*cos(angle_per1);
        pertx2     = pix_perturb*cos(angle_per2);
        perty1     = pix_perturb*sin(angle_per1);
        perty2     = pix_perturb*sin(angle_per2);
        reproject_coor2_pret = [reproject_coor2(1,1)+(0.75+0.25*rand(1))*pertx1, ...
                                reproject_coor2(1,2)+(0.75+0.25*rand(1))*perty1; ...
                                reproject_coor2(2,1)+(0.75+0.25*rand(1))*pertx2, ...
                                reproject_coor2(2,2)+(0.75+0.25*rand(1))*perty2];
        
%         edgel_HYPO1        = [reproject_coor2(1,:),reproject_tgt2(1,:)];
%         edgels_HYPO2       = [reproject_coor2(2,:),reproject_tgt2(2,:)];
        edgel_HYPO1        = [reproject_coor2_pret(1,:),reproject_tgt2(1,:)];
        edgels_HYPO2       = [reproject_coor2_pret(2,:),reproject_tgt2(2,:)];

        [ore_list1bar_cur, ore_list2_cur, ~, ~]= ...
    getOreList_New(edgel_HYPO1, edgels_HYPO2, R_g, T_g, params.HYPO1_VIEW_INDX, params.HYPO2_VIEW_INDX, params, K);
        if(ore_list2_cur(1,1) > max(ore_list1bar_cur) || ore_list2_cur(1,1) < min(ore_list1bar_cur) )
            continue;
        end

        coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
        Apixel = coeffs(1,1);
        Bpixel = coeffs(2,1);
        Cpixel = coeffs(3,1);
        Epipolar_Coeffs.A = Apixel;
        Epipolar_Coeffs.B = Bpixel;
        Epipolar_Coeffs.C = Cpixel;

        [edgels_HYPO2_corrected, edgels_HYPO1_corrected] = edgelsHYPO2correct(Epipolar_Coeffs, edgels_HYPO2, edgel_HYPO1, R_g, T_g, params, K);

        edgels_HYPO1_final = edgels_HYPO1_corrected;
        edgels_HYPO2_final = edgels_HYPO2_corrected;
%                 edgels_HYPO1_final = edgel_HYPO1;
%                 edgels_HYPO2_final = edgels_HYPO2;

        yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
        yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
        epipolarslope = (yMin-yMax)/(0-cols);
        edgeslope     = edgels_HYPO2_final(1,4)/edgels_HYPO2_final(1,3);
        pert1dist = norm(edgels_HYPO1_corrected(:,1:2) - reproject_coor2(1,:));
        pert2dist = norm(edgels_HYPO2_corrected(:,1:2) - reproject_coor2(2,:));
%         pert1dist = norm(edgels_HYPO1_corrected(:,1:2) - reproject_coor2_pret(1,:));
%         pert2dist = norm(edgels_HYPO2_corrected(:,1:2) - reproject_coor2_pret(2,:));
        
%         if(pert1dist > params.delta || pert2dist > params.delta )
%             continue;
%         end
        if(pert1dist > pix_perturb || pert2dist > pix_perturb )
            continue;
        end
        h2_angle = atan((edgels_HYPO2_corrected(1,2) - epipole_pix_view_2(2,1))/...
                        (edgels_HYPO2_corrected(1,1) - epipole_pix_view_2(1,1)));
        h2_angle = rad2deg(h2_angle);
        if(h2_angle > max(ore_list1bar_cur) || h2_angle < min(ore_list1bar_cur) )
            continue;
        end
%         if(pert1dist > params.delta || pert2dist > params.delta )
%             continue;
%         end
%         if(abs(atan(epipolarslope)-atan(edgeslope))<=0.1745)
%             continue;
%         end


        % get the validation view
        VALID_INDX = validateidex(1,3);
        VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));        
        % get edges in the validation view
        TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
        Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
        % get the orientation list between hypo1 and validation view
% %         [ore_31bar, ore_list31_sorted, epipole_pix_view1_3, epipole_pix_view3_1] = ...
% %             getOreList_Vali(edgels_HYPO1_final, TO_Edges_VALID, R_g, T_g, ...
% %             params.HYPO1_VIEW_INDX, K, VALID_INDX, params);

        [ore_31bar_all, ore_list31_sorted, epipole_pix_view1_3, epipole_pix_view3_1]= getOreList_New(edgels_HYPO1_final, TO_Edges_VALID, R_g, T_g, params.HYPO1_VIEW_INDX, VALID_INDX, params, K);
%         % find epipolar angle range for validation view
%         angle_range2 = abs(ore_list31_sorted(1,1)- ...
%             ore_list31_sorted(size(ore_list31_sorted,1),1));
%         % get the range for single epipolar wedge in vali
%         range2       = params.ANGLE_FOR_EPIPOLE; %PERCENTAGE_FOR_EPIPOLE*angle_range2;
        %         [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgel ...
%              (edgel_HYPO1, edgels_HYPO2_final, R_matrix, T_matrix, params, VALID_INDX, K);
reproj_edge_pos_gamma3 = [];
reproj_edge_tgt_gamma3 = [];
        [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgelCorrected ...
             (edgels_HYPO1_final, edgels_HYPO2_final, R_g, T_g, params, VALID_INDX, K);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the orientation list between hypo2 and validation view
% %         [ore_32bar, ore_list32_sorted, epipole_pix_view2_3, epipole_pix_view3_2] = ...
% %             getOreList_Vali(edgels_HYPO2_final, TO_Edges_VALID, R_g, T_g, ...
% %             params.HYPO2_VIEW_INDX, K, VALID_INDX, params);
        [ore_32bar_all, ore_list32_sorted, epipole_pix_view2_3, epipole_pix_view3_2]= getOreList_New(edgels_HYPO2_final, TO_Edges_VALID, R_g, T_g, params.HYPO2_VIEW_INDX, VALID_INDX, params, K);
        % find epipolar angle range for validation view
% %         angle_range3 = abs(ore_list32_sorted(1,1) - ...
% %             ore_list32_sorted(size(ore_list32_sorted,1),1));
% %         % get the range for single epipolar wedge in vali
% %         range3       = params.ANGLE_FOR_EPIPOLE; %PERCENTAGE_FOR_EPIPOLE*angle_range3;
        % investigate all the potential pairs
        % "form" the quadrilateral
        allquad = [];
        isparal = ones(size(edgels_HYPO2_final,1),1);
        hypo2idx = 1;
            % get the epipolar wedge range for the edge in hypo1 on vali
            thresh_ore31_1 = min(ore_31bar_all(hypo2idx,:));
            thresh_ore31_2 = max(ore_31bar_all(hypo2idx,:));
            idxlist_31     = find(ore_list31_sorted(:,1) >= thresh_ore31_1 & ...
                ore_list31_sorted(:,1) <= thresh_ore31_2);
            vali_idx     = ore_list31_sorted(idxlist_31,2);
            edgels_31    = TO_Edges_VALID(vali_idx, :);
            % get the epipolar wedge range for the edge in hypo2 on vali
            thresh_ore32_1 = min(ore_32bar_all(hypo2idx,:));
            thresh_ore32_2 = max(ore_32bar_all(hypo2idx,:));
            % get edges in vali fall inside the epipolar wedge
            idxlist_32 = find(ore_list32_sorted(:,1) >= thresh_ore32_1 & ...
                ore_list32_sorted(:,1) <= thresh_ore32_2);
            vali_idx1 = ore_list32_sorted(idxlist_32,2);
            edgels_32 = TO_Edges_VALID(vali_idx1, :);
            % find the edges falls inside both two wedges
            quad_idx  = intersect(vali_idx, vali_idx1);
            if(isempty(find(quad_idx == idx,1)) == 0)
                quad_idx = idx;
                realvalidinside = [realvalidinside;idx];
            else
                realvalidoutside = [realvalidoutside; idx];
            end
            if(isempty(quad_idx))
                emptyquad = [emptyquad;idx];
                continue;
            end
            quad_edge = TO_Edges_VALID(idx, :);
            commonedgeidx{hypo2idx} = quad_idx;
            edge32{hypo2idx} = edgels_32;
            allquad = [allquad, quad_idx'];
            anglediff    = [abs(thresh_ore31_1 - thresh_ore32_1); ...
                            abs(thresh_ore31_1 - thresh_ore32_2); ...
                            abs(thresh_ore31_2 - thresh_ore32_1); ...
                            abs(thresh_ore31_2 - thresh_ore32_2)];
            maxanglediff = max(anglediff);
            minanglediff = min(anglediff);
            if (maxanglediff <=30)
                isparal(hypo2idx,1) = 0;
            end
            isparal(hypo2idx,1) = 1;
%         if(isempty(allquad))
%             indice = [];
%             supported_indices_stack       = [supported_indices_stack; indice];
%             supported_link                = [supported_link, zeros(size(commonedgeidx,2),1)];
%         else
            
            Supports_orient = getSupportedEdgels_Orientation(Tangents_VALID, ...
                          reproj_edge_tgt_gamma3, params, commonedgeidx, isparal);
%         supported_edges_hypo2vali{vi} = Supports_orient.indices;
            supported_indices_stack       = [supported_indices_stack; Supports_orient.indices];
            supported_link                = [supported_link, Supports_orient.link];
            supported_simi                = [supported_simi, Supports_orient.simi];
            reprojerror = [reprojerror; norm(quad_edge(1,1:2) - reproj_edge_pos_gamma3(1,1:2))];
            if(norm(quad_edge(1,1:2) - reproj_edge_pos_gamma3(1,1:2))>=5)
                norm(quad_edge(1,1:2) - reproj_edge_pos_gamma3(1,1:2));
            end
%}
%         end

%     figure(1);
%     subplot(2,2,3)
%     hold on;
%     plot(reproj_edge_pos_gamma3(1,1), reproj_edge_pos_gamma3(1,2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
%     edgel_tgt2  = reproj_edge_tgt_gamma3(1, 1:2);
%     hold on;
%     x1 = reproj_edge_pos_gamma3(1, 1)+len*edgel_tgt2(1,1);
%     x2 = reproj_edge_pos_gamma3(1, 1)-len*edgel_tgt2(1,1);
%     y1 = reproj_edge_pos_gamma3(1, 2)+len*edgel_tgt2(1,2);
%     y2 = reproj_edge_pos_gamma3(1, 2)-len*edgel_tgt2(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.25);
%{
    figure(2);
    subplot(2,2,1)
    imshow(imageArray{params.HYPO1_VIEW_INDX});
    hold on;
    plot(TO_Edges_HYPO1(:,1), TO_Edges_HYPO1(:,2), 'y.', 'MarkerSize',2, 'LineWidth', 1);
    plot(reproject_coor2(1,1), reproject_coor2(1,2), 'gx', 'MarkerSize',10, 'LineWidth', 1);
%     plot(reproject_coor2_pret(1,1), reproject_coor2_pret(1,2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
    plot(edgels_HYPO1_final(1,1), edgels_HYPO1_final(1,2), 'bx', 'MarkerSize',10, 'LineWidth', 1);
    edgel_tgt2  = reproject_tgt2(1, 1:2);
    hold on;
    x1 = reproject_coor2(1, 1)+len*edgel_tgt2(1,1);
    x2 = reproject_coor2(1, 1)-len*edgel_tgt2(1,1);
    y1 = reproject_coor2(1, 2)+len*edgel_tgt2(1,2);
    y2 = reproject_coor2(1, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.25);
    hold on;
%     x1 = reproject_coor2_pret(1, 1)+len*edgel_tgt2(1,1);
%     x2 = reproject_coor2_pret(1, 1)-len*edgel_tgt2(1,1);
%     y1 = reproject_coor2_pret(1, 2)+len*edgel_tgt2(1,2);
%     y2 = reproject_coor2_pret(1, 2)-len*edgel_tgt2(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.25);
    hold on;
    x1 = edgels_HYPO1_final(1, 1)+len*edgel_tgt2(1,1);
    x2 = edgels_HYPO1_final(1, 1)-len*edgel_tgt2(1,1);
    y1 = edgels_HYPO1_final(1, 2)+len*edgel_tgt2(1,2);
    y2 = edgels_HYPO1_final(1, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'b', 'LineWidth', 0.25);
    
    hypo31_all = edgels_HYPO1_final(1,1:2)' - epipole_pix_view1_3(1:2,:);
    b31_all    = hypo31_all(2,:)./hypo31_all(1,:);
    c31_all    = edgels_HYPO1_final(1,2) - b31_all * edgels_HYPO1_final(1,1);
    y31_min1 = c31_all;
    y31_max1 = b31_all*params.cols + c31_all;
    hold on;
    line([0, params.cols], [y31_min1, y31_max1], 'Color', 'c', 'LineWidth', 0.5);
    hold on;
    plot(TO_Edges_HYPO1(:,1), TO_Edges_HYPO1(:,2), 'y.', 'MarkerSize',4, 'LineWidth', 1);
    title(['view ', num2str(params.HYPO1_VIEW_INDX), ' consider as Hypo1'])

    subplot(2,2,2)
    imshow(imageArray{params.HYPO2_VIEW_INDX});
    hold on;
    plot(TO_Edges_HYPO2(:,1), TO_Edges_HYPO2(:,2), 'y.', 'MarkerSize',2, 'LineWidth', 1);
    plot(reproject_coor2(2,1), reproject_coor2(2,2), 'gx', 'MarkerSize',10, 'LineWidth', 1);
%     plot(reproject_coor2_pret(2,1), reproject_coor2_pret(2,2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
    plot(edgels_HYPO2_final(1,1), edgels_HYPO2_final(1,2), 'bx', 'MarkerSize',10, 'LineWidth', 1);
    edgel_tgt2  = reproject_tgt2(2, 1:2);
    hold on;
    x1 = reproject_coor2(2, 1)+len*edgel_tgt2(1,1);
    x2 = reproject_coor2(2, 1)-len*edgel_tgt2(1,1);
    y1 = reproject_coor2(2, 2)+len*edgel_tgt2(1,2);
    y2 = reproject_coor2(2, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.25);
    hold on;
%     x1 = reproject_coor2_pret(2, 1)+len*edgel_tgt2(1,1);
%     x2 = reproject_coor2_pret(2, 1)-len*edgel_tgt2(1,1);
%     y1 = reproject_coor2_pret(2, 2)+len*edgel_tgt2(1,2);
%     y2 = reproject_coor2_pret(2, 2)-len*edgel_tgt2(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.25);
    hold on;
    x1 = edgels_HYPO2_final(1, 1)+len*edgel_tgt2(1,1);
    x2 = edgels_HYPO2_final(1, 1)-len*edgel_tgt2(1,1);
    y1 = edgels_HYPO2_final(1, 2)+len*edgel_tgt2(1,2);
    y2 = edgels_HYPO2_final(1, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'b', 'LineWidth', 0.25);
    hypo32_all = edgels_HYPO2_final(1,1:2)' - epipole_pix_view2_3(1:2,:);
    b32_all    = hypo32_all(2,:)./hypo32_all(1,:);
    c32_all    = edgels_HYPO2_final(1,2) - b32_all * edgels_HYPO2_final(1,1);
    y32_min1   = c32_all;
    y32_max1   = b32_all*params.cols + c32_all;
    hold on;
    line([0, params.cols], [y32_min1, y32_max1], 'Color', 'm', 'LineWidth', 0.5);
    hold on;
    plot(TO_Edges_HYPO2(:,1), TO_Edges_HYPO2(:,2), 'y.', 'MarkerSize',4, 'LineWidth', 2);
    title(['view ', num2str(params.HYPO2_VIEW_INDX), ' consider as Hypo2'])

    subplot(2,2,3)
    imshow(imageArray{VALID_INDX});
    hold on;
    plot(TO_Edges_VALID(:,1), TO_Edges_VALID(:,2), 'y.', 'MarkerSize',2, 'LineWidth', 1);
    plot(reproj_edge_pos_gamma3(1,1), reproj_edge_pos_gamma3(1,2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
    edgel_tgt2  = quad_edge(1, 3:4);
    hold on;
    x1 = reproj_edge_pos_gamma3(1, 1)+len*edgel_tgt2(1,1);
    x2 = reproj_edge_pos_gamma3(1, 1)-len*edgel_tgt2(1,1);
    y1 = reproj_edge_pos_gamma3(1, 2)+len*edgel_tgt2(1,2);
    y2 = reproj_edge_pos_gamma3(1, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.25);
    edgel_tgt2  = reproj_edge_tgt_gamma3(1, 1:2);
    hold on;
    x1 = reproj_edge_pos_gamma3(1, 1)+len*edgel_tgt2(1,1);
    x2 = reproj_edge_pos_gamma3(1, 1)-len*edgel_tgt2(1,1);
    y1 = reproj_edge_pos_gamma3(1, 2)+len*edgel_tgt2(1,2);
    y2 = reproj_edge_pos_gamma3(1, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.25);
    k1 = tan(thresh_ore31_1/180*pi);
    b1 = epipole_pix_view3_1(2,1) - k1*epipole_pix_view3_1(1,1);
    yMax1 = k1*params.cols  + b1;
    % y = k2*x + b2
    k2 = tan(thresh_ore31_2/180*pi);
    b2 = epipole_pix_view3_1(2,1) - k2*epipole_pix_view3_1(1,1);
    yMax2 = k2*params.cols  + b2;
    hold on;
    line([0, params.cols ], [b1, yMax1], 'Color', 'c', 'LineWidth', 0.5);
    hold on;
    line([0, params.cols ], [b2, yMax2], 'Color', 'c', 'LineWidth', 0.5);

    k1 = tan(thresh_ore32_1/180*pi);
    b1 = epipole_pix_view3_2(2,1) - k1*epipole_pix_view3_2(1,1);
    yMax1 = k1*params.cols + b1;
    % y = k2*x + b2
    k2 = tan(thresh_ore32_2/180*pi);
    b2 = epipole_pix_view3_2(2,1) - k2*epipole_pix_view3_2(1,1);
    yMax2 = k2*params.cols + b2;
    hold on;
    line([0, params.cols], [b1, yMax1], 'Color', 'm', 'LineWidth', 0.5);
    hold on;
    line([0, params.cols], [b2, yMax2], 'Color', 'm', 'LineWidth', 0.5);
    hold on;
    plot(TO_Edges_VALID(commonedgeidx{1},1), TO_Edges_VALID(commonedgeidx{1},2), 'b.', 'MarkerSize',4, 'LineWidth', 2);
    edgel_tgt2  = TO_Edges_VALID(commonedgeidx{1},3:4);
    hold on;
    x1 = TO_Edges_VALID(commonedgeidx{1}, 1)+len*edgel_tgt2(1,1);
    x2 = TO_Edges_VALID(commonedgeidx{1}, 1)-len*edgel_tgt2(1,1);
    y1 = TO_Edges_VALID(commonedgeidx{1}, 2)+len*edgel_tgt2(1,2);
    y2 = TO_Edges_VALID(commonedgeidx{1}, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'b', 'LineWidth', 0.25);
    title(['view ', num2str(VALID_INDX)])
    sgtitle (['epipolar wedge with Δ=',num2str(params.delta),' pixels'])
%}
if(supported_link(1,end)==0)
    supported_link(1,end);
end
end
% figure
% histogram(reprojerror);
% title 'distribution of reprojection error'
% xlabel 'reprojection error in pixel'
% ylabel 'number of occurrence'
losevalid = find(supported_link == 0);
ratioofvalid = (size(supported_link,2)-size(losevalid,2))/size(supported_link,2);
notinquad    = (size(realvalidoutside,1))/size(recons_coor2,1);

fprintf('epipolar wedge with Δ= ±%f \n', ...
params.delta);
fprintf('Perturb edge with in a cirlce of radius: %f pixels\n', ...
pix_perturb);
% fprintf('Number of true correspondence being validated: %d \n', ...
% (size(supported_link,2)-size(losevalid,2)));
fprintf('Percentage of true correspondence being validated: %f %%\n', ...
ratioofvalid*100);
% if(isempty(emptyquad) ~= 1)
%     fprintf('Number of edge pair formed empty quadrilateral: %d \n', ...
%         size(emptyquad,1));
%     fprintf('edge idx: ');
%     fprintf('%d ', emptyquad');
%     fprintf('\n');
% end

pix_pert        = [pix_pert;        pix_perturb];
percentage_vali = [percentage_vali; ratioofvalid*100];

% fprintf('Percentage of true correspondence outside the quadrilateral in validation view: %f %%\n', ...
% notinquad*100);

end
figure;
plot(pix_pert, percentage_vali, 'b-o');
xlabel 'radius of 2D edge being perturbed (in pixel)'
ylabel 'percentage of true correspondence being validated'
ylim ([97 100])



%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perturb 2D edges in the view 1 and 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_pert = [];
T_pert = [];
reproject_coor2_pret = [reproject_coor2(:,1)*(0.99 + (0.02)*rand(1)), ...
                        reproject_coor2(:,2)*(0.99 + (0.02)*rand(1))];
reproject_tgt2_pret   = reproject_tgt2;
subnum = 1;
figure(2);
for(vidx = 1:3)
    vi = validateidex(1,vidx);
    VALID_INDX = vi;%validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    subplot(2,2,subnum)
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

refidx = validateidex;
newidx = [validateidex(1,2), validateidex(1,3), validateidex(1,1);...
          validateidex(1,3), validateidex(1,1), validateidex(1,2);];
subidx = [2, 3, 1; 3, 1, 2]; 
colorlin3 = [0 1 1;
             1 0 1;
             0.4660 0.6740 0.1880];
R_relative = [];
T_relative = [];

for(vidx = 1:3)
    idx_ref = refidx(1,vidx);
    idx_new = newidx(1,vidx);
    abs_R1 =  R_g(:,:,idx_ref);
    abs_C1 = -R_g(:,:,idx_ref)' * T_g(:,:,idx_ref);
    abs_R2 =  R_g(:,:,idx_new);
    abs_C2 = -R_g(:,:,idx_new)' * T_g(:,:,idx_new);

    % CPP modification
    R21    = abs_R2 * abs_R1';
    % CPP modification
    T21    = abs_R2 * (abs_C1 - abs_C2);
    %%%
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
    subplot(2,2,subidx(1,vidx))
    hold on;
    line([0, params.cols], [yMin, yMax], 'Color', colorlin3(vidx,:), 'LineWidth', 0.5);

    %%%%%
    idx_new = newidx(2,vidx);
    abs_R2 =  R_g(:,:,idx_new);
    abs_C2 = -R_g(:,:,idx_new)' * T_g(:,:,idx_new);
    % CPP modification
    R21    = abs_R2 * abs_R1';
    % CPP modification
    T21    = abs_R2 * (abs_C1 - abs_C2);
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
    subplot(2,2,subidx(2,vidx))
    hold on;
    line([0, params.cols], [yMin, yMax], 'Color', colorlin3(vidx,:), 'LineWidth', 0.5);
end
sgtitle 'both 2D edges and poses(R and T) being perturbed'

%}