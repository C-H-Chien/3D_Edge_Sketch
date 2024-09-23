% Epipolar wedges method is used in this main code.
clear all; 
close all;

load('MATData/TO_Edges_syn.mat');
load('MATData/imageArray_syn.mat')
load('MATData/R_matrix_syn.mat')
load('MATData/T_matrix_syn.mat')
load('MATData/K_syn.mat')

% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 2;
params.HYPO2_VIEW_INDX          = 6;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.99;
params.MAX_NUM_OF_SUPPORT_VIEWS = 4;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.005*5;
params.ICL_DATA                 = 0;
params.multiK                   = 0;

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

[R21, T21, E, F] = getRelativePose(R_matrix, T_matrix, params, K);

finalEdgePair    = [];
paired_edge      = [];

% get the orientation list between hypo1 and hypo2
[ore_list1bar, ore_list2_sorted, epipole_pix_view1, epipole_pix_view2]= ...
    getOreList(TO_Edges_HYPO1, TO_Edges_HYPO2, R_matrix, T_matrix, params, K);
tic
% pipeline starts
cnt_sup   = 0;
sup_num   = [];
ambi_edge = [];
dist_52_57 = cell(1,8);
dist_to_epipolarline = [];
for edge_idx = 1: size(TO_Edges_HYPO1, 1)
    if mod(edge_idx, 10) == 0,  fprintf(". %d",edge_idx); end
    if mod(edge_idx, 500) == 0, fprintf("\n"); end
    edgel_HYPO1 = TO_Edges_HYPO1(edge_idx, :);
    if edgel_HYPO1(1,1) < 10 || edgel_HYPO1(1,1) > cols-10 || ...
            edgel_HYPO1(1,2) < 10 || edgel_HYPO1(1,2) > rows-10
        continue;
    end
%     if(isempty(find(ambi_edge_prev == edge_idx)) == 1)
%         continue
%     end
    coeffs = F * [edgel_HYPO1(1,1:2)'; 1];
    Apixel = coeffs(1,1);
    Bpixel = coeffs(2,1);
    Cpixel = coeffs(3,1);
    Epipolar_Coeffs.A = Apixel;
    Epipolar_Coeffs.B = Bpixel;
    Epipolar_Coeffs.C = Cpixel;
    % find epipolar angle range for hypo 2
    angle_range1 = abs(ore_list2_sorted(1,1) - ...
        ore_list2_sorted(size(ore_list2_sorted,1),1));
    % get the range for single epipolar wedge in hypo 2
    range1       = params.PERCENTAGE_FOR_EPIPOLE*angle_range1;
    % get the epipolar wedge range for the edge in hypo1 on hypo2
    thresh_ore1  = ore_list1bar(edge_idx,1)-range1;
    thresh_ore2  = ore_list1bar(edge_idx,1)+range1;
    % get edges in hypo2 fall inside the epipolar wedge
    idxlist      = find(ore_list2_sorted(:,1) >= thresh_ore1 & ...
                       ore_list2_sorted(:,1) <= thresh_ore2);
    if(isempty(idxlist))
        continue;
    end
    hypo2_idx    = ore_list2_sorted(idxlist,2);
    edgels_HYPO2 = TO_Edges_HYPO2(sort(hypo2_idx), :);
%     edgels_HYPO2 = TO_Edges_HYPO2;

    dist_52_57matr = [];
    dist_52_57list = [];
    dist_to_line = abs(Epipolar_Coeffs.A.*edgels_HYPO2(:,1)+Epipolar_Coeffs.B.*edgels_HYPO2(:,2)+Epipolar_Coeffs.C)./sqrt(Epipolar_Coeffs.A^2+Epipolar_Coeffs.B^2);
    dist_to_epipolarline = [dist_to_epipolarline; dist_to_line];
%     figure;
    for vi = 1:size(validation_view_indices, 1)
        % get the validation view
        VALID_INDX = validation_view_indices(vi, 1);
        VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
        % get edges in the validation view
        TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
        Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
        [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgel ...
             (edgel_HYPO1, edgels_HYPO2, R_matrix, T_matrix, params, VALID_INDX, K);
        x_1 = reproj_edge_pos_gamma3(:,1);
        y_1 = reproj_edge_pos_gamma3(:,2);
        x_2 = reproj_edge_pos_gamma3(:,3);
        y_2 = reproj_edge_pos_gamma3(:,4);
        dist_52_57_curr = sqrt((x_1-x_2).^2+(y_1-y_2).^2);
        dist_52_57{vi} = [dist_52_57{vi};dist_52_57_curr];
        dist_52_57matr = [dist_52_57matr, dist_52_57_curr];
        dist_52_57list = [dist_52_57list; dist_52_57_curr];
%         subplot(3,3,vi)
%         plot(dist_to_epipolarline, dist_52_57_curr, 'b.', 'MarkerSize',1, 'LineWidth', 1);
%         xlabel 'distance to the epipolar line from \gamma 2'
%         ylabel ({'norm error between \gamma 3 (eq 52)';' and \gamma3 bar (eq 57)'})
%         title(['view ', num2str(validation_view_indices(vi)), '(vali view ', num2str(vi), ')'])
    end

end

for(iv = 1:8)
    subplot(3,3,iv)
    plot(dist_to_epipolarline, dist_52_57{iv}, 'b.', 'MarkerSize',1, 'LineWidth', 1);
    xlabel ('distance to the epipolar line from $\gamma 2$','fontsize',10,'interpreter','latex')
    ylabel ({'norm error between $\gamma 3$ (eq 52)';' and $\bar{\gamma 3}$ bar (eq 57)'},'fontsize',10,'interpreter','latex')
    title(['view ', num2str(validation_view_indices(iv)), '(vali view ', num2str(iv), ')'],'fontsize',10,'interpreter','latex')
end
sgtitle (['Number of $\gamma 1$ =', num2str(size(TO_Edges_HYPO1, 1)),', Number of $\gamma 3$ and $\bar{\gamma 3}$ =', num2str(size(dist_52_57{vi}, 1))],'fontsize',10,'interpreter','latex')