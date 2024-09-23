hypo2idx = 6;
for(vi = 1:48)
    % get the validation view
    VALID_INDX = validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    % get edges in the validation view
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,6};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    % get the orientation list between hypo1 and validation view
    %         [ore_31bar, ore_list31_sorted, epipole_pix_view1_3, epipole_pix_view3_1] = ...
    %             getOreList_Vali(edgels_HYPO1_final, TO_Edges_VALID, R_matrix, T_matrix, ...
    %             params.HYPO1_VIEW_INDX, K, VALID_INDX, params);
    [ore_31bar_all, ...
        ore_list31_sorted, ...
        epipole_pix_view1_3, ...
        epipole_pix_view3_1] = getOreList_New(edgels_HYPO1_final, ...
        TO_Edges_VALID, ...
        R_matrix, T_matrix, ...
        params.HYPO1_VIEW_INDX, ...
        VALID_INDX, params, K);
    % find epipolar angle range for validation view
    %         angle_range2 = abs(ore_list31_sorted(1,1)- ...
    %             ore_list31_sorted(size(ore_list31_sorted,1),1));
    % get the range for single epipolar wedge in vali
    %         range2       = params.ANGLE_FOR_EPIPOLE; %PERCENTAGE_FOR_EPIPOLE*angle_range2;
    %         [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgel ...
    %              (edgel_HYPO1, edgels_HYPO2_final, R_matrix, T_matrix, params, VALID_INDX, K);
    [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgelCorrected ...
        (edgels_HYPO1_final, edgels_HYPO2_final, R_matrix, T_matrix, params, VALID_INDX, K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the orientation list between hypo2 and validation view
    %         [ore_32bar, ore_list32_sorted, epipole_pix_view2_3, epipole_pix_view3_2] = ...
    %             getOreList_Vali(edgels_HYPO2_final, TO_Edges_VALID, R_matrix, T_matrix, ...
    %             params.HYPO2_VIEW_INDX, K, VALID_INDX, params);
    [ore_32bar_all, ...
        ore_list32_sorted, ...
        epipole_pix_view2_3, ...
        epipole_pix_view3_2] = getOreList_New(edgels_HYPO2_final, ...
        TO_Edges_VALID, ...
        R_matrix, T_matrix, ...
        params.HYPO2_VIEW_INDX, ...
        VALID_INDX, params, K);
    % find epipolar angle range for validation view
    %         angle_range3 = abs(ore_list32_sorted(1,1) - ...
    %             ore_list32_sorted(size(ore_list32_sorted,1),1));
    %         % get the range for single epipolar wedge in vali
    %         range3       = params.ANGLE_FOR_EPIPOLE; %PERCENTAGE_FOR_EPIPOLE*angle_range3;
    % investigate all the potential pairs
    % "form" the quadrilateral
    allquad = [];
    isparal = ones(size(edgels_HYPO2_final,1),1);

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
        quad_edge = TO_Edges_VALID(quad_idx, :);
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
        commonedge{1} = commonedgeidx{hypo2idx};
    %         if(isempty(allquad))
    %             indice = [];
    %             supported_indices_stack       = [supported_indices_stack; indice];
    %             supported_link                = [supported_link, zeros(size(commonedgeidx,2),1)];
    %         else
    Supports_orient = getSupportedEdgels_Orientation(Tangents_VALID, ...
        reproj_edge_tgt_gamma3(hypo2idx,:), params, commonedge, isparal(hypo2idx,1));
    %         supported_edges_hypo2vali{vi} = Supports_orient.indices;
%     supported_indices_stack  = [supported_indices_stack; Supports_orient.indices];
%     supported_link          = [supported_link, Supports_orient.link];
%     supported_simi          = [supported_simi, Supports_orient.simi];

    %         end


len = 1;
% VALID_INDX = params.HYPO2_VIEW_INDX;
figure(vi);
% % % imshow(imageArray{VALID_INDX});
hold on;
% plot(TO_Edges_VALID(:,1), TO_Edges_VALID(:,2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
candidateindex = hypo2idx;
plot(reproj_edge_pos_gamma3(candidateindex,1), reproj_edge_pos_gamma3(candidateindex,2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
edgel_tgt2  = reproj_edge_tgt_gamma3(candidateindex, :);
hold on;
x1 = reproj_edge_pos_gamma3(candidateindex, 1)+len*edgel_tgt2(1,1);
x2 = reproj_edge_pos_gamma3(candidateindex, 1)-len*edgel_tgt2(1,1);
y1 = reproj_edge_pos_gamma3(candidateindex, 2)+len*edgel_tgt2(1,2);
y2 = reproj_edge_pos_gamma3(candidateindex, 2)-len*edgel_tgt2(1,2);
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
line([0, params.cols ], [b1, yMax1], 'Color', 'y', 'LineWidth', 0.5);
hold on;
line([0, params.cols ], [b2, yMax2], 'Color', 'y', 'LineWidth', 0.5);
% 31

abs_R1v =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C1v = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
           T_matrix(:,params.HYPO1_VIEW_INDX);
abs_R2v =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C2v = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' *...
           T_matrix(:,params.HYPO2_VIEW_INDX);
abs_R3v =  R_matrix(:,:,VALID_INDX);
abs_C3v = -R_matrix(:,:,VALID_INDX)' *...
           T_matrix(:,VALID_INDX);
R31v    = abs_R3v * abs_R1v';
T31v    = abs_R3v * (abs_C1v - abs_C3v);
R32v    = abs_R3v * abs_R2v';
T32v    = abs_R3v * (abs_C2v - abs_C3v);
%> Calculate Essential matrix
T_x = @(T)[0,      -T(3,1),  T(2,1); ...
           T(3,1),  0,      -T(1,1); ...
          -T(2,1),  T(1,1),  0];
E31v   = T_x(T31v) * R31v;
E32v   = T_x(T32v) * R32v;
% Calculate fundamental matrix
if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    K3 = K(:,:,VALID_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
    invK3 = inv(K3);
    F31v  = invK3'*E31v*invK1;
    F32v  = invK3'*E32v*invK2;
else
    invK = inv(K);
    F31v = invK'*E31v*invK;
    F32v = invK'*E32v*invK;
end
coeffs31v = F31v * [edgels_HYPO1_final(candidateindex,1:2)'; 1];
Apixel1 = coeffs31v(1,1);
Bpixel1 = coeffs31v(2,1);
Cpixel1 = coeffs31v(3,1);
Epipolar_Coeffs31.A = Apixel1;
Epipolar_Coeffs31.B = Bpixel1;
Epipolar_Coeffs31.C = Cpixel1;
yMin = -Epipolar_Coeffs31.C./Epipolar_Coeffs31.B;
yMax = (-Epipolar_Coeffs31.C - Epipolar_Coeffs31.A*cols) ./ Epipolar_Coeffs31.B;
line([0, cols], [yMin, yMax], 'Color', 'y', 'LineWidth', 0.5,'LineStyle','--');
hold on;
plot(TO_Edges_VALID(vali_idx,1), TO_Edges_VALID(vali_idx,2), 'ydiamond', 'MarkerSize',5, 'LineWidth', 1);

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
% hold on;
% plot(TO_Edges_VALID(commonedgeidx{candidateindex},1), TO_Edges_VALID(commonedgeidx{1},2), 'gx', 'MarkerSize',5, 'LineWidth', 2);
% edgel_tgt2  = TO_Edges_VALID(commonedgeidx{candidateindex},3:4);
% hold on;
% x1 = TO_Edges_VALID(commonedgeidx{candidateindex}, 1)+len*edgel_tgt2(1,1);
% x2 = TO_Edges_VALID(commonedgeidx{candidateindex}, 1)-len*edgel_tgt2(1,1);
% y1 = TO_Edges_VALID(commonedgeidx{candidateindex}, 2)+len*edgel_tgt2(1,2);
% y2 = TO_Edges_VALID(commonedgeidx{candidateindex}, 2)-len*edgel_tgt2(1,2);
% hold on;
% plot([x1 x2], [y1 y2], 'b', 'LineWidth', 0.25);



% 32
coeffs32v = F32v * [edgels_HYPO2_final(candidateindex,1:2)'; 1];
Apixel2 = coeffs32v(1,1);
Bpixel2 = coeffs32v(2,1);
Cpixel2 = coeffs32v(3,1);
Epipolar_Coeffs32.A = Apixel2;
Epipolar_Coeffs32.B = Bpixel2;
Epipolar_Coeffs32.C = Cpixel2;
yMin = -Epipolar_Coeffs32.C./Epipolar_Coeffs32.B;
yMax = (-Epipolar_Coeffs32.C - Epipolar_Coeffs32.A*cols) ./ Epipolar_Coeffs32.B;
line([0, cols], [yMin, yMax], 'Color', 'm', 'LineWidth', 0.5,'LineStyle','--');

hold on;
plot(TO_Edges_VALID(vali_idx1,1), TO_Edges_VALID(vali_idx1,2), 'ms', 'MarkerSize',7, 'LineWidth', 1);
hold on;
if(isempty(quad_idx))
    plot([-1,401], [401,-1], 'gx', 'MarkerSize',5, 'LineWidth', 1);
    legend '2D edges from this view' ...
       'γ3' 'γ3 orientation' ...
       'epipolar wedge between hypo1 and this view' '' ...
       'epipolar line of corrected γ1' ...
       'edge from this view inside yellow wedge' ...
       'epipolar wedge between hypo2 and this view' '' ...
       'epipolar line of corrected γ2' ...
       'edge from this view inside pink wedge' ...
       'edge from this view inside both wedges' ''
title(['example of validation view - view ', num2str(VALID_INDX), ', quadrilateral is empty'])
close(vi)
else
% plot(TO_Edges_VALID(quad_idx,1), TO_Edges_VALID(quad_idx,2), 'gx', 'MarkerSize',10, 'LineWidth', 2);

tgt3 = Tangents_VALID;
% for(idxvi = 1: size (quad_idx,1))
%     edgel_tgt3  = tgt3(quad_idx(idxvi), 1:2);
%     edgels_vali = TO_Edges_VALID(quad_idx(idxvi),1:2);
%     hold on;
%     x1 = edgels_vali(1, 1)+len*edgel_tgt3(1,1);
%     x2 = edgels_vali(1, 1)-len*edgel_tgt3(1,1);
%     y1 = edgels_vali(1, 2)+len*edgel_tgt3(1,2);
%     y2 = edgels_vali(1, 2)-len*edgel_tgt3(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.5);
% end
legend '2D edges from this view' ...
       'γ3' 'γ3 orientation' ...
       'epipolar wedge between hypo1 and this view' '' ...
       'epipolar line of corrected γ1' ...
       'edge from this view inside yellow wedge' ...
       'epipolar wedge between hypo2 and this view' '' ...
       'epipolar line of corrected γ2' ...
       'edge from this view inside pink wedge' ...
       'edge from this view inside both wedges' ''
title(['example of validation view - view ', num2str(VALID_INDX)])
end

end