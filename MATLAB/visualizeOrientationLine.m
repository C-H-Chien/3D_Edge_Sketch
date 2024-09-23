chosen_idx = finalPairIndx;
% paired_edge1 = [edge_idx, hypo2_idx(chosen_idx), supported_link(chosen_idx,:)];
cols = size(hypo_img2, 2);
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
% figure
fignum = 2;
subnum = 1;
for (vi = 1:6)
    % get the validation view
    VALID_INDX = validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    % get edges in the validation view
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    % get the orientation list between hypo1 and validation view
    [ore_31bar, ore_list31_sorted, epipole_pix_view1_3, epipole_pix_view3_1] = ...
            getOreList_Vali(edgel_HYPO1, TO_Edges_VALID, R_matrix, T_matrix, ...
            params.HYPO1_VIEW_INDX, K, VALID_INDX, params);
    % find epipolar angle range for validation view
    angle_range2   = abs(ore_list31_sorted(1,1)- ...
        ore_list31_sorted(size(ore_list31_sorted,1),1));
    % get the range for single epipolar wedge in vali
    range2         = params.PERCENTAGE_FOR_EPIPOLE*angle_range2;
    % get the epipolar wedge range for the edge in hypo1 on vali
    thresh_ore31_1 = ore_31bar-range2;
    thresh_ore31_2 = ore_31bar+range2;
    idxlist_31     = find(ore_list31_sorted(:,1) >= thresh_ore31_1 & ...
        ore_list31_sorted(:,1) <= thresh_ore31_2);
    vali_idx       = ore_list31_sorted(idxlist_31,2);
    edgels_31      = TO_Edges_VALID(vali_idx, :);
    [reproj_edge_pos_gamma3, reproj_edge_tgt_gamma3] = getReprojectedEdgel ...
        (edgel_HYPO1, edgels_HYPO2, R_matrix, T_matrix, params, VALID_INDX, K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the orientation list between hypo2 and validation view
    [ore_32bar, ore_list32_sorted, epipole_pix_view2_3, epipole_pix_view3_2] = ...
        getOreList_Vali(edgels_HYPO2, TO_Edges_VALID, R_matrix, T_matrix, ...
        params.HYPO2_VIEW_INDX, K, VALID_INDX, params);
    % find epipolar angle range for validation view
    angle_range3 = abs(ore_list32_sorted(1,1) - ...
        ore_list32_sorted(size(ore_list32_sorted,1),1));
    % get the range for single epipolar wedge in vali
    range3       = params.PERCENTAGE_FOR_EPIPOLE*angle_range3;
    % investigate all the potential pairs
    allquad      = [];
    for(hypo2idx = 1:size(edgels_HYPO2,1))
        % hypo2idx  = chosen_idx;
        % get the epipolar wedge range for the edge in hypo2 on vali
        thresh_ore32_1 = ore_32bar(hypo2idx,1)-range3;
        thresh_ore32_2 = ore_32bar(hypo2idx,1)+range3;
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

        if(hypo2idx  == chosen_idx)
            anglediff    = [abs(thresh_ore31_1 - thresh_ore32_1); ...
                abs(thresh_ore31_1 - thresh_ore32_2); ...
                abs(thresh_ore31_2 - thresh_ore32_1); ...
                abs(thresh_ore31_2 - thresh_ore32_2)];
            maxanglediff = max(anglediff);
            minanglediff = min(anglediff);
        end

    end
    Supports_orient = getSupportedEdgels_Orientation_1(Tangents_VALID, ...
        reproj_edge_tgt_gamma3, params, commonedgeidx,chosen_idx);


%     if(supported_link(chosen_idx,vi) ~= 0)
%         (ceil(vi/4)+18)
%     if(mod(vi,4) == 0)
%         subplot(2,2,4);
%     else
%         subplot(2,2,mod(vi,4));
%     end
%     subplot(5,5,vi);

    figure(fignum)
    subplot(2,3,subnum)
    subnum = subnum+1;

    imshow(rgb2gray(imageArray{VALID_INDX}));
    hold on;
    plot(TO_Edges_VALID(:,1), TO_Edges_VALID(:,2), 'c.', 'MarkerSize',5, 'LineWidth', 1);
    title(['view ', num2str(validation_view_indices(vi)), '(vali view ', num2str(vi), ')'])
    [R31, T31, E31, F31] = getRelativePose1(R_matrix, T_matrix,params.HYPO1_VIEW_INDX,VALID_INDX,params, K);
    [R32, T32, E32, F32] = getRelativePose1(R_matrix, T_matrix,params.HYPO2_VIEW_INDX,VALID_INDX,params, K);
    coeffs31 = F31 * [edgel_HYPO1(1,1:2)'; 1];
    Apixel31 = coeffs31(1,1);
    Bpixel31 = coeffs31(2,1);
    Cpixel31 = coeffs31(3,1);
    Epipolar_Coeffs31.A = Apixel31;
    Epipolar_Coeffs31.B = Bpixel31;
    Epipolar_Coeffs31.C = Cpixel31;
    yMin31 = -Epipolar_Coeffs31.C./Epipolar_Coeffs31.B;
    yMax31 = (-Epipolar_Coeffs31.C - Epipolar_Coeffs31.A*cols) ./ Epipolar_Coeffs31.B;
    line([1, cols], [yMin31, yMax31], 'Color', 'r', 'LineWidth', 1); hold on;
    hold on;

    coeffs32 = F32 * [edgels_HYPO2(chosen_idx, 1:2)'; 1];
    Apixel32 = coeffs32(1,1);
    Bpixel32 = coeffs32(2,1);
    Cpixel32 = coeffs32(3,1);
    Epipolar_Coeffs32.A = Apixel32;
    Epipolar_Coeffs32.B = Bpixel32;
    Epipolar_Coeffs32.C = Cpixel32;
    yMin32 = -Epipolar_Coeffs32.C./Epipolar_Coeffs32.B;
    yMax32 = (-Epipolar_Coeffs32.C - Epipolar_Coeffs32.A*cols) ./ Epipolar_Coeffs32.B;
    line([1, cols], [yMin32, yMax32], 'Color', 'b', 'LineWidth', 1); hold on;
    hold on;


%     if(size(edgels_31,1) ~= 0)
%         hypo31_1 = edgels_31(1,1:2)' - epipole_pix_view3_1(1:2,:);
%         b31_1    = hypo31_1(2,:)./hypo31_1(1,:);
%         c31_1    = edgels_31(1,2) - b31_1 * edgels_31(1,1);
%         y31_min1 = b31_1*1 + c31_1;
%         y31_max1 = b31_1*cols + c31_1;
%         line([1, cols], [y31_min1, y31_max1], 'Color', 'r', 'LineWidth', 0.5);
%         hold on;
%         valiend = size(edgels_31,1);
%         hypo31_2  = edgels_31(valiend,1:2)' - epipole_pix_view3_1(1:2,:);
%         b31_2     = hypo31_2(2,:)./hypo31_2(1,:);
%         c31_2     = edgels_31(valiend,2) - b31_2 * edgels_31(valiend,1);
%         y31_min2 = b31_2*1 + c31_2;
%         y31_max2 = b31_2*cols + c31_2;
%         line([1, cols], [y31_min2, y31_max2], 'Color', 'r', 'LineWidth', 0.5);
%             hold on;
%             set(gcf,'color','w');
%     end
% 
%     edgels_32 = edge32{chosen_idx};
%     if(size(edgels_32,1) ~= 0)
%         hypo32_1 = edgels_32(1,1:2)' - epipole_pix_view3_2(1:2,:);
%         b32_1 = hypo32_1(2,:)./hypo32_1(1,:);
%         c32_1 = edgels_32(1,2) - b32_1 * edgels_32(1,1);
%         y32_min1 = b32_1*1 + c32_1;
%         y32_max1 = b32_1*cols + c32_1;
%         line([1, cols], [y32_min1, y32_max1], 'Color', 'b', 'LineWidth', 0.5);
%         hold on;
%         valiend1 = size(edgels_32,1);
%         hypo32_2 = edgels_32(valiend1,1:2)' - epipole_pix_view3_2(1:2,:);
%         b32_2 = hypo32_2(2,:)./hypo32_2(1,:);
%         c32_2 = edgels_32(valiend1,2) - b32_2 * edgels_32(valiend1,1);
%         y32_min2 = b32_2*1 + c32_2;
%         y32_max2 = b32_2*cols + c32_2;
%         line([1, cols], [y32_min2, y32_max2], 'Color', 'b', 'LineWidth', 0.5);
%         hold on;
%             set(gcf,'color','w');
%     end
% 
    if(subnum > 6)
        if(chosen_idx == finalPairIndx)
            sgtitle ('false match')
        else
            sgtitle ('true match')
        end
        subnum = 1;
        fignum = fignum+1;
    end

    if(isempty(commonedgeidx{chosen_idx}))
        continue
    end
    supported_edge_pos_gamma3 = TO_Edges_VALID(commonedgeidx{chosen_idx}, :);
    supported_edge_tgt_gamma3 = Tangents_VALID(commonedgeidx{chosen_idx}, :);
    len = 3;
% % % %     for(i = 1: size(commonedgeidx{chosen_idx},1))
% % % %         index_current = i;%commonedgeidx{chosen_idx}(i,1);
% % % %         hold on;
% % % %         plot(supported_edge_pos_gamma3(index_current,1), supported_edge_pos_gamma3(index_current,2), 'ys', 'MarkerSize',5, 'LineWidth', 1);
% % % %         hold on;
% % % %         x1 = supported_edge_pos_gamma3(index_current,1)+len*supported_edge_tgt_gamma3(index_current,1);
% % % %         x2 = supported_edge_pos_gamma3(index_current,1)-len*supported_edge_tgt_gamma3(index_current,1);
% % % %         y1 = supported_edge_pos_gamma3(index_current,2)+len*supported_edge_tgt_gamma3(index_current,2);
% % % %         y2 = supported_edge_pos_gamma3(index_current,2)-len*supported_edge_tgt_gamma3(index_current,2);
% % % %         plot([x1 x2], [y1 y2], 'y', 'LineWidth', 0.5);
% % % %     end
    pair_tgt = reproj_edge_tgt_gamma3(chosen_idx,:);
    pair_pt  = reproj_edge_pos_gamma3(chosen_idx,:);
    
    x1 = pair_pt(1,1)+len*pair_tgt(1,1);
    x2 = pair_pt(1,1)-len*pair_tgt(1,1);
    y1 = pair_pt(1,2)+len*pair_tgt(1,2);
    y2 = pair_pt(1,2)-len*pair_tgt(1,2);
if(chosen_idx == finalPairIndx)
    hold on;
    plot(pair_pt(1,1), pair_pt(1,2), 'rs', 'MarkerSize',5, 'LineWidth', 1);
    hold on;
    plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.5);
else
    hold on;
    plot(pair_pt(1,1), pair_pt(1,2), 'bs', 'MarkerSize',5, 'LineWidth', 1);
    hold on;
    plot([x1 x2], [y1 y2], 'b', 'LineWidth', 0.5);
end
    if(supported_link(chosen_idx,vi) ~= 0)
        support_tgt = Tangents_VALID(supported_link(chosen_idx,vi),:);
        support_pt  = TO_Edges_VALID(supported_link(chosen_idx,vi),:);
        hold on;
% % % %         plot(support_pt(1,1), support_pt(1,2), 'ms', 'MarkerSize',5, 'LineWidth', 1);
        hold on;
        x1 = support_pt(1,1)+len*support_tgt(1,1);
        x2 = support_pt(1,1)-len*support_tgt(1,1);
        y1 = support_pt(1,2)+len*support_tgt(1,2);
        y2 = support_pt(1,2)-len*support_tgt(1,2);
% % % %         plot([x1 x2], [y1 y2], 'm', 'LineWidth', 0.5);
        title(['view ', num2str(validation_view_indices(vi)), '(vali view ', num2str(vi), ') support'])
    end
    if (maxanglediff <=30 && minanglediff <=30 )
       title(['parallel-discard'])
%        title(['view ', num2str(validation_view_indices(vi_test)), '(vali view ', num2str(vi_test), ') parallel'])
       fprintf('View %d is having parallel epipolar wedges.\n',validation_view_indices(vi));
       fprintf('Max angle difference(in degree): %f\n',maxanglediff);
       fprintf('Min angle difference(in degree): %f\n',minanglediff);
    end
end