chosen_idx = 1;
% paired_edge1 = [edge_idx, hypo2_idx(chosen_idx), supported_link(chosen_idx,:)];
cols = size(hypo_img2, 2);
hypo_img1 = double(rgb2gray(imageArray{params.HYPO1_VIEW_INDX}));
hypo_img2 = double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX}));
% figure
for (vi_test = 1:48)
    VALID_INDX = validation_view_indices(vi_test, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    % get edges in the validation view
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    % get the orientation list between hypo1 and validation view
        [ore_31bar, ore_list31_sorted, epipole_pix_view1_3, epipole_pix_view3_1] = ...
            getOreList_Vali(edgel_HYPO1, TO_Edges_VALID, R_matrix, T_matrix, ...
            params.HYPO1_VIEW_INDX, K, VALID_INDX, params);
        % find epipolar angle range for validation view
        angle_range2 = abs(ore_list31_sorted(1,1)- ...
            ore_list31_sorted(size(ore_list31_sorted,1),1));
        % get the range for single epipolar wedge in vali
        range2       = params.PERCENTAGE_FOR_EPIPOLE*angle_range2;
        % get the epipolar wedge range for the edge in hypo1 on vali
        thresh_ore31_1 = ore_31bar-range2;
        thresh_ore31_2 = ore_31bar+range2;
        idxlist_31      = find(ore_list31_sorted(:,1) >= thresh_ore31_1 & ...
            ore_list31_sorted(:,1) <= thresh_ore31_2);
        vali_idx     = ore_list31_sorted(idxlist_31,2);
        edgels_31    = TO_Edges_VALID(vali_idx, :);
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
        % "form" the quadrilateral
        allquad = [];
        isparal = ones(size(edgels_HYPO2,1),1);
%         for(hypo2idx = 1:size(edgels_HYPO2,1))
            hypo2idx  = chosen_idx;
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
    anglediff    = [abs(thresh_ore31_1 - thresh_ore32_1); ...
                            abs(thresh_ore31_1 - thresh_ore32_2); ...
                            abs(thresh_ore31_2 - thresh_ore32_1); ...
                            abs(thresh_ore31_2 - thresh_ore32_2)];
            maxanglediff = max(anglediff);
            minanglediff = min(anglediff);
    
%     subplot(6,8,vi_test);
figure
    imshow(uint8(VALID_IMG));
%     hold on;
%     title(['view ', num2str(validation_view_indices(vi_test)), '(vali view ', num2str(vi_test), ')'])
    hold on;
    plot(TO_Edges_VALID(:, 1), TO_Edges_VALID(:, 2), 'c.', 'MarkerSize',0.25);
    % y = bx + c
    if(size(edgels_31,1) ~= 0)
        hypo31_1 = edgels_31(1,1:2)' - epipole_pix_view3_1(1:2,:);
        b31_1    = hypo31_1(2,:)./hypo31_1(1,:);
        c31_1    = edgels_31(1,2) - b31_1 * edgels_31(1,1);
        y31_min1 = b31_1*1 + c31_1;
        y31_max1 = b31_1*cols + c31_1;
        line([1, cols], [y31_min1, y31_max1], 'Color', 'r', 'LineWidth', 0.5);
        hold on;
        valiend = size(edgels_31,1);
        hypo31_2  = edgels_31(valiend,1:2)' - epipole_pix_view3_1(1:2,:);
        b31_2     = hypo31_2(2,:)./hypo31_2(1,:);
        c31_2     = edgels_31(valiend,2) - b31_2 * edgels_31(valiend,1);
        y31_min2 = b31_2*1 + c31_2;
        y31_max2 = b31_2*cols + c31_2;
        line([1, cols], [y31_min2, y31_max2], 'Color', 'r', 'LineWidth', 0.5);
            hold on;
            set(gcf,'color','w');
    end

    edgels_32 = edge32{chosen_idx};
    if(size(edgels_32,1) ~= 0)
        hypo32_1 = edgels_32(1,1:2)' - epipole_pix_view3_2(1:2,:);
        b32_1 = hypo32_1(2,:)./hypo32_1(1,:);
        c32_1 = edgels_32(1,2) - b32_1 * edgels_32(1,1);
        y32_min1 = b32_1*1 + c32_1;
        y32_max1 = b32_1*cols + c32_1;
        line([1, cols], [y32_min1, y32_max1], 'Color', 'w', 'LineWidth', 0.5);
        hold on;
        valiend1 = size(edgels_32,1);
        hypo32_2 = edgels_32(valiend1,1:2)' - epipole_pix_view3_2(1:2,:);
        b32_2 = hypo32_2(2,:)./hypo32_2(1,:);
        c32_2 = edgels_32(valiend1,2) - b32_2 * edgels_32(valiend1,1);
        y32_min2 = b32_2*1 + c32_2;
        y32_max2 = b32_2*cols + c32_2;
        line([1, cols], [y32_min2, y32_max2], 'Color', 'w', 'LineWidth', 0.5);
        hold on;
            set(gcf,'color','w');
    end
    title(['view ', num2str(validation_view_indices(vi_test)), '(vali view ', num2str(vi_test), ')'])
    quad_edge1 =[];
    if(isempty(commonedgeidx{chosen_idx})~= true)
        quad_edge1 = TO_Edges_VALID(commonedgeidx{chosen_idx}, :);
        plot(quad_edge1(:,1), quad_edge1(:,2), 'y.', 'MarkerSize',3.5);
        set(gcf,'color','w');
% % % %         if(paired_edge1(1,vi_test+2) >0)
% % % %             support_edge = collect_third_order_edges{validation_view_indices(vi_test),1}(paired_edge1(1,vi_test+2), :);
% % % %             plot(support_edge(1, 1), support_edge(1, 2), 'rx', 'MarkerSize',7, 'LineWidth', 1);
% % % %             plot(support_edge(1, 1), support_edge(1, 2), 'ro', 'MarkerSize',7, 'LineWidth', 1);
% % % %             %         plot(support_edge(1, 1), support_edge(1, 2), 'y.', 'LineWidth', 1);
% % % %             title(['view ', num2str(validation_view_indices(vi_test)), '(vali view ', num2str(vi_test), ') support'])
% % % %         end
        % plot(quad_edge(:,1), quad_edge(:,2), 'yo','MarkerSize',2);
        %     set(gcf,'color','w');
    end

   if (maxanglediff <=30 && minanglediff <=30 )
       title(['parallel-discard'])
%        title(['view ', num2str(validation_view_indices(vi_test)), '(vali view ', num2str(vi_test), ') parallel'])
       fprintf('View %d is having parallel epipolar wedges.\n',validation_view_indices(vi_test));
       fprintf('Max angle difference(in degree): %f\n',maxanglediff);
       fprintf('Min angle difference(in degree): %f\n',minanglediff);
   end
% % % %    if(chosen_idx == finalPairIndx)
% % % %        sgtitle ('Wrong Pair')
% % % %    else
% % % %        sgtitle ('Correct Pair')
% % % %    end
end