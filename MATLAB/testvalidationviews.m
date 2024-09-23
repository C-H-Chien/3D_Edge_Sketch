for vi = 1:size(validation_view_indices, 1)
    % get the validation view
    VALID_INDX = validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    % get edges in the validation view
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
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
    for(hypo2idx = 2:size(edgels_HYPO2_final,1))
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
    end
    %         if(isempty(allquad))
    %             indice = [];
    %             supported_indices_stack       = [supported_indices_stack; indice];
    %             supported_link                = [supported_link, zeros(size(commonedgeidx,2),1)];
    %         else
    Supports_orient = getSupportedEdgels_Orientation(Tangents_VALID, ...
        reproj_edge_tgt_gamma3, params, commonedgeidx, isparal);
    %         supported_edges_hypo2vali{vi} = Supports_orient.indices;
    supported_indices_stack  = [supported_indices_stack; Supports_orient.indices];
    supported_link          = [supported_link, Supports_orient.link];
    supported_simi          = [supported_simi, Supports_orient.simi];

    %         end
end