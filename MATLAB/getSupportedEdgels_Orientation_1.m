function Supports = getSupportedEdgels_Orientation_1 ...
         (Tangents_VALID, reproj_edge_tgt_gamma3, params, inliner,chosen_idx)

    %> Code Description: 

    supported_edge_tgt_gamma3 = [];
    supported_pair_indx       = [];
    supported_link_indx       = [];
    
    %> Loop over all reprojected edgels and get supports from distance and orientation metrics
    all_oren = [];
    for ri = 1:size(reproj_edge_tgt_gamma3, 1)
        inliner_cell = inliner{ri};
        if(size(inliner_cell,2) ~= 0)

            %> Orientation test
            edgels_tgt_reproj = [reproj_edge_tgt_gamma3(ri,1); reproj_edge_tgt_gamma3(ri,2)];
            prev_prod         = 0;
            edge_tgt_gamma3   = [];
            pair_indx         = [];
            all_oren          = [];
            for(inlineridx = 1 : size(inliner_cell,1))
                target_edges      = Tangents_VALID(inliner_cell(inlineridx,1), :);
                abs_dot_prod      = abs(edgels_tgt_reproj(1,1).*target_edges(:,1) + edgels_tgt_reproj(2,1).*target_edges(:,2));
                all_oren = [all_oren; abs_dot_prod];
                if abs_dot_prod > params.SUPPORT_OREN_THRESH
                    if(abs_dot_prod > prev_prod)
                        prev_prod = abs_dot_prod;
                        edge_tgt_gamma3        = edgels_tgt_reproj';
                        pair_indx              = [ri, inliner_cell(inlineridx,1)];
                    end
                end
            end
            if(ri == chosen_idx)
                maxoren = max(all_oren)
            end
            if(size(pair_indx,2) ~= 0)
                supported_edge_tgt_gamma3 = [supported_edge_tgt_gamma3; edge_tgt_gamma3];
                supported_pair_indx       = [supported_pair_indx; pair_indx];
                supported_link_indx       = [supported_link_indx; pair_indx(1,2)];
            else
                supported_link_indx       = [supported_link_indx; 0];
            end
        else
            supported_link_indx       = [supported_link_indx; 0];
            continue;
        end
    end
    
    
    Supports.tangents  = supported_edge_tgt_gamma3;
    Supports.indices   = supported_pair_indx;
    Supports.link      = supported_link_indx;
end