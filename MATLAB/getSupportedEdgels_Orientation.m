function Supports = getSupportedEdgels_Orientation ...
         (Tangents_VALID, reproj_edge_tgt_gamma3, params, inliner,isparal)

    %> Code Description: 

    abs_dot_prod = 0;
    supported_edge_tgt_gamma3 = [];
    supported_pair_indx       = [];
    supported_link_indx       = [];
    supported_pair_simi       = [];
    supported_link_curr       = [];
    
    %> Loop over all reprojected edgels and get supports from distance and orientation metrics
    all_oren = [];
    for ri = 1:size(reproj_edge_tgt_gamma3, 1)
        if(isparal(ri,1) == 0)
            supported_link_indx       = [supported_link_indx; {0}];
            continue;
        end
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
                if abs_dot_prod >= abs(cos(deg2rad(params.SUPPORT_OREN_THRESH)))
%                     if(abs_dot_prod >= prev_prod)
%                         prev_prod              = abs_dot_prod;
%                         edge_tgt_gamma3        = edgels_tgt_reproj';
%                         pair_indx              = [ri, inliner_cell(inlineridx,1)];
%                     end
                        prev_prod              = abs_dot_prod;
                        edge_tgt_gamma3        = edgels_tgt_reproj';
                        pair_indx              = [pair_indx, inliner_cell(inlineridx,1)];
                end
            end
            if(size(pair_indx,2) ~= 0)
                supported_edge_tgt_gamma3 = [supported_edge_tgt_gamma3; edge_tgt_gamma3];
                supported_pair_indx       = [supported_pair_indx; ri, pair_indx(1,1)];
                supported_link_indx       = [supported_link_indx; {pair_indx}];
                supported_pair_simi       = [supported_pair_simi; prev_prod];
            else
                supported_link_indx       = [supported_link_indx; {0}];
            end
        else
            supported_link_indx       = [supported_link_indx; {0}];
            continue;
        end
    end
    
    
    Supports.tangents  = supported_edge_tgt_gamma3;
    Supports.indices   = supported_pair_indx;
    Supports.link      = supported_link_indx;
    Supports.simi      = supported_pair_simi;
end