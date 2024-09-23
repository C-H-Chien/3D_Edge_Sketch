function [edge_pos_gamma3, edge_tgt_gamma3] = getReprojectedEdgelCorrected ...
         (edgels_HYPO1_final, edgels_HYPO2, R_matrix, T_matrix, params, VALID_INDX, K)

    %> Code Description: 
    if(params.multiK == 1)
        K1 = K(:,:,params.HYPO1_VIEW_INDX);
        K2 = K(:,:,params.HYPO2_VIEW_INDX);
        K3 = K(:,:,VALID_INDX);
        invK1 = inv(K1);
        invK2 = inv(K2);
        invK3 = inv(K3);
    else
        invK = inv(K);
    end

    %>> Calculate relative pose (R, T) for reference & new, reference & support
    % for dataset "amsterdam_house_full", calculation 2 will be used
    % % Calculation 1
%     R21 = transpose(R_matrix(:,:,params.HYPO2_VIEW_INDX)) *  R_matrix(:,:,params.HYPO1_VIEW_INDX);
%     R31 = transpose(R_matrix(:,:,VALID_INDX))             *  R_matrix(:,:,params.HYPO1_VIEW_INDX);
%     T21 = transpose(R_matrix(:,:,params.HYPO2_VIEW_INDX)) * (T_matrix(:,:,params.HYPO1_VIEW_INDX) - T_matrix(:,:,params.HYPO2_VIEW_INDX));
%     T31 = transpose(R_matrix(:,:,VALID_INDX))             * (T_matrix(:,:,params.HYPO1_VIEW_INDX) - T_matrix(:,:,VALID_INDX));

    abs_R1 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
    abs_C1 = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' * T_matrix(:,params.HYPO1_VIEW_INDX);
    abs_R2 =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
    abs_C2 = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' * T_matrix(:,params.HYPO2_VIEW_INDX);
    abs_R3 =  R_matrix(:,:,VALID_INDX);
    abs_C3 = -R_matrix(:,:,VALID_INDX)' * T_matrix(:,VALID_INDX);

    R21 = abs_R2 * abs_R1';
    T21 = abs_R2 * (abs_C1 - abs_C2);
    R31 = abs_R3 * abs_R1';
    T31 = abs_R3 * (abs_C1 - abs_C3);

    e1 = [1;0;0];
    e3 = [0;0;1];

    edge_pos_gamma3   = zeros(size(edgels_HYPO2,1), 4);
    edge_tgt_gamma3   = zeros(size(edgels_HYPO2,1), 2);

%     point1            = [edgels_HYPO1_final(:,1:2)'; 1];
%     if(params.multiK == 1)
%         gamma1            = invK1*[edgels_HYPO1_final(:,1:2)'; 1];
%         tangent1          = [cos(edgels_HYPO1_final(:,3)); ...
%             sin(edgels_HYPO1_final(:,3))];
%         pt1_tgt_to_pixels = point1(1:2,1) + tangent1;
%         pt1_tgt_to_meters = invK1 * [pt1_tgt_to_pixels; 1];
%     else
%         gamma1            = invK*[edgels_HYPO1_final(:,1:2)'; 1];
%         tangent1          = [cos(edgels_HYPO1_final(:,3)); ...
%             sin(edgels_HYPO1_final(:,3))];
%         pt1_tgt_to_pixels = point1(1:2,1) + tangent1;
%         pt1_tgt_to_meters = invK * [pt1_tgt_to_pixels; 1];
%     end
%     tgt1_meters       = pt1_tgt_to_meters - gamma1;
%     edge_pos_gamma3   = zeros(size(edgels_HYPO2,1), 4);
%     edge_tgt_gamma3   = zeros(size(edgels_HYPO2,1), 2);

    for pi = 1 : size(edgels_HYPO2,1)
        point1            = [edgels_HYPO1_final(pi,1:2)'; 1];
        if(params.syn == 0)
            tangent1          = [cos(edgels_HYPO1_final(pi,3)); ...
                                 sin(edgels_HYPO1_final(pi,3))];
            tangent2          = [cos(edgels_HYPO2(pi,3)); ...
                                 sin(edgels_HYPO2(pi,3))];
        else
            tangent1          = [edgels_HYPO1_final(pi,3); ...
                                 edgels_HYPO1_final(pi,4)];
            tangent2          = [edgels_HYPO2(pi,3); ...
                                 edgels_HYPO2(pi,4)];
        end
        if(params.multiK == 1)
            gamma1            = invK1*[edgels_HYPO1_final(pi,1:2)'; 1];
            pt1_tgt_to_pixels = point1(1:2,1) + tangent1;
            pt1_tgt_to_meters = invK1 * [pt1_tgt_to_pixels; 1];
        else
            gamma1            = invK*[edgels_HYPO1_final(pi,1:2)'; 1];
            pt1_tgt_to_pixels = point1(1:2,1) + tangent1;
            pt1_tgt_to_meters = invK * [pt1_tgt_to_pixels; 1];
        end
        tgt1_meters       = pt1_tgt_to_meters - gamma1;
        
        
        %> fetch hypothesis pair in the second hypothesis view
        point2 = [edgels_HYPO2(pi, 1:2)'; 1];
        if(params.multiK == 1)
            gamma2 = invK2*[edgels_HYPO2(pi, 1:2)'; 1];
        else
            gamma2 = invK*[edgels_HYPO2(pi, 1:2)'; 1];
        end
        %> Get the position of the reprojected point
        rho1   = (e1'*T21 - (e3' * T21) * (e1'*gamma2)) / ...
                 ((e3'*R21*gamma1) * (e1'*gamma2) - (e1'*R21*gamma1));
        rho3   = rho1 * (e3'*R31*gamma1) + e3'*T31;
        gamma3 = 1 / rho3 * (R31*rho1*gamma1 + T31);
        if(params.multiK == 1)
            point3 = K3 * gamma3;
        else
            point3 = K * gamma3;
        end

        %> Analyticially find gamma3 from gamma1, gamma2, R12, T12, R13, and T13
        e1 = [1;0;0];
        e2 = [0;1;0];
        e3 = [0;0;1];
        numerator   = (e1'*T21 - (e3'*T21)*(e1'*gamma2))*R31*gamma1 + ((e3'*R21*gamma1)*(e1'*gamma2) - (e1'*R21*gamma1))*T31;
        denominator = (e3'*R31*gamma1)*(e1'*T21-(e3'*T21)*(e1'*gamma2)) + ((e3'*R21*gamma1)*(e1'*gamma2)-(e1'*R21*gamma1))*e3'*T31;
        gamma3 = numerator ./ denominator;
        if(params.multiK == 1)
            point3 = K3 * gamma3;
        else
            point3 = K * gamma3;
        end
        [~, ~, E31, ~] = getRelativePose1(R_matrix, T_matrix,params.HYPO1_VIEW_INDX,VALID_INDX,params, K);
        [~, ~, E32, ~] = getRelativePose1(R_matrix, T_matrix,params.HYPO2_VIEW_INDX,VALID_INDX,params, K);
        %> Alternate approach: analytically find gamma3 from the intersection of
        %  two epipolar lines
        a1 = E31(1,:)*gamma1;
        b1 = E31(2,:)*gamma1;
        c1 = E31(3,:)*gamma1;
        a2 = E32(1,:)*gamma2;
        b2 = E32(2,:)*gamma2;
        c2 = E32(3,:)*gamma2;
        bar_gamma3_x = -(b2*c1-b1*c2)/(b2*a1-b1*a2);
        bar_gamma3_y = -(a2*c1-a1*c2)/(a2*b1-a1*b2);
        
        if(params.multiK == 1)
            bar_point3 = K3 * [bar_gamma3_x; bar_gamma3_y; 1];
        else
            bar_point3 = K * [bar_gamma3_x; bar_gamma3_y; 1];
        end
%         point3     = bar_point3;
        



        edge_pos_gamma3(pi, :) = [point3(1:2, 1)', bar_point3(1:2, 1)'];
        
        %> Get the orientation of the reprojected point
        pt2_tgt_to_pixels = point2(1:2,1) + tangent2;
        if(params.multiK == 1)
            pt2_tgt_to_meters = invK2 * [pt2_tgt_to_pixels; 1];
        else
            pt2_tgt_to_meters = invK * [pt2_tgt_to_pixels; 1];
        end
        tgt2_meters       = pt2_tgt_to_meters - gamma2;

        P1   = gamma1;
        P2   = gamma2;
        t1   = tgt1_meters;
        t2   = tgt2_meters;
        n1   = cross(t1, P1);
        n2   = R21' * cross(t2, P2);
        T_v1 = cross(n1,n2) ./ norm(cross(n1,n2));

        T_v3 = R31 * T_v1;
        t_v3 = T_v3 - (e3'*T_v3)*gamma3;
        t_v3 = t_v3 ./ norm(t_v3);

        if(params.ICL_DATA == 1)
            t_v3(2,1) = -1*t_v3(2,1);
        end
        
        tangent3              = t_v3(1:2,1);
        edge_tgt_gamma3(pi,:) = tangent3';
    end


end