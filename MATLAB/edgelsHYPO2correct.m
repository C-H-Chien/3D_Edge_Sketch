function [edgels_HYPO2_corrected, edgels_HYPO1_corrected, exclude_idx] = edgelsHYPO2correct(Epipolar_Coeffs, edgels_HYPO2, edgel_HYPO1, R_matrix, T_matrix, params, K)
edgels_HYPO1_corrected = zeros(size(edgels_HYPO2));
edgels_HYPO2_corrected = edgels_HYPO2;
exclude_idx = [];

abs_R1 =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C1 = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' *...
          T_matrix(:,params.HYPO2_VIEW_INDX);
abs_R2 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2 = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
          T_matrix(:,params.HYPO1_VIEW_INDX);
R21    = abs_R2 * abs_R1';
T21    = abs_R2 * (abs_C1 - abs_C2);
% Calculation 2
% R21 =  R_matrix(:,:,new_idx    ) * R_matrix(:,:,reference_idx)';
% T21 = -R_matrix(:,:,new_idx    ) * R_matrix(:,:,reference_idx)' * T_matrix(:,:,reference_idx) + T_matrix(:,:,new_idx    );
%> Calculate Essential matrix
T_x = @(T)[0,      -T(3,1),  T(2,1); ...
           T(3,1),  0,      -T(1,1); ...
          -T(2,1),  T(1,1),  0];
E   = T_x(T21) * R21;
% Calculate fundamental matrix
if(params.multiK == 1)
    K1 = K(:,:,params.HYPO2_VIEW_INDX);
    K2 = K(:,:,params.HYPO1_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
    F12   = invK2'*E*invK1;
else
    invK = inv(K);
    F12   = invK'*E*invK;
end

% match1 = [img1_points, ones(size(img1_points, 1), 1)]';
% match2 = [img2_points, ones(size(img1_points, 1), 1)]';

% InitialReprojectionErrors = zeros(size(edgels_HYPO2_corrected, 1), 1);
% OptimizedReprojectionErrors = zeros(size(edgels_HYPO2_corrected, 1), 1);

a1_line = -Epipolar_Coeffs.A./Epipolar_Coeffs.B;
b1_line = -1;
c1_line = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;

slop_edgeH1  = tan(edgel_HYPO1(1,3));
tgt_edgeH1   = [cos(edgel_HYPO1(1,3)), ...
                sin(edgel_HYPO1(1,3))];
a_edgeH1 = slop_edgeH1;
b_edgeH1 = -1;
c_edgeH1 = -(slop_edgeH1*edgel_HYPO1(1,1)-edgel_HYPO1(1,2));

for i = 1:size(edgels_HYPO2_corrected, 1)
    hypo2_current = edgels_HYPO2(i,:);
    if (params.syn == 1)
        slop_current  = hypo2_current(1,4)/hypo2_current(1,3);
    else
        slop_current  = tan(hypo2_current(1,3));
    end
    tgt_current   = [cos(hypo2_current(1,3)), ...
                     sin(hypo2_current(1,3))];
%     x1 = hypo2_current(1,1)+len*tgt_current(1,1);
%     x2 = hypo2_current(1,1)-len*tgt_current(1,1);
%     y1 = hypo2_current(1,2)+len*tgt_current(1,2);
%     y2 = hypo2_current(1,2)-len*tgt_current(1,2);

    a_current = slop_current;
    b_current = -1;
    c_current = -(slop_current*hypo2_current(1,1)-hypo2_current(1,2));

    correctedpt = [(b1_line*c_current-b_current*c1_line)/(a1_line*b_current-a_current*b1_line), ...
                   (c1_line*a_current-c_current*a1_line)/(a1_line*b_current-a_current*b1_line)];
    edgels_HYPO2_corrected_current = (correctedpt+hypo2_current(1,1:2))/2;

    % back to hypothesis view 1
    coeffs = F12 * [edgels_HYPO2_corrected_current(1,1:2)'; 1];
    a2_line = -coeffs(1,1)./coeffs(2,1);
    b2_line = -1;
    c2_line = -coeffs(3,1)./coeffs(2,1);
    correctedptH1 = [(b2_line*c_edgeH1-b_edgeH1*c2_line)/(a2_line*b_edgeH1-a_edgeH1*b2_line), ...
                     (c2_line*a_edgeH1-c_edgeH1*a2_line)/(a2_line*b_edgeH1-a_edgeH1*b2_line)];
    edgels_HYPO1_corrected_current = correctedptH1;

    edgels_HYPO2_corrected(i, 1:2) = edgels_HYPO2_corrected_current;
    edgels_HYPO1_corrected(i, :)   = [edgels_HYPO1_corrected_current, edgel_HYPO1(1,3:4)];
    dist_hypo2 = sqrt((edgels_HYPO2_corrected_current(1,1) - hypo2_current(1,1))^2 + ...
                      (edgels_HYPO2_corrected_current(1,2) - hypo2_current(1,2))^2);
    dist_hypo1 = sqrt((edgels_HYPO1_corrected_current(1,1) - edgel_HYPO1(1,1))^2 + ...
                      (edgels_HYPO1_corrected_current(1,2) - edgel_HYPO1(1,2))^2);
    if(dist_hypo1> params.circle || dist_hypo2> params.circle)
        exclude_idx = [exclude_idx; i];
    end
end

end