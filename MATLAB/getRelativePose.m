function [R21, T21, E, F] = getRelativePose ...
    (R_matrix, T_matrix, params, K)

%> Code Description:



%     R21 = transpose(R_matrix(:,:,params.HYPO2_VIEW_INDX)) *  R_matrix(:,:,params.HYPO1_VIEW_INDX);
%     T21 = transpose(R_matrix(:,:,params.HYPO2_VIEW_INDX)) * (T_matrix(:,:,params.HYPO1_VIEW_INDX) - T_matrix(:,:,params.HYPO2_VIEW_INDX));
abs_R1 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C1 = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
          T_matrix(:,params.HYPO1_VIEW_INDX);
abs_R2 =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C2 = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' *...
          T_matrix(:,params.HYPO2_VIEW_INDX);

% CPP modification
R21    = abs_R2 * abs_R1';
% CPP modification
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
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,params.HYPO2_VIEW_INDX);
    invK1 = inv(K1);
    invK2 = inv(K2);
    F   = invK2'*E*invK1;
else
    invK = inv(K);
    F   = invK'*E*invK;
end
end
