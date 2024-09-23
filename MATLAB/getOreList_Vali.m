function [ore_31bar, ore_list31_sorted, epipole_pix_view1, epipole_pix_view3]= getOreList_Vali(edgel_HYPO1, TO_Edges_VALID, R_matrix, T_matrix, HYPO1_VIEW_INDX, K, VALID_INDX,params)

abs_R1 =  R_matrix(:,:,HYPO1_VIEW_INDX);
abs_C1 = -R_matrix(:,:,HYPO1_VIEW_INDX)' * T_matrix(:,HYPO1_VIEW_INDX);
abs_R3 =  R_matrix(:,:,VALID_INDX);
abs_C3 = -R_matrix(:,:,VALID_INDX)' * T_matrix(:,VALID_INDX);

R31 = abs_R3 * abs_R1';
T31 = abs_R3 * (abs_C1 - abs_C3);

T_x = @(T)[0, -T(3,1), T(2,1); T(3,1), 0, -T(1,1); -T(2,1), T(1,1), 0];
E  = T_x(T31) * R31;
if(params.multiK == 1)
    K1 = K(:,:,HYPO1_VIEW_INDX);
    K3 = K(:,:,VALID_INDX);
    invK1 = inv(K1);
    invK3 = inv(K3);
    F   = invK3'*E*invK1;
else
    invK = inv(K);
    F   = invK'*E*invK;
end

e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% Calculation for epipole
epipole_met_view1 = [(e1'*R31'*T31) / (e3'*R31'*T31); (e2'*R31'*T31) / (e3'*R31'*T31)];
epipole_met_view3 = [(e1'*T31) / (e3'*T31); (e2'*T31) / (e3'*T31)];
if(params.multiK == 1)
    epipole_pix_view1 = K1 * [epipole_met_view1; 1];
    epipole_pix_view3 = K3 * [epipole_met_view3; 1];
else
    epipole_pix_view1 = K * [epipole_met_view1; 1];
    epipole_pix_view3 = K * [epipole_met_view3; 1];
end

% for hypo 1
hypo1       = edgel_HYPO1(:,1:2)' - epipole_pix_view1(1:2,:);
slope_hypo1 = hypo1(2,:)./hypo1(1,:);
% for validation
hypo3       = TO_Edges_VALID(:,1:2)' - epipole_pix_view3(1:2,:);
slope_hypo3 = hypo3(2,:)./hypo3(1,:);

ore_31    = slope_hypo1';
ore_list31    = slope_hypo3';
ore_list31    = rad2deg(atan(ore_list31));
idx_less0 = find(ore_list31<0);
ore_list31(idx_less0,1) = ore_list31(idx_less0,1) + 180;
% for(idx_hypo3 = 1:size(hypo3,2))
%     if(hypo3(1,idx_hypo3)>=0 && hypo3(2,idx_hypo3)>=0)
%         ore_list31(idx_hypo3,1)    = rad2deg(atan(hypo3(2,idx_hypo3)./hypo3(1,idx_hypo3)));
%     elseif(hypo3(1,idx_hypo3)<0 && hypo3(2,idx_hypo3)>=0)
%         ore_list31(idx_hypo3,1)    = rad2deg(atan(hypo3(2,idx_hypo3)./hypo3(1,idx_hypo3)))+180;
%     elseif(hypo3(1,idx_hypo3)<0 && hypo3(2,idx_hypo3)<0)
%         ore_list31(idx_hypo3,1)    = rad2deg(atan(hypo3(2,idx_hypo3)./hypo3(1,idx_hypo3)))+180;
%     else
%         ore_list31(idx_hypo3,1)    = rad2deg(atan(hypo3(2,idx_hypo3)./hypo3(1,idx_hypo3)))+360;
%     end
% end
ore_31bar = -(F(1,1)+F(1,2)*ore_31)./(F(2,1)+F(2,2)*ore_31);
ore_31bar    = rad2deg(atan(ore_31bar));
idx_less0_1 = find(ore_31bar<0);
ore_31bar(idx_less0_1,1) = ore_31bar(idx_less0_1,1) + 180;

[ore_list31_sorted(:,1),    ore_list31_sorted(:,2)   ] = sort(ore_list31);
end