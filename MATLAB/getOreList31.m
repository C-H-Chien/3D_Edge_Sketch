function [ore_list1bar, ore_list2_sorted, epipole_pix_view1, epipole_pix_view2]= getOreList31(edgels_HYPO1, edgels_HYPO2, R_matrix, T_matrix, params, K, idx2)

abs_R1 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C1 = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' * T_matrix(:,params.HYPO1_VIEW_INDX);
abs_R2 =  R_matrix(:,:,idx2);
abs_C2 = -R_matrix(:,:,idx2)' * T_matrix(:,idx2);

R21 = abs_R2 * abs_R1';
T21 = abs_R2 * (abs_C1 - abs_C2);

T_x = @(T)[0, -T(3,1), T(2,1); T(3,1), 0, -T(1,1); -T(2,1), T(1,1), 0];
E  = T_x(T21) * R21;

if(params.multiK == 1)
    K1 = K(:,:,params.HYPO1_VIEW_INDX);
    K2 = K(:,:,idx2);
    invK1 = inv(K1);
    invK2 = inv(K2);
    F   = invK2'*E*invK1;
else
    invK = inv(K);
    F   = invK'*E*invK;
end

e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% Calculation for epipole
epipole_met_view1 = [(e1'*R21'*T21) / (e3'*R21'*T21); (e2'*R21'*T21) / (e3'*R21'*T21)];
epipole_met_view2 = [(e1'*T21) / (e3'*T21); (e2'*T21) / (e3'*T21)];
if(params.multiK == 1)
    epipole_pix_view1 = K1 * [epipole_met_view1; 1];
    epipole_pix_view2 = K2 * [epipole_met_view2; 1];
else
    epipole_pix_view1 = K * [epipole_met_view1; 1];
    epipole_pix_view2 = K * [epipole_met_view2; 1];
end

% for hypo 1
hypo1       = edgels_HYPO1(:,1:2)' - epipole_pix_view1(1:2,:);
slope_hypo1 = hypo1(2,:)./hypo1(1,:);
% for hypo 2
hypo2       = edgels_HYPO2(:,1:2)' - epipole_pix_view2(1:2,:);
slope_hypo2 = hypo2(2,:)./hypo2(1,:);

ore_list1    = slope_hypo1';
ore_list2    = slope_hypo2';
ore_list2    = rad2deg(atan(ore_list2));
idx_less0 = find(ore_list2<0);
ore_list2(idx_less0,1) = ore_list2(idx_less0,1) + 180;
% for(idx_hypo2 = 1:size(hypo2,2))
%     if(hypo2(1,idx_hypo2)>=0 && hypo2(2,idx_hypo2)>=0)
%         ore_list2(idx_hypo2,1)    = rad2deg(atan(hypo2(2,idx_hypo2)./hypo2(1,idx_hypo2)));
%     elseif(hypo2(1,idx_hypo2)<0 && hypo2(2,idx_hypo2)>=0)
%         ore_list2(idx_hypo2,1)    = rad2deg(atan(hypo2(2,idx_hypo2)./hypo2(1,idx_hypo2)))+180;
%     elseif(hypo2(1,idx_hypo2)<0 && hypo2(2,idx_hypo2)<0)
%         ore_list2(idx_hypo2,1)    = rad2deg(atan(hypo2(2,idx_hypo2)./hypo2(1,idx_hypo2)))+180;
%     else
%         ore_list2(idx_hypo2,1)    = rad2deg(atan(hypo2(2,idx_hypo2)./hypo2(1,idx_hypo2)))+360;
%     end
% end
%%%
ore_list1bar = -(F(1,1)+F(1,2)*ore_list1)./(F(2,1)+F(2,2)*ore_list1);
ore_list1bar    = rad2deg(atan(ore_list1bar));
idx_less0_1 = find(ore_list1bar<0);
ore_list1bar(idx_less0_1,1) = ore_list1bar(idx_less0_1,1) + 180;

[ore_list1_sorted(:,1),    ore_list1_sorted(:,2)   ] = sort(ore_list1);
[ore_list2_sorted(:,1),    ore_list2_sorted(:,2)   ] = sort(ore_list2);
[ore_list1bar_sorted(:,1), ore_list1bar_sorted(:,2)] = sort(ore_list1bar);
end