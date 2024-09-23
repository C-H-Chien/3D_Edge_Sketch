function [ore_p1p2_all, ore_list2_sorted, epipole_pix_view1, epipole_pix_view2]= getOreList_New(edgel_HYPO1, edgels_HYPO2, R_matrix, T_matrix, V1IDX, V2IDX, params, K)

abs_R1 =  R_matrix(:,:,V1IDX);
abs_C1 = -R_matrix(:,:,V1IDX)' * T_matrix(:,V1IDX);
abs_R2 =  R_matrix(:,:,V2IDX);
abs_C2 = -R_matrix(:,:,V2IDX)' * T_matrix(:,V2IDX);

R21 = abs_R2 * abs_R1';
T21 = abs_R2 * (abs_C1 - abs_C2);

T_x = @(T)[0, -T(3,1), T(2,1); T(3,1), 0, -T(1,1); -T(2,1), T(1,1), 0];
E  = T_x(T21) * R21;

p1_all = [];
p2_all = [];

if(params.multiK == 1)
    K1 = K(:,:,V1IDX);
    K2 = K(:,:,V2IDX);
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

hypo1_ptall = edgel_HYPO1(:,1:2);

dxdy_hypo1all             = hypo1_ptall'-epipole_pix_view1(1:2,1);
dist_hypo1pt_epipoleall   = sqrt(dxdy_hypo1all(2,:).^2 + ... 
                                 dxdy_hypo1all(1,:).^2);
dist_hypo1pt12_epipoleall = sqrt(dist_hypo1pt_epipoleall.^2 - ... 
                                 (params.delta)^2);
thetahypo1all             = asin(params.delta./dist_hypo1pt_epipoleall);
anglehypo1all             = atan(dxdy_hypo1all(2,:)./dxdy_hypo1all(1,:));
angle_theta_hypo1all      = [anglehypo1all+thetahypo1all; ... 
                             anglehypo1all-thetahypo1all];
dxdyhypo1 = [(cos(angle_theta_hypo1all(1,:)).*dist_hypo1pt12_epipoleall)', ...
             (sin(angle_theta_hypo1all(1,:)).*dist_hypo1pt12_epipoleall)'];

dxdyhypo2 = [(cos(angle_theta_hypo1all(2,:)).*dist_hypo1pt12_epipoleall)', ...
             (sin(angle_theta_hypo1all(2,:)).*dist_hypo1pt12_epipoleall)'];

i_n1 = find(dxdy_hypo1all(1,:) < 0);
i_n2 = find(dxdy_hypo1all(2,:) < 0);

dxdyhypo1 = abs(dxdyhypo1);
dxdyhypo2 = abs(dxdyhypo2);

dxdyhypo1(i_n1',1) = -dxdyhypo1(i_n1',1);
dxdyhypo1(i_n2',2) = -dxdyhypo1(i_n2',2);
dxdyhypo2(i_n1',1) = -dxdyhypo2(i_n1',1);
dxdyhypo2(i_n2',2) = -dxdyhypo2(i_n2',2);

% dxdyhypo1 = -dxdyhypo1([hypo1_icn',hypo1_icp'],:);
% dxdyhypo1 = -dxdyhypo1([hypo1_icn',hypo1_icp'],:);
% dxdyhypo1 = -dxdyhypo1([hypo1_icn',hypo1_icp'],:);

% hypo1_icn = find(dxdyhypo1(i_n,1) >0);
% hypo1_icp = find(dxdyhypo1(i_p,1) <0);
% hypo2_icn = find(dxdyhypo2(i_n,1) >0);
% hypo2_icp = find(dxdyhypo2(i_p,1) <0);

% dxdyhypo1([hypo1_icn',hypo1_icp'],:) = -dxdyhypo1([hypo1_icn',hypo1_icp'],:);
% dxdyhypo2([hypo2_icn',hypo2_icp'],:) = -dxdyhypo2([hypo2_icn',hypo2_icp'],:);

hypo1_pt_1all             = dxdyhypo1 + epipole_pix_view1(1:2,1)';
hypo1_pt_2all             = dxdyhypo2 + epipole_pix_view1(1:2,1)';
hypo1_pt_all              = [hypo1_pt_1all, hypo1_pt_2all];

% for hypo 1
coeffspt1 = F * [hypo1_pt_all(:,1:2)'; ones(1,size(hypo1_pt_all,1))];
Apixel_1 = coeffspt1(1,:)';
Bpixel_1 = coeffspt1(2,:)';
Cpixel_1 = coeffspt1(3,:)';
coeffspt2 = F * [hypo1_pt_all(:,3:4)'; ones(1,size(hypo1_pt_all,1))];
Apixel_2 = coeffspt2(1,:)';
Bpixel_2 = coeffspt2(2,:)';
Cpixel_2 = coeffspt2(3,:)';

% for hypo 2
hypo2       = edgels_HYPO2(:,1:2)' - epipole_pix_view2(1:2,:);
slope_hypo2 = hypo2(2,:)./hypo2(1,:);
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

slope_hypo_pt1 = Bpixel_1./Apixel_1;
ore_list1bar_1 = slope_hypo_pt1;
ore_list1bar_1 = rad2deg(atan(ore_list1bar_1));
idx_less0_1    = find(ore_list1bar_1<0);
ore_list1bar_1(idx_less0_1,1) = ore_list1bar_1(idx_less0_1,1) + 180;

slope_hypo_pt2 = Bpixel_2./Apixel_2;
ore_list1bar_2 = slope_hypo_pt2;
ore_list1bar_2 = rad2deg(atan(ore_list1bar_2));
idx_less0_2    = find(ore_list1bar_2<0);
ore_list1bar_2(idx_less0_2,1) = ore_list1bar_2(idx_less0_2,1) + 180;

ore_list1bar_12 = [ore_list1bar_1, ore_list1bar_2];
[~, oremin_idx]    = min(ore_list1bar_12,[],2);
[~, oremax_idx]    = max(ore_list1bar_12,[],2);
ore_list1bar_all = ore_list1bar_12(:,[unique(oremin_idx), unique(oremax_idx)]);
% a_hypo2_all = tan(deg2rad(ore_list1bar_all));
% c_hypo2_all = epipole_pix_view2(2,1) - a_hypo2_all*epipole_pix_view2(1,1);
p1x_hypo_all = [zeros(size(Apixel_1,1),1), ...
                ones(size(Apixel_1,1),1)*params.cols, ...
                -Cpixel_1./Apixel_1, ...
                -(Cpixel_1 + Bpixel_1*params.rows)./Apixel_1];
p1y_hypo_all = [-Cpixel_1./Bpixel_1, ...
                -(Cpixel_1 + Apixel_1*params.cols)./Bpixel_1, ...
                zeros(size(Apixel_1,1),1), ...
                ones(size(Apixel_1,1),1)*params.rows];
p2x_hypo_all = [zeros(size(Apixel_2,1),1), ...
                ones(size(Apixel_2,1),1)*params.cols, ...
                -Cpixel_2./Apixel_2, ...
                -(Cpixel_2 + Bpixel_2*params.rows)./Apixel_2];
p2y_hypo_all = [-Cpixel_2./Bpixel_2, ...
                -(Cpixel_2 + Apixel_2*params.cols)./Bpixel_2, ...
                zeros(size(Apixel_2,1),1), ...
                ones(size(Apixel_2,1),1)*params.rows];
p1_idx_all = (1:size(edgel_HYPO1,1))'.*ones(size(edgel_HYPO1));
p1_1 = find(p1y_hypo_all(:,1) < 0 | p1y_hypo_all(:,1) > params.rows);
p1_2 = find(p1y_hypo_all(:,2) < 0 | p1y_hypo_all(:,2) > params.rows);
p1_3 = find(p1x_hypo_all(:,3) < 0 | p1x_hypo_all(:,3) > params.cols);
p1_4 = find(p1x_hypo_all(:,4) < 0 | p1x_hypo_all(:,4) > params.cols);
p1_idx_all(p1_1,1) = 0;
p1_idx_all(p1_2,2) = 0;
p1_idx_all(p1_3,3) = 0;
p1_idx_all(p1_4,4) = 0;
[p1_r,p1_c,~] = find(p1_idx_all);
if(size(p1_r,2) ~= 1)
    p1_r = p1_r';
    p1_c = p1_c';
end
p1_cr = [p1_r,p1_c];
[~, cr1_idx] = sort(p1_cr(:,1));
p1_cr = p1_cr(cr1_idx,:);
for(idxp1 = 1 : 2 : size(p1_cr,1))
p1_all = [p1_all;...
          p1x_hypo_all(p1_cr(idxp1,  1),p1_cr(idxp1,  2)),...
          p1y_hypo_all(p1_cr(idxp1,  1),p1_cr(idxp1,  2)), ...
          p1x_hypo_all(p1_cr(idxp1+1,1),p1_cr(idxp1+1,2)),...
          p1y_hypo_all(p1_cr(idxp1+1,1),p1_cr(idxp1+1,2))];
end
p2_idx_all = (1:size(edgel_HYPO1,1))'.*ones(size(edgel_HYPO1));
p2_1 = find(p2y_hypo_all(:,1) < 0 | p2y_hypo_all(:,1) > params.rows);
p2_2 = find(p2y_hypo_all(:,2) < 0 | p2y_hypo_all(:,2) > params.rows);
p2_3 = find(p2x_hypo_all(:,3) < 0 | p2x_hypo_all(:,3) > params.cols);
p2_4 = find(p2x_hypo_all(:,4) < 0 | p2x_hypo_all(:,4) > params.cols);
p2_idx_all(p2_1,1) = 0;
p2_idx_all(p2_2,2) = 0;
p2_idx_all(p2_3,3) = 0;
p2_idx_all(p2_4,4) = 0;
[p2_r,p2_c,~] = find(p2_idx_all);
if(size(p2_r,2) ~= 1)
    p2_r = p2_r';
    p2_c = p2_c';
end
p2_cr = [p2_r,p2_c];
[~, cr2_idx] = sort(p2_cr(:,1));
p2_cr = p2_cr(cr2_idx,:);
for(idxp2 = 1 : 2 : size(p2_cr,1))
p2_all = [p2_all;...
          p2x_hypo_all(p2_cr(idxp2,  1),p2_cr(idxp2,  2)),...
          p2y_hypo_all(p2_cr(idxp2,  1),p2_cr(idxp2,  2)), ...
          p2x_hypo_all(p2_cr(idxp2+1,1),p2_cr(idxp2+1,2)),...
          p2y_hypo_all(p2_cr(idxp2+1,1),p2_cr(idxp2+1,2))];
end
p1_dxdy(:,1:2) = p1_all(:,1:2) -  epipole_pix_view2(1:2,1)';
p1_dxdy(:,3:4) = p1_all(:,3:4) -  epipole_pix_view2(1:2,1)';
p1_dist(:,1)   = sqrt(p1_dxdy(:,1).^2 + p1_dxdy(:,2).^2);
p1_dist(:,2)   = sqrt(p1_dxdy(:,3).^2 + p1_dxdy(:,4).^2);
[~, p1_idx]    = min(p1_dist,[],2);
for(p1I = 1 : size(p1_idx,1))
    p1_xy(p1I,:) = p1_all(p1I,[p1_idx(p1I,1)*2-1, p1_idx(p1I,1)*2]);
end
% p1_xy(:,2) = p1_xy(:,2) - ...
%              abs(params.delta*cos(deg2rad(ore_list1bar_all(:,1))));

p2_dxdy(:,1:2) = p2_all(:,1:2) -  epipole_pix_view2(1:2,1)';
p2_dxdy(:,3:4) = p2_all(:,3:4) -  epipole_pix_view2(1:2,1)';
p2_dist(:,1)   = sqrt(p2_dxdy(:,1).^2 + p2_dxdy(:,2).^2);
p2_dist(:,2)   = sqrt(p2_dxdy(:,3).^2 + p2_dxdy(:,4).^2);
[~, p2_idx]    = min(p2_dist,[],2);
for(p2I = 1 : size(p2_idx,1))
    p2_xy(p2I,:) = p2_all(p2I,[p2_idx(p2I,1)*2-1, p2_idx(p2I,1)*2]);
end
% p2_xy(:,2) = p2_xy(:,2) + ...
%              abs(params.delta*cos(deg2rad(ore_list1bar_all(:,2))));
% [ore_list1_sorted(:,1),    ore_list1_sorted(:,2)   ] = sort(ore_list1);


% for hypo 2
p1_final_dxdyp      = p1_xy(:,1:2)' - epipole_pix_view2(1:2,:);
p1_final_dxdyn      = p1_final_dxdyp;
slope_p1p           = [];
slope_p1n           = [];
for(p1_final_idx = 1:size(p1_final_dxdyp,2))
    if(p1_xy(p1_final_idx,1) ==0 || ...
       p1_xy(p1_final_idx,1) == params.cols)
        p1_final_dxdyp(2,p1_final_idx) = p1_final_dxdyp(2,p1_final_idx) - ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p1_final_idx,1))))';
        slope_p1p(p1_final_idx,1) = p1_final_dxdyp(2,p1_final_idx)./...
                                    p1_final_dxdyp(1,p1_final_idx);
        p1_final_dxdyn(2,p1_final_idx) = p1_final_dxdyn(2,p1_final_idx) + ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p1_final_idx,1))))';
        slope_p1n(p1_final_idx,1) = p1_final_dxdyn(2,p1_final_idx)./...
                                    p1_final_dxdyn(1,p1_final_idx);
    else
        p1_final_dxdyp(1,p1_final_idx) = p1_final_dxdyp(1,p1_final_idx) - ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p1_final_idx,1))))';
        slope_p1p(p1_final_idx,1) = p1_final_dxdyp(2,p1_final_idx)./p1_final_dxdyp(1,p1_final_idx);
        p1_final_dxdyn(1,p1_final_idx) = p1_final_dxdyn(1,p1_final_idx) + ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p1_final_idx,1))))';
        slope_p1n(p1_final_idx,1) = p1_final_dxdyn(2,p1_final_idx)./p1_final_dxdyn(1,p1_final_idx);
    end
end
p2_final_dxdyp      = p2_xy(:,1:2)' - epipole_pix_view2(1:2,:);
p2_final_dxdyn      = p2_final_dxdyp;
slope_p2p           = [];
slope_p2n           = [];
for(p2_final_idx = 1:size(p2_final_dxdyp,2))
    if(p2_xy(p2_final_idx,1) ==0 || ...
       p2_xy(p2_final_idx,1) == params.cols)
        p2_final_dxdyp(2,p2_final_idx) = p2_final_dxdyp(2,p2_final_idx) - ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p2_final_idx,2))))';
        slope_p2p(p2_final_idx,1) = p2_final_dxdyp(2,p2_final_idx)./...
                                    p2_final_dxdyp(1,p2_final_idx);
        p2_final_dxdyn(2,p2_final_idx) = p2_final_dxdyn(2,p2_final_idx) + ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p2_final_idx,2))))';
        slope_p2n(p2_final_idx,1) = p2_final_dxdyn(2,p2_final_idx)./...
                                    p2_final_dxdyn(1,p2_final_idx);
    else
        p2_final_dxdyp(1,p2_final_idx) = p2_final_dxdyp(1,p2_final_idx) - ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p2_final_idx,2))))';
        slope_p2p(p2_final_idx,1) = p2_final_dxdyp(2,p2_final_idx)./p2_final_dxdyp(1,p2_final_idx);
        p2_final_dxdyn(1,p2_final_idx) = p2_final_dxdyn(1,p2_final_idx) + ...
                                         abs(params.delta*cos(deg2rad(...
                                         ore_list1bar_all(p2_final_idx,2))))';
        slope_p2n(p2_final_idx,1) = p2_final_dxdyn(2,p2_final_idx)./p2_final_dxdyn(1,p2_final_idx);
    end
end

ore_list_p1p = slope_p1p;
ore_list_p1p = rad2deg(atan(ore_list_p1p));
idx_less0   = find(ore_list_p1p<0);
ore_list_p1p(idx_less0,1) = ore_list_p1p(idx_less0,1) + 180;

ore_list_p1n = slope_p1n;
ore_list_p1n = rad2deg(atan(ore_list_p1n));
idx_less0   = find(ore_list_p1n<0);
ore_list_p1n(idx_less0,1) = ore_list_p1n(idx_less0,1) + 180;

ore_list_p2p = slope_p2p;
ore_list_p2p = rad2deg(atan(ore_list_p2p));
idx_less0   = find(ore_list_p2p<0);
ore_list_p2p(idx_less0,1) = ore_list_p2p(idx_less0,1) + 180;

ore_list_p2n = slope_p2n;
ore_list_p2n = rad2deg(atan(ore_list_p2n));
idx_less0   = find(ore_list_p2n<0);
ore_list_p2n(idx_less0,1) = ore_list_p2n(idx_less0,1) + 180;

ore_p1p2_4 = [ore_list_p1p, ore_list_p1n, ore_list_p2p, ore_list_p2n];

[~, p12_min]    = min(ore_p1p2_4,[],2);
[~, p12_max]    = max(ore_p1p2_4,[],2);
for(p12I = 1 : size(p1_idx,1))
    ore_p1p2_all(p12I,:) = [ore_p1p2_4(p12I,p12_min(p12I,1)), ...
                            ore_p1p2_4(p12I,p12_max(p12I,1)),];
end

[ore_list2_sorted(:,1),    ore_list2_sorted(:,2)   ] = sort(ore_list2);
% [ore_list1bar_sorted(:,1), ore_list1bar_sorted(:,2)] = sort(ore_list1bar);
end