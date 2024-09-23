figure
% hypo2idx = finalPairIndx;
len = 0.5;
[p1_all, p2_all]= getOreList_New_green(TO_Edges_HYPO1(edge_idx,:), TO_Edges_HYPO2, ...
                                     R_matrix, T_matrix, ...
                                     params.HYPO1_VIEW_INDX, ...
                                     params.HYPO2_VIEW_INDX, params, K);
% edgels_HYPO2_final = edgels_HYPO2_corrected;
% edgels_HYPO1_final = edgels_HYPO1_corrected;
tgt1 = TO_Edges_Tangents{params.HYPO1_VIEW_INDX,1};
edgel_tgt1  = tgt1(edge_idx, 1:2);
tgt2 = TO_Edges_Tangents{params.HYPO2_VIEW_INDX,1};
subplot(1,2,1);
imshow(uint8(hypo_img1));
hold on;
idx_exclude = find(TO_Edges_HYPO1(:,1) >= 5 & TO_Edges_HYPO1(:,1) <= cols-5 & ...
            TO_Edges_HYPO1(:,2) >= 5 & TO_Edges_HYPO1(:,2) <= rows-5);
plot(TO_Edges_HYPO1(idx_exclude, 1), TO_Edges_HYPO1(idx_exclude, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
hold on;
x1 = edgel_HYPO1(1,1)+len*edgel_tgt1(1,1);
x2 = edgel_HYPO1(1,1)-len*edgel_tgt1(1,1);
y1 = edgel_HYPO1(1,2)+len*edgel_tgt1(1,2);
y2 = edgel_HYPO1(1,2)-len*edgel_tgt1(1,2);
hold on;
plot([x1 x2], [y1 y2], 'r', 'LineWidth', 1);
hold on;
hypo31_all = edgel_HYPO1(1,1:2)' - epipole_pix_view1(1:2,:);
b31_all    = hypo31_all(2,:)./hypo31_all(1,:);
c31_all    = edgel_HYPO1(1,2) - b31_all * edgel_HYPO1(1,1);
y31_min1 = b31_all*1 + c31_all;
y31_max1 = b31_all*cols + c31_all;
hold on;
abs_R1 =  R_matrix(:,:,params.HYPO2_VIEW_INDX);
abs_C1 = -R_matrix(:,:,params.HYPO2_VIEW_INDX)' *...
    T_matrix(:,params.HYPO2_VIEW_INDX);
abs_R2 =  R_matrix(:,:,params.HYPO1_VIEW_INDX);
abs_C2 = -R_matrix(:,:,params.HYPO1_VIEW_INDX)' *...
    T_matrix(:,params.HYPO1_VIEW_INDX);
R21    = abs_R2 * abs_R1';
T21    = abs_R2 * (abs_C1 - abs_C2);
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
k1 = (edgel_HYPO1(1,2)-epipole_pix_view1(2,1))/(edgel_HYPO1(1,1)-epipole_pix_view1(1,1));
b1 = epipole_pix_view1(2,1) - k1*epipole_pix_view1(1,1);
yMax1 = k1*cols + b1;
hold on;
line([0, cols], [b1, yMax1], 'Color', 'g', 'LineWidth', 0.5);
% for(idxH1 = 1: size (edgels_HYPO1_final,1))
%     coeffs = F12 * [edgels_HYPO2_final(idxH1,1:2)'; 1];
%     yMin = -coeffs(3,1)./coeffs(2,1);
%     yMax = (-coeffs(3,1) - coeffs(1,1)*cols) ./ coeffs(2,1);
%     line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 0.5);
%     hold on;
%     plot(edgels_HYPO1_final(idxH1, 1), edgels_HYPO1_final(idxH1, 2), 'gx', 'MarkerSize',7, 'LineWidth', 1);
% end
hold on;
plot(edgel_HYPO1(1, 1), edgel_HYPO1(1, 2), 'rx', 'MarkerSize',10, 'LineWidth', 1);
% distancecorrectedH1 = sqrt((edgels_HYPO1_final(idxH1,1)-edgel_HYPO1(1,1))^2+(edgels_HYPO1_final(idxH1,2)-edgel_HYPO1(1,2))^2);
title 'hypothesis view 1'
legend '2D edges from this view' 'γ1' 'γ1 orientation' 'epipolar line'
% legend '2D edges from this view' 'γ1' 'γ1 orientation' '' 'corrected γ1' 'epipolar line of corrected γ2'

subplot(1,2,2);
imshow(uint8(hypo_img2));
title 'hypothesis view 2'
hold on;
idx_exclude2 = find(TO_Edges_HYPO2(:,1) >= 5 & TO_Edges_HYPO2(:,1) <= cols-5 & ...
            TO_Edges_HYPO2(:,2) >= 5 & TO_Edges_HYPO2(:,2) <= rows-5);
plot(TO_Edges_HYPO2(idx_exclude2, 1), TO_Edges_HYPO2(idx_exclude2, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);
hold on;
yMin = -Epipolar_Coeffs.C./Epipolar_Coeffs.B;
yMax = (-Epipolar_Coeffs.C - Epipolar_Coeffs.A*cols) ./ Epipolar_Coeffs.B;
line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 0.5,'LineStyle','-');
% y = k1*x + b1
k1 = tan(thresh_ore1/180*pi);
b1 = epipole_pix_view2(2,1) - k1*epipole_pix_view2(1,1);
yMax1 = k1*cols + b1;
% y = k2*x + b2
k2 = tan(thresh_ore2/180*pi);
b2 = epipole_pix_view2(2,1) - k2*epipole_pix_view2(1,1);
yMax2 = k2*cols + b2;
hold on;
line([0, cols], [b1, yMax1], 'Color', 'm', 'LineWidth', 0.5);
hold on;
line([0, cols], [b2, yMax2], 'Color', 'm', 'LineWidth', 0.5);
% hold on; 
% line([p1_all(1,1), p1_all(1,3)], [p1_all(1,2), p1_all(1,4)], 'Color', 'g', 'LineWidth', 0.5);
% hold on; 
% line([p2_all(1,1), p2_all(1,3)], [p2_all(1,2), p2_all(1,4)], 'Color', 'g', 'LineWidth', 0.5);
% plot(edgels_HYPO2(:, 1), edgels_HYPO2(:, 2), 'mx', 'MarkerSize',7, 'LineWidth', 1);
len = 0.5;
for(Indx = 1: size(edgels_HYPO2,1))
    %     Indx = 1;
    edgel_tgt2  = tgt2(hypo2_idxsorted(Indx), 1:2);
    if edgels_HYPO2(Indx,1) < 5 || edgels_HYPO2(Indx,1) > cols-5 || ...
            edgels_HYPO2(Indx,2) < 5 || edgels_HYPO2(Indx,2) > rows-5
        continue;
    end
    % hold on;
    % x1 = edgels_HYPO2_final(Indx, 1)+len*edgel_tgt2(1,1);
    % x2 = edgels_HYPO2_final(Indx, 1)-len*edgel_tgt2(1,1);
    % y1 = edgels_HYPO2_final(Indx, 2)+len*edgel_tgt2(1,2);
    % y2 = edgels_HYPO2_final(Indx, 2)-len*edgel_tgt2(1,2);
    % hold on;
    % plot([x1 x2], [y1 y2], 'm', 'LineWidth', 0.25);
    hold on;
    plot(edgels_HYPO2(Indx, 1), edgels_HYPO2(Indx, 2), 'rx', 'MarkerSize',7, 'LineWidth', 1);
    x1 = edgels_HYPO2(Indx, 1)+len*edgel_tgt2(1,1);
    x2 = edgels_HYPO2(Indx, 1)-len*edgel_tgt2(1,1);
    y1 = edgels_HYPO2(Indx, 2)+len*edgel_tgt2(1,2);
    y2 = edgels_HYPO2(Indx, 2)-len*edgel_tgt2(1,2);
    hold on;
    plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.5);
%     hold on;
%     plot(edgels_HYPO2_final(Indx, 1), edgels_HYPO2_final(Indx, 2), 'rx', 'MarkerSize',7, 'LineWidth', 1);
%     x1 = edgels_HYPO2_final(Indx, 1)+len*edgel_tgt2(1,1);
%     x2 = edgels_HYPO2_final(Indx, 1)-len*edgel_tgt2(1,1);
%     y1 = edgels_HYPO2_final(Indx, 2)+len*edgel_tgt2(1,2);
%     y2 = edgels_HYPO2_final(Indx, 2)-len*edgel_tgt2(1,2);
%     hold on;
%     plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.5);
%     if(Indx == 1)
%         hold on; plot(edgels_HYPO2_final(hypo2idx, 1), edgels_HYPO2_final(hypo2idx, 2), 'go', 'MarkerSize',13, 'LineWidth', 1);
%     end
%     distancecorrectedH2 = sqrt((edgels_HYPO2_final(Indx,1)-edgels_HYPO2(Indx,1))^2+(edgels_HYPO2_final(Indx,2)-edgels_HYPO2(Indx,2))^2);
end
% hold on; 
% line([p1_all(1,1), p1_all(1,3)], [p1_all(1,2), p1_all(1,4)], 'Color', 'g', 'LineWidth', 0.5);
% hold on; 
% line([p2_all(1,1), p2_all(1,3)], [p2_all(1,2), p2_all(1,4)], 'Color', 'g', 'LineWidth', 0.5);

legend '2D edges from this view' 'epipolar line of γ1' '' '' 'γ2' 'γ2 orientation'
% legend '2D edges from this view' 'epipolar line of γ1' 'epipolar wedge(purple)' '' 'epipolar wedge (green)' '' 'γ2' 'γ2 orientation'
% legend '2D edges from this view' 'epipolar line of γ1' 'epipolar wedge(purple)' '' 'epipolar wedge (green)' '' 'γ2' 'γ2 orientation' 'corrected γ2''corrected γ2 orientation' 'corrected γ2 chosen to visualize on vali views'
% ptdelta = [121.6434  188.6648  122.1446  188.3348];
% coeffs1 = F * [ptdelta(1,1:2)'; 1];
%     Apixel1 = coeffs1(1,1);
%     Bpixel1 = coeffs1(2,1);
%     Cpixel1 = coeffs1(3,1);
%     Epipolar_Coeffs1.A = Apixel1;
%     Epipolar_Coeffs1.B = Bpixel1;
%     Epipolar_Coeffs1.C = Cpixel1;
% hold on;
% yMin = -Epipolar_Coeffs1.C./Epipolar_Coeffs1.B;
% yMax = (-Epipolar_Coeffs1.C - Epipolar_Coeffs1.A*cols) ./ Epipolar_Coeffs1.B;
% line([0, cols], [yMin, yMax], 'Color', 'g', 'LineWidth', 0.5,'LineStyle','--');