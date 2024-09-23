load('MATData/TO_Edges_syn.mat');
load('MATData/imageArray_syn.mat')
% R_matrixbeforeBA = load('MATData/R_matrix_10.mat').R_matrix;
% T_matrixbeforeBA = load('MATData/T_matrix_10.mat').T_matrix;
R_matrixafterBA  = load('MATData/R_matrix_syn.mat').R_matrix;
T_matrixafterBA  = load('MATData/T_matrix_syn.mat').T_matrix;
load('MATData/K_syn.mat')
% set all the parameters
params.NUM_OF_IMGS              = size(collect_third_order_edges, 1);
params.HYPO1_VIEW_INDX          = 4;
params.HYPO2_VIEW_INDX          = 8;
params.SUPPORT_DIST_THRESH      = 2;
params.SUPPORT_OREN_THRESH      = 0.99;
params.MAX_NUM_OF_SUPPORT_VIEWS = 2;
params.SHOW_SINGLE_VALIDATION   = 1;
params.SHOW_ALL_VALIDATIONS     = 1;
params.PERCENTAGE_FOR_EPIPOLE   = 0.005*5;
params.ANGLE_FOR_EPIPOLE        = 0.005;
params.ICL_DATA                 = 0;
params.multiK                   = 0;
params.cols                     = size(double(rgb2gray(imageArray{params.HYPO2_VIEW_INDX})), 2);
%> Convert from orientation to tangent vector for each edgel
TO_Edges_Tangents = cell(params.NUM_OF_IMGS, 1);
for i = 1:params.NUM_OF_IMGS
    TO_Edgels_Orientations = collect_third_order_edges{i,1}(:,3);
    
    Tgts = zeros(size(TO_Edgels_Orientations, 1), 2);
    for j = 1:size(TO_Edgels_Orientations, 1)
        Tgts(j,:) = [cos(TO_Edgels_Orientations(j,1)), ...
                     sin(TO_Edgels_Orientations(j,1))];
    end
    TO_Edges_Tangents{i,1} = Tgts;
end


load('GenData\gamma1_for_reproj_syn.mat')
recons_coor2 = Gamma1s';
% [pair2_sort,indexpair2] = sortrows(pair2,1);
% compare2pair=(pair1(:,2) == pair2_sort(:,2));
% recons_coor2 = recons_coor2(indexpair2,:);
% comparison12 = (recons_coor1 == recons_coor2);
% diff_idx = find(comparison12(:,1) == 0);
% load('GenData\gamma1_for_reproj_tgt_007.mat')
% GammaTangent_3D1 = GammaTangent_3D;
% load('GenData\gamma1_for_reproj_hnh_tgt_007.mat')
% GammaTangent_3D2 = GammaTangent_3D;
fignum = 1;
subnum = 1;
e3 = [0;0;1];
% abs_R1beforeBA =  R_matrixbeforeBA(:,:,params.HYPO1_VIEW_INDX);
% abs_C1beforeBA = -R_matrixbeforeBA(:,:,params.HYPO1_VIEW_INDX)' *...
%               T_matrixbeforeBA(:,params.HYPO1_VIEW_INDX);
abs_R1afterBA  =  R_matrixafterBA(:,:,params.HYPO1_VIEW_INDX);
abs_C1afterBA  = -R_matrixafterBA(:,:,params.HYPO1_VIEW_INDX)' *...
              T_matrixafterBA(:,params.HYPO1_VIEW_INDX);

reproject_error1_all = [];
reproject_error2_all = [];
validateidex = [3,7,10];
rangeplot = [1 2 3 5 6 8 9 10 11 12 13 14 15];
% for(vi = 1:3)
for(vidx = 1:3)
    vi = validateidex(1,vidx);
    reproject_coor1  = [];
    reproject_coor2  = [];
    reproject_tgt1   = [];
    reproject_tgt2   = [];
    reproject_error1 = [];    
    reproject_error2 = [];
    VALID_INDX = vi;%validation_view_indices(vi, 1);
    VALID_IMG  = double(rgb2gray(imageArray{VALID_INDX}));
    TO_Edges_VALID = collect_third_order_edges{VALID_INDX,1};
    Tangents_VALID = TO_Edges_Tangents{VALID_INDX,1};
    Kmatrix    = K;
    invK       = inv(Kmatrix);
%     abs_R2beforeBA =  R_matrixbeforeBA(:,:,VALID_INDX);
%     abs_C2beforeBA = -R_matrixbeforeBA(:,:,VALID_INDX)' *...
%                   T_matrixbeforeBA(:,VALID_INDX);
%     R21beforeBA    = abs_R2beforeBA * abs_R1beforeBA';
%     T21beforeBA    = abs_R2beforeBA * (abs_C1beforeBA - abs_C2beforeBA);
% 
%     extrbeforeBA           = zeros(3,4);
%     extrbeforeBA(1:3, 1:3) = R21beforeBA;
%     extrbeforeBA(1:3, 4)   = T21beforeBA;
%     PmatrixbeforeBA        = Kmatrix*extrbeforeBA;
% 
%     
%     for(i = 1 : size(recons_coor1,1))
%         reprojection = PmatrixbeforeBA*[recons_coor1(i,1);recons_coor1(i,2);recons_coor1(i,3);1];
%         x_coor       = reprojection(1,1)/reprojection(3,1);
%         y_coor       = reprojection(2,1)/reprojection(3,1);
%         if( x_coor > 0 && y_coor> 0 && x_coor <= size(imageArray{1},2) && y_coor <= size(imageArray{1},1))
%             reproject_coor1 = [reproject_coor1; [x_coor, y_coor]];
%             distance_error1 = ((TO_Edges_VALID(:,1)-x_coor).^2+(TO_Edges_VALID(:,2)-y_coor).^2).^(0.5);
%             
% %             T_v3 = R21 * GammaTangent_3D1(:,i);
% %             t_v3 = T_v3 - (e3'*T_v3)*invK*[reproject_coor1(1,:),1]';
% %             t_v3_reproj = t_v3 ./ norm(t_v3);
% %             reproject_tgt1  = [reproject_tgt1; t_v3_reproj(1,1), t_v3_reproj(2,1)];
% %             while(1)
% %                 [val_dis, index_dis] = mink(distance_error1,1);
% %                 target_edges      = Tangents_VALID(index_dis, :);
% %                 abs_dot_prod      = abs(t_v3_reproj(1,1).*target_edges(:,1) + t_v3_reproj(2,1).*target_edges(:,2));
% %                 if abs_dot_prod > params.SUPPORT_OREN_THRESH
%                     reproject_error1 = [reproject_error1; min(distance_error1)];
% %                     break;
% %                 else
% %                     distance_error1(index_dis,1) = inf;
% %                 end
% %             end
%         end
%     end

    abs_R2afterBA =  R_matrixafterBA(:,:,VALID_INDX);
    abs_C2afterBA = -R_matrixafterBA(:,:,VALID_INDX)' *...
                 T_matrixafterBA(:,VALID_INDX);
    R21afterBA    = abs_R2afterBA * abs_R1afterBA';
    T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);

    extrafterBA           = zeros(3,4);
    extrafterBA(1:3, 1:3) = R21afterBA;
    extrafterBA(1:3, 4)   = T21afterBA;
    PmatrixafterBA        = Kmatrix*extrafterBA;

    
%     for(i = 1 : size(recons_coor2,1))% index_end)%
     for(id = 1 : size(rangeplot,2))% index_end)%
         i = rangeplot(1,id);
        reprojection = PmatrixafterBA*[recons_coor2(i,1);recons_coor2(i,2);recons_coor2(i,3);1];
        x_coor       = reprojection(1,1)/reprojection(3,1);
        y_coor       = reprojection(2,1)/reprojection(3,1);
        if( x_coor > 0 && y_coor> 0 && x_coor <= size(imageArray{1},2) && y_coor <= size(imageArray{1},1))
            reproject_coor2 = [reproject_coor2; [x_coor, y_coor]];
            distance_error2 = ((TO_Edges_VALID(:,1)-x_coor).^2+(TO_Edges_VALID(:,2)-y_coor).^2).^(0.5);

%             T_v3 = R21 * GammaTangent_3D2(:,i);
%             t_v3 = T_v3 - (e3'*T_v3)*invK*[reproject_coor2(1,:),1]';
%             t_v3_reproj = t_v3 ./ norm(t_v3);
%             reproject_tgt2  = [reproject_tgt2; t_v3_reproj(1,1), t_v3_reproj(2,1)];
            reproject_error2 = [reproject_error2; min(distance_error2)];
%             while(1)
%                 [val_dis, index_dis] = mink(distance_error2,1);
%                 target_edges      = Tangents_VALID(index_dis, :);
%                 abs_dot_prod      = abs(t_v3_reproj(1,1).*target_edges(:,1) + t_v3_reproj(2,1).*target_edges(:,2));
%                 if abs_dot_prod > params.SUPPORT_OREN_THRESH
%                     reproject_error2 = [reproject_error2; min(distance_error2)];
%                     break;
%                 else
%                     distance_error2(index_dis,1) = inf;
%                 end
%             end

            
        end
    end
    
    %{}
figure(1);
     subplot(1,3,subnum)
%     subnum = subnum+1;
    imshow(imageArray{VALID_INDX});
    hold on;
    plot(TO_Edges_VALID(:,1), TO_Edges_VALID(:,2), 'y.', 'MarkerSize',3, 'LineWidth', 1);
%     plot(reproject_coor1(:,1), reproject_coor1(:,2), 'rx', 'MarkerSize',5, 'LineWidth', 1);
%     len = 2;
%     for(Indx = 1: size(reproject_coor1,1))
%         edgel_tgt2  = reproject_tgt1((Indx), 1:2);
%         hold on;
%         x1 = reproject_coor1(Indx, 1)+len*edgel_tgt2(1,1);
%         x2 = reproject_coor1(Indx, 1)-len*edgel_tgt2(1,1);
%         y1 = reproject_coor1(Indx, 2)+len*edgel_tgt2(1,2);
%         y2 = reproject_coor1(Indx, 2)-len*edgel_tgt2(1,2);
%         hold on;
%         plot([x1 x2], [y1 y2], 'r', 'LineWidth', 0.25);
%         hold on;
%     end
    plot(reproject_coor2(:,1), reproject_coor2(:,2), 'gx', 'MarkerSize',5, 'LineWidth', 1);
%     len = 2;
%     for(Indx = 1: size(reproject_coor2,1))
%         edgel_tgt2  = reproject_tgt2((Indx), 1:2);
%         hold on;
%         x1 = reproject_coor2(Indx, 1)+len*edgel_tgt2(1,1);
%         x2 = reproject_coor2(Indx, 1)-len*edgel_tgt2(1,1);
%         y1 = reproject_coor2(Indx, 2)+len*edgel_tgt2(1,2);
%         y2 = reproject_coor2(Indx, 2)-len*edgel_tgt2(1,2);
%         if(reproject_error2(Indx,1) >= 15)
%             hold on;
%             plot([x1 x2], [y1 y2], 'm', 'LineWidth', 0.25);
%             hold on;
%             plot(reproject_coor2(Indx,1), reproject_coor2(Indx,2), 'mx', 'MarkerSize',5, 'LineWidth', 2);
%         else
%             hold on;
%             plot([x1 x2], [y1 y2], 'g', 'LineWidth', 0.25);
%         end
%     end
    title(['view ', num2str(vi)])
    
%}    
    
    reproject_error1_all = [reproject_error1_all; reproject_error1];
    reproject_error2_all = [reproject_error2_all; reproject_error2];
%     figure(fignum);
%     subplot(2,3,subnum)
    subnum = subnum+1;
%     % histogram(reproject_error1,'BinWidth',1)
%     hold on;
%     histogram(reproject_error2,'BinWidth',1)
%     % legend('H2 only (95 pairs)', 'half and half (107 pairs)')
%     xlabel ('Reprojection Error', 'FontSize',10)
%     ylabel ('Times of occurrence', 'FontSize',10)
%     title(['view ', num2str(vi)], 'FontSize',10)
%     if(subnum > 3)
%         sgtitle ('Distribution of reprojection error for one image')
%         subnum = 1;
%         fignum = fignum+1;
%     end
end

figure;
histogram(reproject_error1_all,'BinWidth',1)
hold on;
histogram(reproject_error2_all,'BinWidth',1)
legend('1 round result before BA', '1 round result after BA')
xlabel ('Reprojection Error', 'FontSize',10)
ylabel ('Times of occurrence', 'FontSize',10)
title ('Distribution of reprojection error for all 50 images', 'FontSize',10)

figure;
subplot(1,2,1)
histogram(reproject_error1_all,'BinWidth',1)
%     ylim([0 2200])
%     xlim([-1 10])
xlabel ('Reprojection Error', 'FontSize',10)
ylabel ('Times of occurrence', 'FontSize',10)
title ('Distribution of reprojection error for all 50 images (1 round result before BA)', 'FontSize',10)
subplot(1,2,2)
histogram(reproject_error2_all,'BinWidth',1)
%     ylim([0 2200])
%     xlim([-1 10])
xlabel ('Reprojection Error', 'FontSize',10)
ylabel ('Times of occurrence', 'FontSize',10)
title ('Distribution of reprojection error for all 50 images (1 round result after BA)', 'FontSize',10)

idx11 = find(reproject_error1_all<=2);
idx22 = find(reproject_error2_all<=2);
ratio1 = size(idx11,1)/size(reproject_error1_all,1)
ratio2 = size(idx22,1)/size(reproject_error2_all,1)

% figure;
% plot3(recons_coor1(:,1), -recons_coor1(:,2), -recons_coor1(:,3), 'r.');
% hold on;
% plot3(recons_coor1(diff_idx,1), -recons_coor1(diff_idx,2), -recons_coor1(diff_idx,3), 'b.');
% axis equal;
% % title (['3D reconstruction result (1 round, wedge = ± ', num2str(params.ANGLE_FOR_EPIPOLE),'°)']) 
% title (['3D reconstruction result (single round)']) 
% 
% figure;
% plot3(recons_coor2(1:index_end,1,:), -recons_coor2(1:index_end,2), -recons_coor2(1:index_end,3), 'g.');
% hold on;
% plot3(recons_coor2(diff_idx,1), -recons_coor2(diff_idx,2), -recons_coor2(diff_idx,3), 'b.');
% axis equal;
% % title (['3D reconstruction result (accumulative, wedge = ± ', num2str(params.ANGLE_FOR_EPIPOLE),'°)']) 
% title (['3D reconstruction result (multiple round)']) 
rangeplot = 50:2:100;
figure;
plot3(recons_coor1(rangeplot,1), -recons_coor1(rangeplot,2), -recons_coor1(rangeplot,3), 'r.');
% hold on;
% plot3(recons_coor1(diff_idx,1), -recons_coor1(diff_idx,2), -recons_coor1(diff_idx,3), 'b.');
axis equal;
% title (['3D reconstruction result (1 round, wedge = ± ', num2str(params.ANGLE_FOR_EPIPOLE),'°)']) 
title (['3D reconstruction result (single round before BA)']) 

figure;
plot3(recons_coor2(rangeplot,1), -recons_coor2(rangeplot,2), -recons_coor2(rangeplot,3), 'g.');
% hold on;
% plot3(recons_coor2(diff_idx,1), -recons_coor2(diff_idx,2), -recons_coor2(diff_idx,3), 'b.');
axis equal;
% title (['3D reconstruction result (accumulative, wedge = ± ', num2str(params.ANGLE_FOR_EPIPOLE),'°)']) 
title (['3D reconstruction result (single round after BA)']) 