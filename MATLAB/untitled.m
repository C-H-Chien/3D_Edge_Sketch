dist2Dall_all = [];
for(valiimg = 1:50)
% % % %     figure(fig)
% % % %     subplot(sbr,sbc,sub)
% % % % %     valididx = validation_view_indices(valiimg,1);
valididx = valiimg;
% % % %     imshow(imageArray{valididx});
% % % %     hold on;
% % % %     TO_Edges_1 = collect_third_order_edges{valididx,1};
% % % %     plot(TO_Edges_1(:,1), TO_Edges_1(:,2), 'c.', 'MarkerSize',4, 'LineWidth', 1);
    abs_R2afterBA =  R_matrixafterBA(:,:,valididx);
    abs_C2afterBA = -R_matrixafterBA(:,:,valididx)' *...
        T_matrixafterBA(:,valididx);
    R21afterBA    = abs_R2afterBA * abs_R1afterBA';
    T21afterBA    = abs_R2afterBA * (abs_C1afterBA - abs_C2afterBA);
    Kmatrix       = K(:,:,valididx);
    extrafterBA           = zeros(3,4);
    extrafterBA(1:3, 1:3) = R21afterBA;
    extrafterBA(1:3, 4)   = T21afterBA;
    PmatrixafterBA        = Kmatrix*extrafterBA;
    reprojection = PmatrixafterBA*[recons_coor1(:,1)';...
        recons_coor1(:,2)';...
        recons_coor1(:,3)';...
        ones(1,size(recons_coor1,1))];
    x_coor       = reprojection(1,:)./reprojection(3,:);
    y_coor       = reprojection(2,:)./reprojection(3,:);

    dx = ones(1,size(TO_Edges_1,1)).*x_coor' - TO_Edges_1(:,1)';
    dy = ones(1,size(TO_Edges_1,1)).*y_coor' - TO_Edges_1(:,2)';
    dist2Dall = sqrt(dx.^2 + dy.^2);
    [dist2D,dist2D_idx]=mink(dist2Dall,1,2);
    fp_pt  = find(dist2D(:,1) >= 2);
%     edge_left = dist2D_idx(fp_pt,:);
%     leftedge{valiimg} = TO_Edges_1(edge_left,:);
%     leftedgeidx{valiimg} = edge_left;

    tp_pt = 1:size(x_coor,2);
    tp_pt(fp_pt) = [];
    edge_ignore = dist2D_idx(tp_pt,:);
% % %     hold on;
% % %     plot(x_coor(1,tp_pt), y_coor(1,tp_pt), 'g.', 'MarkerSize',4, 'LineWidth', 1);
% % %     hold on;
% % %     plot(x_coor(1,fp_pt), y_coor(1,fp_pt), 'r.', 'MarkerSize',4, 'LineWidth', 1);
% % % %     legend ({'original edges' 'reprojected reconstruction edges tp' 'reprojected reconstruction edges fp'},...
% % % %         'FontSize',12);
% % %     title (['view ',num2str(valididx)])
    if(sub == 7 || sub == 10)
% % %         legend ({'original edges' 'reprojected reconstruction edges (TP)' 'reprojected reconstruction edges (FP)'},'FontSize',12);
    sub = sub;
    end
    sub = sub +1;
    if(sub>(sbr*sbc))
        sub = 1;
        fig = fig+1;
    end

% % % % % % % %     dist2Dall = dist2Dall';
% % % % % % % %     [dist2D1,dist2D_idx1]=mink(dist2Dall,1,2);
% % % % % % % %     edge_left  = find(dist2D1(:,1) >= 2);
%     edge_left = dist2D_idx(fp_pt,:);
%     leftedge{valiimg} = TO_Edges_1(edge_left,:);
%     leftedgeidx{valiimg} = edge_left;
    dist2Dall_all = [dist2Dall_all; dist2D];
end