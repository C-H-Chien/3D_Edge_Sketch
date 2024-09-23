% Epipolar wedges method is used in this main code.
clear all;
dataidx = ["0006";...
    "0077";...
    "0325";...
    "0568";...
    "2211"];

figure(1)
figure(2)
for(dataset = 1:5)
    if(dataset == 1)
        load('MATData/TO_Edges_ABC00000006.mat');
        load('MATData/imageArray_ABC00000006.mat')
        load('MATData/R_matrix_ABC00000006.mat')
        load('MATData/T_matrix_ABC00000006.mat')
        load('MATData/K_ABC00000006.mat')

        fileID = fopen('GenData/pairededge_ABC0006_6n8_t32to1excludehypo1n2.txt','r');
        % fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec);
        B = reshape(A,50,size(A,1)/50);
        B=B+1;
        paired_edge = B';
        H1_idx = 7;
        H2_idx = 9;
    elseif(dataset == 2)
        load('MATData/TO_Edges_ABC00000077.mat');
        load('MATData/imageArray_ABC00000077.mat')
        load('MATData/R_matrix_ABC00000077.mat')
        load('MATData/T_matrix_ABC00000077.mat')
        load('MATData/K_ABC00000077.mat')

        fileID = fopen('GenData/pairededge_ABC0077_4n23_t32to1excludehypo1n2.txt','r');
        % fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec);
        B = reshape(A,50,size(A,1)/50);
        B=B+1;
        paired_edge = B';
        H1_idx = 5;
        H2_idx = 24;
    elseif(dataset == 3)
        load('MATData/TO_Edges_ABC00000325.mat');
        load('MATData/imageArray_ABC00000325.mat')
        load('MATData/R_matrix_ABC00000325.mat')
        load('MATData/T_matrix_ABC00000325.mat')
        load('MATData/K_ABC00000325.mat')

        fileID = fopen('GenData/pairededge_ABC0325_3n11_t32to1excludehypo1n2.txt','r');
        % fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec);
        B = reshape(A,50,size(A,1)/50);
        B=B+1;
        paired_edge = B';
        H1_idx = 4;
        H2_idx = 12;
    elseif(dataset == 4)
        load('MATData/TO_Edges_ABC00000568.mat');
        load('MATData/imageArray_ABC00000568.mat')
        load('MATData/R_matrix_ABC00000568.mat')
        load('MATData/T_matrix_ABC00000568.mat')
        load('MATData/K_ABC00000568.mat')

        fileID = fopen('GenData/pairededge_ABC0568_21n39_t32to1excludehypo1n2.txt','r');
        % fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec);
        B = reshape(A,50,size(A,1)/50);
        B=B+1;
        paired_edge = B';
        H1_idx = 22;
        H2_idx = 40;
    elseif(dataset == 5)
        load('MATData/TO_Edges_ABC00002211.mat');
        load('MATData/imageArray_ABC00002211.mat')
        load('MATData/R_matrix_ABC00002211.mat')
        load('MATData/T_matrix_ABC00002211.mat')
        load('MATData/K_ABC00002211.mat')

        fileID = fopen('GenData/pairededge_ABC2211_4n18_t32to1excludehypo1n2.txt','r');
        % fileID = fopen('GenData/pairededge_ABC_7n9_t32to1excludehypo1n2.txt','r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec);
        B = reshape(A,50,size(A,1)/50);
        B=B+1;
        paired_edge = B';
        H1_idx = 5;
        H2_idx = 19;
    end

    x_img    = double(rgb2gray(imageArray{H1_idx}));
    x_xcor   = collect_third_order_edges{H1_idx,size(collect_third_order_edges,2)}(paired_edge(:,1),1);
    x_ycor   = collect_third_order_edges{H1_idx,size(collect_third_order_edges,2)}(paired_edge(:,1),2);
    x_oren   = collect_third_order_edges{H1_idx,size(collect_third_order_edges,2)}(paired_edge(:,1),3);
    x_tan    = [cos(x_oren), sin(x_oren)];
    x_tx     = [x_xcor+x_tan(:,1), x_xcor-x_tan(:,1)];
    x_ty     = [x_ycor+x_tan(:,2), x_ycor-x_tan(:,2)];
    x_oren   = x_oren + pi/2;
    x_cossin = [cos(x_oren), sin(x_oren)]*sqrt(2);
    x_x1     = x_xcor+x_cossin(:,1);
    x_y1     = x_ycor+x_cossin(:,2);
    x_x2     = x_xcor-x_cossin(:,1);
    x_y2     = x_ycor-x_cossin(:,2);
    x_co     = [];
    for(idxH1 = 1:size(x_oren,1))
        if(x_x1(idxH1,1)<x_xcor(idxH1,1))
            x_cocur = [abs(x_img(floor(x_y1(idxH1,1)),floor(x_x1(idxH1,1))) - ...
                x_img(floor(x_y2(idxH1,1)), floor(x_x2(idxH1,1)))); ...
                abs(x_img(ceil(x_y1(idxH1,1)),ceil(x_x1(idxH1,1))) - ...
                x_img(floor(x_y2(idxH1,1)), floor(x_x2(idxH1,1)))); ...
                abs(x_img(floor(x_y1(idxH1,1)),floor(x_x1(idxH1,1))) - ...
                x_img(ceil(x_y2(idxH1,1)), ceil(x_x2(idxH1,1)))); ...
                abs(x_img(ceil(x_y1(idxH1,1)),ceil(x_x1(idxH1,1))) - ...
                x_img(ceil(x_y2(idxH1,1)), ceil(x_x2(idxH1,1))))];
        else
            x_cocur = [abs(x_img(floor(x_y2(idxH1,1)),floor(x_x2(idxH1,1))) - ...
                x_img(floor(x_y1(idxH1,1)), floor(x_x1(idxH1,1)))); ...
                abs(x_img(ceil(x_y2(idxH1,1)),ceil(x_x2(idxH1,1))) - ...
                x_img(floor(x_y1(idxH1,1)), floor(x_x1(idxH1,1)))); ...
                abs(x_img(floor(x_y2(idxH1,1)),floor(x_x2(idxH1,1))) - ...
                x_img(ceil(x_y1(idxH1,1)), ceil(x_x1(idxH1,1)))); ...
                abs(x_img(ceil(x_y2(idxH1,1)),ceil(x_x2(idxH1,1))) - ...
                x_img(ceil(x_y1(idxH1,1)), ceil(x_x1(idxH1,1))))];
        end
        %     if(max(x_cocur)<= 5)
        %         figure;
        %         imshow(x_img,[]);
        %         hold on;
        %         plot(x_xcor(:,1),x_ycor(:,1), 'c.', 'MarkerSize',3);
        %         hold on;
        %         plot(x_xcor(idxH1,1),x_ycor(idxH1,1), 'mx', 'MarkerSize',7);
        %         hold on;
        %         plot(x_xcor(idxH1,1),x_ycor(idxH1,1), 'mx', 'MarkerSize',7);
        %         hold on;
        %         plot(x_x1(idxH1,1),x_y1(idxH1,1), 'm+', 'MarkerSize',7);
        %         hold on;
        %         plot(x_x2(idxH1,1),x_y2(idxH1,1), 'm+', 'MarkerSize',7);
        %         hold on;
        %         plot(x_tx(idxH1,:),x_ty(idxH1,:), 'm', 'LineWidth',1);
        %         hold off;
        %     end
        x_co = [x_co; max(x_cocur)];
    end




    y_img    = double(rgb2gray(imageArray{H2_idx}));
    y_xcor   = collect_third_order_edges{H2_idx,size(collect_third_order_edges,2)}(paired_edge(:,2),1);
    y_ycor   = collect_third_order_edges{H2_idx,size(collect_third_order_edges,2)}(paired_edge(:,2),2);
    y_oren   = collect_third_order_edges{H2_idx,size(collect_third_order_edges,2)}(paired_edge(:,2),3);
    y_tan    = [cos(x_oren), sin(x_oren)];
    y_tx     = [y_xcor+y_tan(:,1), y_xcor-y_tan(:,1)];
    y_ty     = [y_ycor+y_tan(:,2), y_ycor-y_tan(:,2)];
    y_oren   = y_oren + pi/2;
    y_cossin = [cos(y_oren), sin(y_oren)]*sqrt(2);
    y_x1     = y_xcor+y_cossin(:,1);
    y_y1     = y_ycor+y_cossin(:,2);
    y_x2     = y_xcor-y_cossin(:,1);
    y_y2     = y_ycor-y_cossin(:,2);
    y_co     = [];
    for(idxH2 = 1:size(y_oren,1))
        if(y_x1(idxH2,1)<y_xcor(idxH2,1))
            y_cocur = [abs(y_img(floor(y_y1(idxH2,1)),floor(y_x1(idxH2,1))) - ...
                y_img(floor(y_y2(idxH2,1)), floor(y_x2(idxH2,1)))); ...
                abs(y_img(ceil(y_y1(idxH2,1)),ceil(y_x1(idxH2,1))) - ...
                y_img(floor(y_y2(idxH2,1)), floor(y_x2(idxH2,1)))); ...
                abs(y_img(floor(y_y1(idxH2,1)),floor(y_x1(idxH2,1))) - ...
                y_img(ceil(y_y2(idxH2,1)), ceil(y_x2(idxH2,1)))); ...
                abs(y_img(ceil(y_y1(idxH2,1)),ceil(y_x1(idxH2,1))) - ...
                y_img(ceil(y_y2(idxH2,1)), ceil(y_x2(idxH2,1))))];
        else
            y_cocur = [abs(y_img(floor(y_y2(idxH2,1)),floor(y_x2(idxH2,1))) - ...
                y_img(floor(y_y1(idxH2,1)), floor(y_x1(idxH2,1)))); ...
                abs(y_img(ceil(y_y2(idxH2,1)),ceil(y_x2(idxH2,1))) - ...
                y_img(floor(y_y1(idxH2,1)), floor(y_x1(idxH2,1)))); ...
                abs(y_img(floor(y_y2(idxH2,1)),floor(y_x2(idxH2,1))) - ...
                y_img(ceil(y_y1(idxH2,1)), ceil(y_x1(idxH2,1)))); ...
                abs(y_img(ceil(y_y2(idxH2,1)),ceil(y_x2(idxH2,1))) - ...
                y_img(ceil(y_y1(idxH2,1)), ceil(y_x1(idxH2,1))))];
        end
        %     if(max(y_cocur)<= 5)
        %         figure;
        %         imshow(y_img,[]);
        %         hold on;
        %         plot(y_xcor(:,1),y_ycor(:,1), 'c.', 'MarkerSize',3);
        %         hold on;
        %         plot(y_xcor(idxH2,1),y_ycor(idxH2,1), 'mx', 'MarkerSize',7);
        %         hold on;
        %         plot(y_xcor(idxH2,1),y_ycor(idxH2,1), 'mx', 'MarkerSize',7);
        %         hold on;
        %         plot(y_x1(idxH2,1),y_y1(idxH2,1), 'm+', 'MarkerSize',7);
        %         hold on;
        %         plot(y_x2(idxH2,1),y_y2(idxH2,1), 'm+', 'MarkerSize',7);
        %         hold on;
        %         plot(y_tx(idxH2,:),y_ty(idxH2,:), 'm', 'LineWidth',1);
        %         hold off;
        %     end
        y_co = [y_co; max(y_cocur)];
    end

    x_data = x_co;
    y_data = y_co;


    % x_data = collect_third_order_edges{H1_idx,size(collect_third_order_edges,2)}(paired_edge(:,1),4);
    % y_data = collect_third_order_edges{H2_idx,size(collect_third_order_edges,2)}(paired_edge(:,2),4);



    figure(1)
    subplot(2,3,dataset)
    plot(x_data,y_data, 'bx', 'MarkerSize',7);
    xlabel 'contrast on the edge in H1'
    ylabel 'contrast on the edge in H2'
    title ('Object '+dataidx(dataset,1));
    xlim([0,255])
    ylim([0,255])
    figure(2)
    subplot(2,3,dataset)
    imshow(x_img,[]);
    title ('Object '+dataidx(dataset,1));
    figure(3)
    subplot(2,3,dataset)
    x_img_list = reshape(x_img,[size(x_img,1)*size(x_img,2),1]);
    histogram(x_img_list);
    title ('Object '+dataidx(dataset,1));
    xlim([0,260])
    figure(4)
    subplot(2,3,dataset)
    x_img_list = reshape(x_img,[size(x_img,1)*size(x_img,2),1]);
    histogram(x_co);
    title ('Object '+dataidx(dataset,1));
    xlim([0,260])
end
figure(1)
sgtitle 'contrast on the first edge vs. the contrast on the second edge'
figure(2)
sgtitle 'visualization of all 5 H1 images'
figure(3)
sgtitle 'distribution of intensity for all 5 H1 images'
figure(4)
sgtitle 'distribution of contrast for H1 images'