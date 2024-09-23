load('MATData/TO_Edges_ABC00000006.mat');
load('MATData/imageArray_ABC00000006.mat');

% ImageMatri = zeros(400,400);
% 
% TO_Edges_HYPO1 = collect_third_order_edges{1,size(collect_third_order_edges,2)};
Arrayedge = [];
for (imgidx = 1: 50)
    ImageMatri = zeros(800,800);
    TO_Edges_HYPO1 = collect_third_order_edges{imgidx,size(collect_third_order_edges,2)};
    for edge_idx = 1:size(TO_Edges_HYPO1, 1)
        pt = TO_Edges_HYPO1(edge_idx,1:2);
        pt = floor(pt);
        ImageMatri(pt(1,1),pt(1,2)) = ImageMatri(pt(1,1),pt(1,2))+1;
    end
    Arrayedge_current = reshape(ImageMatri,800*800,1);
    Arrayedge_current = nonzeros(Arrayedge_current');
    Arrayedge = [Arrayedge; Arrayedge_current];
end
histogram(Arrayedge,'BinWidth', 1)
xlabel ('Number of edges per pixel', 'FontSize',20)
ylabel ('Number of pixels', 'FontSize',20)
% title ('Distribution of number of edges per pixel', 'FontSize',15)

hypo_img1 = double(rgb2gray(imageArray{7}));
TO_Edges_HYPO1  = collect_third_order_edges{7,size(collect_third_order_edges,2)};

figure
imshow(uint8(hypo_img1));
hold on;
plot(TO_Edges_HYPO1(:, 1), TO_Edges_HYPO1(:, 2), 'c.', 'MarkerSize',5, 'LineWidth', 2);