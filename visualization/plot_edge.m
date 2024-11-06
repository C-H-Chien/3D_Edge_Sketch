clear;
close all;

gamma1 = [520.554,  426.4240, 1]';
gamma2 = [519.6100, 407.3082, 1]';
k = [1111.11136542426,	0,	399.500000000000; 0,	1111.11136542426,	399.500000000000;  0,	0,	1];
k_inv = inv(k);


R_file = importdata('/gpfs/data/bkimia/Datasets/ABC-NEF/00000006/RnT/R_matrix.txt');
T_file = importdata('/gpfs/data/bkimia/Datasets/ABC-NEF/00000006/RnT/T_matrix.txt');

R_hyp1 = R_file(3*(48-1)+1:3*(48-1)+3, :);
T_hyp1 = T_file(3*(48-1)+1:3*(48-1)+3);

R_hyp2 = R_file(3*(43-1)+1:3*(43-1)+3, :);
T_hyp2 = T_file(3*(43-1)+1:3*(43-1)+3);

gamma_bar_1 = k_inv * gamma1;
gamma_bar_2 = k_inv * gamma2;

b3 = transpose([0, 0, 1]);

for i = 2:2
    index = i-1;

    img_path = sprintf('/gpfs/data/bkimia/Datasets/ABC-NEF/00000006/train_img/%d_colors.png', index);
    img = imread(img_path);
    figure;
    imshow(img);
    hold on;

    R_val = R_file(3*(i-1)+1:3*(i-1)+3, :);
    T_val = T_file(3*(i-1)+1:3*(i-1)+3);

    T21 = T_hyp2 - R_hyp2*R_hyp1'*T_hyp1;
    R21 = R_hyp2 * transpose(R_hyp1);

    T31 = T_val - R_val*R_hyp1'*T_hyp1;
    R31 = R_val * transpose(R_hyp1);
    
    top = [(b3'*T21) * (b3'*R21'*gamma_bar_2) - (b3'*R21'*T21)]*R31*gamma_bar_1 + [1-(b3'*R21*gamma_bar_1)*(b3'*R21'*gamma_bar_2)]*T31;
    bottom = [(b3'*T21)*(b3'*R21'*gamma_bar_2) - (b3'*R21'*T21)]*(b3'*R31*gamma_bar_1) + [1-(b3'*R21*gamma_bar_1)*(b3'*R21'*gamma_bar_2)]*(b3'*T31);
    gamma3 = top / bottom;

    gamma3_pixel = k * gamma3;
    plot(gamma3_pixel(1), gamma3_pixel(2), 'ro', 'MarkerSize', 4);
    hold off;

end

    

