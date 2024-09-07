close all;
clear all;

%> This code aims for loading .cemv file
data_path = "/home/chchien/datasets/ABC_NEF/00000006/train_img/";
file_name = "4_colors";
cemv_postfix = ".cemv";
fig_postfix = ".png";

fid = fopen(strcat(data_path, file_name, cemv_postfix));
fid_fig_str = strcat(data_path, file_name, fig_postfix);
if fid < 0, error('CEMV file not found.'); end
% if fid_fig < 0, error('Picture not found'); end

%> Cell array storing contour information
cem = cell(3,1);

cf_idx = cell(1,0);
contours = cell(1,0);
num_contours = 1;

%scan the file to read the contour information
while ~feof(fid)
    %> Read the line of the file
    lineBuffer = fgetl(fid);
    % if ~ischar(lineBuffer), error, end

    %> Ignore comment lines and empty lines
    if (length(lineBuffer) < 2 || lineBuffer(1)=='#')
        continue;
    end

    if strncmp(lineBuffer, 'CONTOUR_COUNT=', length('CONTOUR_COUNT=')) || ...
       strncmp(lineBuffer, 'TOTAL_EDGE_COUNT=', length('TOTAL_EDGE_COUNT='))
        lineBuffer = fgetl(fid);
        continue;
    end

    %
    if strncmp(lineBuffer, '[BEGIN CONTOUR]', length('[BEGIN CONTOUR]'))
        lineBuffer = fgetl(fid);
        num_of_edges = strread(lineBuffer, 'EDGE_COUNT=%d');
        contours{num_contours} = zeros(num_of_edges, 2);
        continue;
    end

    for i = 1:num_of_edges
        ldata = textscan(lineBuffer, '[%d, %d]\t%f\t%f\t[%f, %f]\t%f\t%f', 'CollectOutput', true);
        edges = ldata{1,2};
        contours{num_contours}(i,:) = edges(3:4);
        lineBuffer = fgetl(fid);
    end

    if strncmp(lineBuffer, '[END CONTOUR]', length('[BEGIN CONTOUR]'))
        num_contours = num_contours + 1;
    else
        error("Somehting's wrong!\n");
    end
end

sz = [num_contours 3];
color_code = unifrnd(0,1,sz);
img = imread(fid_fig_str);
figure(1);
imshow(img); hold on;
for ci = 1:num_contours-1
    plot(contours{ci}(:,1), contours{ci}(:,2), 'Color', color_code(ci,:), 'LineWidth', 5);
    hold on;
end

%> Close the file
fclose(fid);
