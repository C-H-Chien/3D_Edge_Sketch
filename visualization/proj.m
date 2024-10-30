edgesFilePath = '/gpfs/data/bkimia/zqiwu/3D_Edge_Sketch/build/bin/projected_vs_observed_edges.txt';
imagesFolderPath = '/gpfs/data/bkimia/Datasets/ABC-NEF/00000568/train_img/';

imageFiles = dir(fullfile(imagesFolderPath, '*_colors.png'));

fid = fopen(edgesFilePath, 'r');
if fid == -1
    error('Failed to open file: %s', edgesFilePath);
end

observedEdges = {};
projectedEdges = {};
claimedProjectedEdges = {};
frameIdx = 1;

currentSection = '';
while ~feof(fid)
    line = fgetl(fid);
    
    % Check for a frame header
    if contains(line, 'Frame')
        frameIdx = sscanf(line, 'Frame %d');
        observedEdges{frameIdx} = [];
        projectedEdges{frameIdx} = [];
        claimedProjectedEdges{frameIdx} = [];
    elseif contains(line, 'Projected Edges')
        currentSection = 'projected';
    elseif contains(line, 'Observed Edges')
        currentSection = 'observed';
    elseif ~isempty(line) && ~contains(line, 'Frame')
        % Read the edge data based on the current section
        edgeData = sscanf(line, '%f %f');
        if strcmp(currentSection, 'projected') && ~isempty(edgeData)
            projectedEdges{frameIdx} = [projectedEdges{frameIdx}; edgeData'];
        elseif strcmp(currentSection, 'observed') && ~isempty(edgeData)
            observedEdges{frameIdx} = [observedEdges{frameIdx}; edgeData'];
        end
    end
end

fclose(fid);

maxDistance = 5.0;
for i = 1:length(projectedEdges)
    if ~isempty(projectedEdges{i}) && ~isempty(observedEdges{i})
        for j = 1:size(projectedEdges{i}, 1)
            for k = 1:size(observedEdges{i}, 1)
                distance = sqrt((observedEdges{i}(k, 1) - projectedEdges{i}(j, 1))^2 + ...
                                (observedEdges{i}(k, 2) - projectedEdges{i}(j, 2))^2);

                % If the current projected edge is within maxDistance to an observed edge
                if distance < maxDistance
                    claimedProjectedEdges{i} = [claimedProjectedEdges{i}; projectedEdges{i}(j, :)];
                    break; % Move to the next projected edge as soon as a match is found
                end
            end
        end
    end
end

for i = 1:length(observedEdges)
    expectedImageName = sprintf('%d_colors.png', i - 1);
    imageIndex = find(strcmp({imageFiles.name}, expectedImageName));
    
    if ~isempty(imageIndex)
        img = imread(fullfile(imagesFolderPath, imageFiles(imageIndex).name));
        disp(imageFiles(imageIndex).name);
    else
        warning('No corresponding image found for Frame %d', i);
        img = []; % Empty image if not found
    end

    figure;
    if ~isempty(img)
        imshow(img);
        hold on;
    end
    
    title(sprintf('Frame %d: Projected vs Observed vs Claimed Projected Edges', i));
    
    if ~isempty(projectedEdges{i})
        plot(projectedEdges{i}(:, 1), projectedEdges{i}(:, 2), 'o', 'MarkerSize', 1, 'MarkerEdgeColor', 'b', 'DisplayName', 'Projected Edges');
    end
    
    if ~isempty(observedEdges{i})
        plot(observedEdges{i}(:, 1), observedEdges{i}(:, 2), 'o', 'MarkerSize', 1, 'MarkerEdgeColor', 'r', 'DisplayName', 'Observed Edges');
    end
    
    if ~isempty(claimedProjectedEdges{i})
        plot(claimedProjectedEdges{i}(:, 1), claimedProjectedEdges{i}(:, 2), 'o', 'MarkerSize', 1, 'MarkerEdgeColor', 'g', 'DisplayName', 'Claimed Projected Edges');
    end
    
    legend;
    xlabel('X');
    ylabel('Y');
    hold off;
    pause(0.5);
end
