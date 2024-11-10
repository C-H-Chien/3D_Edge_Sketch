% MATLAB script to plot projected and observed edges
clc;
clear;
close all;

% Open and read the data file
filename = '/gpfs/data/bkimia/zqiwu/3D_Edge_Sketch/build/bin/projected_vs_observed_edges.txt';
fid = fopen(filename);

% Initialize frame counter
frame_counter = 0;

while ~feof(fid)
    % Read the frame header
    line = fgetl(fid);
    
    if startsWith(line, 'Frame')
        frame_counter = frame_counter + 1;
        fprintf('Processing %s\n', line);
        
        % Read projected edges
        fgetl(fid);  % Skip 'Projected Edges:' line
        projectedEdges = [];
        while true
            line = fgetl(fid);
            if isempty(line) || startsWith(line, 'Frame')
                break;
            end
            data = sscanf(line, '%f %f');
            projectedEdges = [projectedEdges; data'];
        end
        
        % Read observed edges
        fgetl(fid);  % Skip 'Observed Edges:' line
        observedEdges = [];
        while true
            line = fgetl(fid);
            if isempty(line) || feof(fid)
                break;
            end
            data = sscanf(line, '%f %f');
            observedEdges = [observedEdges; data'];
        end

        % Plot the edges
        figure;
        plot(projectedEdges(:,1), projectedEdges(:,2), 'bo', 'DisplayName', 'Projected Edges');
        hold on;
        plot(observedEdges(:,1), observedEdges(:,2), 'ro', 'DisplayName', 'Observed Edges');
        legend;
        title(['Frame ', num2str(frame_counter)]);
        xlabel('X');
        ylabel('Y');
    end
end

fclose(fid);
