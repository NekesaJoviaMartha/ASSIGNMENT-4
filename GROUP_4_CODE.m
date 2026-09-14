   %% LEAF STRUCTURAL ARRAY ASSIGNMENT

clear; 
clc; 
close all;

inputFolder = 'C:\Users\DELL\Desktop\GROUP_4_PICKED_LEAVES';
outputFolder = fullfile(fileparts(mfilename('fullpath')), 'GROUP_4_leaf_analysis_results');
if ~isfolder(outputFolder), mkdir(outputFolder);
end

imageFiles = dir(fullfile(inputFolder, 'leaf *.jpeg'));

% Natural-sort the names: leaf 1, leaf 2, ..., leaf 19.
numbers = cellfun(@(x) sscanf(x, 'leaf %d.jpeg'), {imageFiles.name});
[~, order] = sort(numbers);
imageFiles = imageFiles(order);

leafData = struct([]);

for k = 1:numel(imageFiles)
    fileName = imageFiles(k).name;
    leafNumber = sscanf(fileName, 'leaf %d.jpeg');
    rgb = imread(fullfile(inputFolder, fileName));
    rgb = im2uint8(rgb);

    % Resize only for analysis: preserves the original image in leafData.
    maxDimension = 1200;
    scale = min(1, maxDimension / max(size(rgb,1), size(rgb,2)));
    analysisRGB = imresize(rgb, scale);
    [mask, hsvImage, labImage] = makeLeafMask(analysisRGB);

    % Processed versions required by the assignment.
    grayImage = rgb2gray(analysisRGB);
    enhancedGray = adapthisteq(grayImage);
    edgeImage = edge(enhancedGray, 'Canny');
    maskedRGB = analysisRGB;
    maskedRGB(repmat(~mask, 1, 1, 3)) = 0;

    stats = regionprops(mask, 'Area', 'Perimeter', 'Centroid', ...
        'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', ...
        'Orientation', 'BoundingBox', 'Extent', 'Solidity', 'EquivDiameter');
    stats = stats(1);

    % Measurements are in pixels because no physical scale was photographed.
    props = struct();
    props.Area_pixels = stats.Area;
    props.Perimeter_pixels = stats.Perimeter;
    props.Centroid_xy_pixels = stats.Centroid;
    props.MajorAxis_pixels = stats.MajorAxisLength;
    props.MinorAxis_pixels = stats.MinorAxisLength;
    props.AspectRatio = stats.MajorAxisLength / max(stats.MinorAxisLength, eps);
    props.Eccentricity = stats.Eccentricity;
    props.Orientation_degrees = stats.Orientation;
    props.Extent = stats.Extent;
    props.Solidity = stats.Solidity;
    props.EquivalentDiameter_pixels = stats.EquivDiameter;
    props.Roundness = 4*pi*stats.Area / max(stats.Perimeter^2, eps);
    props.BoundingBox_xywh_pixels = stats.BoundingBox;

    R = analysisRGB(:,:,1); G = analysisRGB(:,:,2); B = analysisRGB(:,:,3);
    props.MeanRGB = [mean(R(mask)), mean(G(mask)), mean(B(mask))];
    props.MeanHSV = [mean(hsvImage(:,:,1), 'all', 'omitnan'), ...
                     mean(hsvImage(:,:,2), 'all', 'omitnan'), ...
                     mean(hsvImage(:,:,3), 'all', 'omitnan')];
    % Recalculate HSV means inside mask (the previous values are images).
    h = hsvImage(:,:,1); s = hsvImage(:,:,2); v = hsvImage(:,:,3);
    props.MeanHSV = [mean(h(mask)), mean(s(mask)), mean(v(mask))];
    props.MeanLab = [mean(labImage(:,:,1), 'all'), mean(labImage(:,:,2), 'all'), mean(labImage(:,:,3), 'all')];
    L = labImage(:,:,1); a = labImage(:,:,2); b = labImage(:,:,3);
    props.MeanLab = [mean(L(mask)), mean(a(mask)), mean(b(mask))];
    props.ExcessGreen = mean(2*double(G(mask)) - double(R(mask)) - double(B(mask)));
    props.Texture = getTextureFeatures(enhancedGray, mask);
    props.MorphologyDescription = describeLeaf(props);

    leafData(k).LeafID = leafNumber;
    leafData(k).FileName = fileName;
    leafData(k).SourcePath = fullfile(inputFolder, fileName);
    leafData(k).OriginalImage = rgb;
    leafData(k).AnalysisImage = analysisRGB;
    leafData(k).GrayImage = grayImage;
    leafData(k).EnhancedGrayImage = enhancedGray;
    leafData(k).HSVImage = hsvImage;
    leafData(k).LabImage = labImage;
    leafData(k).BinaryLeafMask = mask;
    leafData(k).EdgeImage = edgeImage;
    leafData(k).MaskedLeafImage = maskedRGB;
    leafData(k).Properties = props;
    % Save viewable outputs for each photograph.
    imwrite(mask, fullfile(outputFolder, sprintf('leaf_%02d_mask.png', leafNumber)));
    imwrite(maskedRGB, fullfile(outputFolder, sprintf('leaf_%02d_segmented.png', leafNumber)));
    imwrite(edgeImage, fullfile(outputFolder, sprintf('leaf_%02d_edges.png', leafNumber)));

    fig = figure('Visible','off','Color','w','Position',[100 100 1200 700]);
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    nexttile; imshow(analysisRGB); title(sprintf('Leaf %d: original',leafNumber));
    nexttile; imshow(mask); title('Binary leaf mask');
    nexttile; imshow(maskedRGB); title('Segmented leaf');
    nexttile; imshow(enhancedGray); title('Contrast-enhanced grayscale');
    nexttile; imshow(edgeImage); title('Canny edges');
    nexttile; imshow(analysisRGB); hold on;
    visboundaries(mask,'Color','y','LineWidth',1);
    plot(stats.Centroid(1),stats.Centroid(2),'r+','MarkerSize',12,'LineWidth',2);
    title(sprintf('%s | AR %.2f | Roundness %.2f', props.MorphologyDescription, props.AspectRatio, props.Roundness));
    exportgraphics(fig, fullfile(outputFolder, sprintf('leaf_%02d_report.png',leafNumber)), 'Resolution', 160);
    close(fig);
end

% Save the requested structural array in a MATLAB data file.
save(fullfile(outputFolder, 'leafData_structural_array.mat'), 'leafData', '-v7.3');

% Make a compact table that can be opened in MATLAB or spreadsheet software.
summaryTable = struct2table(arrayfun(@(x) flattenLeaf(x), leafData));
writetable(summaryTable, fullfile(outputFolder, 'leaf_measurements.csv'));
save(fullfile(outputFolder, 'leafData_summary_table.mat'), 'summaryTable');

fprintf('Finished %d leaf images. Results saved in:\n%s\n', numel(leafData), outputFolder);
disp(summaryTable);

%% Local functions
function [mask, hsvImage, labImage] = makeLeafMask(rgb)
% Separates the leaf from the brown wooden background. It combines green
% information, colour distance from the image-border background, and cleanup.
    hsvImage = rgb2hsv(rgb);
    labImage = rgb2lab(rgb);
    [rows, cols, ~] = size(rgb);
    border = false(rows, cols);
    border([1:round(.08*rows), end-round(.08*rows)+1:end], :) = true;
    border(:, [1:round(.08*cols), end-round(.08*cols)+1:end]) = true;

    L = labImage(:,:,1); a = labImage(:,:,2); b = labImage(:,:,3);
    bgLab = [median(L(border)), median(a(border)), median(b(border))];
    colourDistance = sqrt((L-bgLab(1)).^2 + (a-bgLab(2)).^2 + (b-bgLab(3)).^2);
    distanceMask = colourDistance > max(9, prctile(colourDistance(:), 68));

    R = im2double(rgb(:,:,1)); G = im2double(rgb(:,:,2)); B = im2double(rgb(:,:,3));
    excessGreen = 2*G - R - B;
    greenMask = excessGreen > max(0.025, graythresh(excessGreen)*0.35);
    hue = hsvImage(:,:,1); sat = hsvImage(:,:,2);
    vegetationMask = hue > 0.16 & hue < 0.52 & sat > 0.12;
    candidate = (distanceMask & (greenMask | vegetationMask)) | (vegetationMask & sat > .2);

    candidate = imclose(candidate, strel('disk', 8));
    candidate = imopen(candidate, strel('disk', 2));
    candidate = imfill(candidate, 'holes');
    candidate = bwareafilt(candidate, 1); % the photographed leaf is largest object
    candidate = imclearborder(candidate); % remove any background touching border
    if ~any(candidate(:))
        candidate = bwareafilt(distanceMask, 1); % safe fallback for non-green leaves
    end
    mask = imfill(imclose(candidate, strel('disk', 5)), 'holes');
end

function texture = getTextureFeatures(grayImage, mask)
    crop = grayImage;
    crop(~mask) = 0;
    glcm = graycomatrix(crop, 'Offset',[0 1; -1 1; -1 0; -1 -1], 'Symmetric',true);
    g = graycoprops(glcm, {'Contrast','Correlation','Energy','Homogeneity'});
    texture = struct('Contrast',mean(g.Contrast), 'Correlation',mean(g.Correlation), ...
        'Energy',mean(g.Energy), 'Homogeneity',mean(g.Homogeneity), ...
        'GrayStd',std(double(grayImage(mask))));
end

function description = describeLeaf(p)
    if p.Solidity < 0.78
        description = 'lobed or compound';
    elseif p.AspectRatio > 4
        description = 'very narrow / linear';
    elseif p.AspectRatio > 2.4
        description = 'elongated / lanceolate';
    elseif p.Roundness > 0.72
        description = 'nearly round / broad';
    else
        description = 'elliptic or ovate';
    end
end

function row = flattenLeaf(x)
    p = x.Properties; t = p.Texture;
    row = struct('LeafID',x.LeafID,'FileName',x.FileName,'Morphology',p.MorphologyDescription, ...
        'Area_pixels',p.Area_pixels,'Perimeter_pixels',p.Perimeter_pixels, ...
        'AspectRatio',p.AspectRatio,'Roundness',p.Roundness,'Solidity',p.Solidity, ...
        'Eccentricity',p.Eccentricity,'MeanRed',p.MeanRGB(1),'MeanGreen',p.MeanRGB(2), ...
        'MeanBlue',p.MeanRGB(3),'ExcessGreen',p.ExcessGreen, ...
        'TextureContrast',t.Contrast,'TextureEnergy',t.Energy,'TextureHomogeneity',t.Homogeneity);
end

