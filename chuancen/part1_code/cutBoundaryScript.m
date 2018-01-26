% Demo code for alignment of Prokudin-Gorskii images
%
% Your code should be run after you implement alignChannels

% Path to your data directory
dataDir = fullfile('..','output','prokudin-gorskii');

% Path to your output directory (change this to your output directory)
outDir = fullfile('..', 'output', 'alignment-prokudin-gorskii-cutedge');

% List of images
imageNames = {'00125-aligned.jpg',	'00153-aligned.jpg', '00398-aligned.jpg', '00149-aligned.jpg', '00351-aligned.jpg',	'01112-aligned.jpg'};

% Display variable
display = true;

% Set maximum shift to check alignment for
maxShift = [15 15];

% Loop over images, untile them into images, align
for i = 1:length(imageNames),
    
    % Read image
    im = imread(fullfile(dataDir, imageNames{i}));
    
    % Convert to double
    im = im2double(im);
    
    % Images are stacked vertically
    % From top to bottom are B, G, R channels (and not RGB)
    imageHeight = size(im,1);
    imageWidth  = size(im,2);
    
    % Allocate memory for the image 
%     channels = zeros(imageHeight, imageWidth, 3);
    
%     % We are loading the color channels from top to bottom
%     % Note the ordering of indices
%     channels(:,:,3) = im(1:imageHeight,:);
%     channels(:,:,2) = im(imageHeight+1:2*imageHeight,:);
%     channels(:,:,1) = im(2*imageHeight+1:3*imageHeight,:);
% 
%     % Align the blue and red channels to the green channel
%     [colorIm, predShift] = alignChannels(channels, maxShift);
%     
%     % Print the alignments
%     fprintf('%2i %s shift: B (%2i,%2i) R (%2i,%2i)\n', i, imageNames{i}, predShift(1,:), predShift(2,:));

    %im_gray = rgb2gray(im);
    
    [im_cut] = cut_boundary_effect(im);


    % Write image output
    outimageName = sprintf([imageNames{i}(1:end-4) '-cutted.jpg']);
    outimageName = fullfile(outDir, outimageName);
    imwrite(im_cut, outimageName);
   
    outimageName2 = sprintf([imageNames{i}(1:end-4) '-cutted-comp.jpg']);
    outimageName2 = fullfile(outDir, outimageName2);
    % Optionally display the results
    if display
        f=figure;
        %figure(1); clf;
        subplot(1,2,1); imagesc(im); axis image off; colormap gray
        title('Input image');
        subplot(1,2,2); imagesc(im_cut); axis image off;
        title('Cutted image');
        saveas(f,outimageName2,'jpg');
        pause(1.5);
    end
end







