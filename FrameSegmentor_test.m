% FrameSegmentor_test - This is a frame segmentor example for identifying
% cell texture in brightfield imaging and splitting into components using a
% number of nuclear fluorescent channel markers.
% 
% hulls = FrameSegmentor(chanIm, primaryChan, t, seRadius)
% INPUTS: 
%   chanIm - A cell array each chanIm{i} contains the image intensity date 
%   for the i-th channel at frame t. Some cells may be empty if they were 
%   unavailable or not imaged at frame t.
% 
%   primaryChan - A number between 1 and CONSTANTS.numChannels indicating
%   the primary channel for the segmentation. Specific algorithms may use
%   information from other available channels as well.
% 
%   t - The frame number of the image data being passed into the
%   segmentation algorithm.
% 
%   seRadius - Radius of a neaighborhood element for the brightfield
%   texture filter. Increasing the radius will generally connect
%   segmentations.
% 
%
% OUTPUTS:
%   hulls - The hulls output should be a structure array with one entry per
%   segmentation result.
%   At a minimum, each hull must contain an 'indexPixels' field, a list of
%   linear indices representing the interior pixels of segmentation results.
%   For example, the PixelIdxList result from a call to regionprops could
%   be used to fill the hulls indexPixels field.
% 
function hulls = FrameSegmentor_test(chanIm, primaryChan, t, seRadius)
    hulls = [];
    
    otherChan = setdiff(1:length(chanIm),primaryChan);
    %% Brightfield texture segmentation
    bfIm = chanIm{primaryChan};
    
    se = strel('disk',seRadius);
    stdIm = entropyfilt(mat2gray(bfIm),se.getnhood);
    
    hFilt = fspecial('gaussian', 6,2);
    
    smoothStd = imfilter(stdIm,hFilt);
    boundLoss = imfilter(ones(size(stdIm)),hFilt);
    boundLoss(boundLoss == 0) = 1.0;
    
    smoothStd = smoothStd ./ boundLoss;
    
    stdImNrm = mat2gray(smoothStd);
    stdThresh = graythresh(stdImNrm);
    bwImBF = im2bw(stdImNrm, stdThresh);
    
    %% Fluorescent channel segmentation
    imAlpha = 1.2;
    
    bwIm = false(size(bfIm));
    for c=1:length(otherChan)
        nzIm = chanIm{otherChan(c)};
        nzIm = nzIm(chanIm{otherChan(c)} > 0);

        fluorThresh = graythresh(nzIm);
        bwIm = bwIm | im2bw(chanIm{otherChan(c)}, imAlpha*fluorThresh);
    end
    
    % Discard small segmentations
    stats = regionprops(bwIm, 'Area', 'PixelIdxList');
    for i=1:length(stats)
        if ( stats(i).Area < 30 )
            bwIm(stats(i).PixelIdxList) = false;
        end
    end
    
    labelIm = bwlabel(bwIm);
    
    components = {};
    
    %% Assign brightfield segmentation pixels to the nearest (Euclidean) nucleus
    cc = bwconncomp(bwImBF);
    for i=1:cc.NumObjects
        isectLabels = unique(labelIm(cc.PixelIdxList{i}));
        isectLabels = isectLabels(isectLabels > 0);
        
        if ( isempty(isectLabels) )
            continue;
        end
        
        compIm = ismember(labelIm,isectLabels);
        [distIm,idxIm] = bwdist(compIm);
        
        bwComp = false(size(labelIm));
        bwComp(cc.PixelIdxList{i}) = true;
        
        assignIm = labelIm(idxIm);
        assignIm(~bwComp) = 0;
        for j=1:length(isectLabels)
            components = [components; {find(assignIm == isectLabels(j))}];
        end
    end

    %% Build a hull structure for each assigned component.
    for i=1:length(components)
        [r,c] = ind2sub(size(bwImBF), components{i});
        chIdx = Helper.ConvexHull(c,r);
        if ( isempty(chIdx) )
            continue;
        end
        
        nh = struct('indexPixels',{components{i}}, 'points',{[c(chIdx), r(chIdx)]});
        hulls = [hulls nh];
    end
end
