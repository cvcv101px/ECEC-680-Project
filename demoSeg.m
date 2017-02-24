%% Read unzipped seq 1 data into a 4D image array
inPath = 'APAPData_DKFZ_AVlasov/WellC04_Seq0001';
outpath = 'testData';

im = [];
for c = 1:4
    imPath = fullfile(inPath,['WellC04_Seq0001_XY1C' num2str(c,'%01d') '.tif']);
    imInfo = imfinfo(imPath);
    
    if ( c==1 )
        im = zeros(imInfo(1).Height, imInfo(1).Width, length(imInfo), 4);
    end
    
    for t=1:length(imInfo)
        im(:,:,t,c) = double(imread(imPath, t)) / (2^14-1);
% Uncomment to write out images in LEVER file name format.
%         imwrite(im(:,:,t,c),fullfile(outpath,sprintf('test_c%02d_t%04d_z%04d.tif', c, t, 1)));
    end
end

%% 
figure;
for t=1:size(im,3)
    imagesc(im(:,:,t,4));colormap(gray);hold on;
    
    %% Brightfield texture segmentation
    se = strel('disk',3);
    stdIm = entropyfilt(mat2gray(im(:,:,t,4)),se.getnhood);
    
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
    bwIm = false(size(im,1),size(im,2));
    for c=1:3
        nzIm = im(:,:,t,c);
        nzIm = nzIm(im(:,:,t,c) > 0);

        fluorThresh = graythresh(nzIm);
        bwIm = bwIm | im2bw(im(:,:,t,c), imAlpha*fluorThresh);
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
    
    %%
    [B,labelIm] = bwboundaries(bwImBF);
    for i=1:length(B)
        plot(B{i}(:,2),B{i}(:,1), '-r');
    end
    
    cmap = lines(8);
    for i=1:length(components)
        cmapIdx = randi(8,1);
        [r,c] = ind2sub(size(bwImBF), components{i});
        plot(c,r, '.', 'Color',cmap(cmapIdx,:))
    end
    
    hold off;
end
