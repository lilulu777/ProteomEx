
%% local renormalization approach (no mask)
%blur out all features, leaving just a local scaling factor
sd = 200;
gFilter = fspecial('gaussian', 4*sd, sd);
blurred = imfilter(imgToFilter,gFilter,'replicate');
imgRenormed = imgToFilter*max(blurred(:))./blurred;
figure;subplot(1,3,1);imagesc(imgToFilter);subplot(1,3,2);imagesc(blurred);subplot(1,3,3);imagesc(imgRenormed);
% doesn't work, try again.

%% try again
WindowHalfWidth = 150;

grainSize = 10;

imgRenormed = 0*imgToFilter;
maxPix = max(imgToFilter(:));

sizeImg = size(imgToFilter);

localMax = 0;

for i = 1:sizeImg(1)
    if mod(i,grainSize) == 0;i;end
    for j = 1:sizeImg(2)
        ymin = i - WindowHalfWidth;
        ymax = i + WindowHalfWidth;
        xmin = j - WindowHalfWidth;
        xmax = j + WindowHalfWidth;

        if j ==1 || i == 1 || (mod(i,grainSize) == 0 && mod(j,grainSize) == 0)
            localMax = max(max(imgToFilter([max(1,ymin):min(ymax,sizeImg(1))],[max(1,xmin):min(xmax,sizeImg(2))])));
        %else keep localMax from previous window
        end
        
        imgRenormed(i,j) = imgToFilter(i,j)*maxPix/localMax;
    end
end
figure;subplot(1,2,1);imagesc(imgToFilter);subplot(1,2,2);imagesc(imgRenormed);
%also doesn't work because the variation in average intensity is on too
%fine a scale relative to the maximum distance between macro-scale features



