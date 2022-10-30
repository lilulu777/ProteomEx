
imgToFilter = double(IafterWithZeros);

gValues = [0 0.5 1 2 4 8 16];

N = 200;

histTraces = zeros(N,length(gValues),2);
figure
%subplot(1,2,1);
hold on
%set(gca,'Ylim',[0 1e5])
%set(gca,'Xlim',[0 100])

%subplot(1,2,2);
hold on

i = 1;

load 'cmapblackandred.mat'

[histTraces(:,i,2) histTraces(:,i,1)] = hist(imgToFilter(:),N);
%subplot(1,2,1)
%plot(histTraces(:,i,1),histTraces(:,i,2),'DisplayName',['Gaussian blur:' num2str(gValues(i))]);
%subplot(1,2,2)
imagesc(imgToFilter,'DisplayName',['Gaussian blur:' num2str(gValues(i))]);
colormap(cmapblackandred);

for i = 2:length(gValues)
    gFilter = fspecial('gaussian', gValues(i)*4, gValues(i));
    blurred = imfilter(imgToFilter,gFilter,'replicate');
    [histTraces(:,i,2) histTraces(:,i,1)] = hist(blurred(:),N);
%    subplot(1,2,1)
%    plot(histTraces(:,i,1), histTraces(:,i,2),'DisplayName',['Gaussian blur:' num2str(gValues(i))]);
%    subplot(1,2,2)
    imagesc(blurred,'DisplayName',['Gaussian blur:' num2str(gValues(i))]);
    colormap(cmapblackandred);
end
    

propertyeditor('on')


figure;plot(histTraces(:,:,1),histTraces(:,:,2))


%% integrate histograms
histTraces(1,:,3) = histTraces(1,:,2);
for i = 2:N
    histTraces(i,:,3) = histTraces(i-1,:,3) + histTraces(i,:,2);
end

figure;semilogx(histTraces(:,:,1),histTraces(:,:,3))


%% Actual algorithm begins here

sd = 4; %Use 8 pixel s.d. blurring
gFilter = fspecial('gaussian', sd*4, sd);
blurred = imfilter(imgToFilter,gFilter,'replicate');

figure; imagesc(blurred); colormap(cmapblackandred);


%% make mask
threshold = 4100;

sizeImg = size(blurred);
mask = zeros(sizeImg);

for i = 1:sizeImg(1)
    for j = 1:sizeImg(2)
        if blurred(i,j) < 4100; mask(i,j) = 0;
        else mask(i,j) = 1;
        end
    end
end

masked = imgToFilter .* mask;
figure;subplot(1,3,1);imagesc(imgToFilter);subplot(1,3,2);imagesc(mask);subplot(1,3,3);imagesc(masked);

