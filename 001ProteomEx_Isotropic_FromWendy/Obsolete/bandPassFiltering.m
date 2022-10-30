
%%
if 0
b = firpm(100,[0 0.1 0.11 1],[0 1 1 1]);
h = ftrans2(b);
[H,w] = freqz(b,1,64,'whole');
figure;
set(gcf,'Position',[    200         500        1600         440]);
colormap(jet(64))
subplot(1,3,1); imagesc(imgPair.moving.input)
subplot(1,3,2);plot(w/pi-1,fftshift(abs(H)))
%subplot(1,3,2); freqz2(h,[32 32])

subplot(1,3,3);imagesc(imfilter(imgPair.moving.input,h))
end

%%
if 1
widthN = 200;
neighborhood = zeros(2*widthN + 1,2);
img = imgPair.moving.input;
imgPadded = zeros(size(img) + 2*widthN);
imgPadded((1:size(img,1))+widthN,(1:size(img,2))+widthN) = img;
maxInNeighborhood = zeros(size(img));

patchSize = 30;
patch = ones(2*patchSize,2*patchSize);

for i = 1+patchSize:2*patchSize:size(img,1)
    for j = 1+patchSize:2*patchSize:size(img,2)
        neighborhood = imgPadded((1:widthN)+i-1,(1:widthN)+j-1);
        patch = 0*patch + max(neighborhood(:));
        maxInNeighborhood((1:size(patch,1))+i-patchSize,(1:size(patch,2))+j-patchSize) = patch;
    end
end
imgProc = img./maxInNeighborhood;
figure;imagesc(imgProc)
end


