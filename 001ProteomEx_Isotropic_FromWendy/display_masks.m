function display_masks(imgPair)

figure;clf;
subplot(1,2,1); imagesc(imgPair.static.mask);
subplot(1,2,2); imagesc(imgPair.moving.mask);