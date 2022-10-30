function [] = transformQuiver(processedImagePair)

%% Figure out axis convention
% x = [10:10:50];
% y = [10:10:30];
% [x,y] = meshgrid(x,y);
% figure;
% quiver(x,y,2*ones(size(x)),ones(size(x)),0)
% set(gca,'YDir','reverse');


%% real stuff
%decimate transform field
spacing = 10;
mytform = processedImagePair.ImgMoving.F;
sizeTform = size(mytform);
Tformdec = zeros(floor(sizeTform(1)/spacing),floor(sizeTform(2)/spacing));
for i = 1:sizeTform(1)/spacing
    for j = 1:sizeTform(2)/spacing
        Tformdec(i,j,1) = mytform(i*spacing,j*spacing,1);
        Tformdec(i,j,2) = mytform(i*spacing,j*spacing,2);
    end
end
    
y = [spacing:spacing:sizeTform(1)];
x = [spacing:spacing:sizeTform(2)];
[x,y] = meshgrid(x,y);

figure

set(gcf,'Position',[756   533   614   247]);

subplot(1,2,1)
set(gca,'Position',[0.1300    0.1100    0.3347    0.8150]);
imshowpair(processedImagePair.ImgStatic.windowed,processedImagePair.ImgMoving.windowed)
axis on
hold on
quiver(x,y,Tformdec(:,:,2),Tformdec(:,:,1),0,'Color','white');
set(gca,'YDir','reverse');
subplot(1,2,2)
set(gca,'Position',[0.5703    0.1100    0.3347    0.8150]);

imshowpair(processedImagePair.ImgStatic.windowed,processedImagePair.ImgMoving.registered)
axis on

