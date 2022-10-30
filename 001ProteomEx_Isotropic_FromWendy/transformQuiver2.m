function [] = transformQuiver2(imgPair)

%% Figure out axis convention
% x = [10:10:50];
% y = [10:10:30];
% [x,y] = meshgrid(x,y);
% figure;
% quiver(x,y,2*ones(size(x)),ones(size(x)),0)
% set(gca,'YDir','reverse');


%% real stuff
%decimate transform field
spacing = 20;
mytform = imgPair.moving.F(256:512, 256:512,:);
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

meanStatic = mean(imgPair.static.windowed(:));
meanMoving = mean(imgPair.moving.windowed(:));
if meanStatic > meanMoving
    imStatic = imgPair.static.windowed;
    imMoving = min(meanStatic/meanMoving*imgPair.moving.windowed,1);
else
    imStatic = min(meanMoving/meanStatic*imgPair.static.windowed,1);
    imMoving = imgPair.moving.windowed;
end

subplot(1,2,1)
set(gca,'Position',[0.1300    0.1100    0.3347    0.8150]);
imshowpair(imStatic(256:512, 256:512),imMoving(256:512, 256:512))
axis off
hold on
scale = 2;
hq1 = quiver(x,y,scale*Tformdec(:,:,2),scale*Tformdec(:,:,1),0,'Color','white', 'LineWidth',1);


set(gca,'YDir','reverse');
subplot(1,2,2)
set(gca,'Position',[0.5703    0.1100    0.3347    0.8150]);

meanStatic = mean(imgPair.static.windowed(:));
meanMoving = mean(imgPair.moving.registered(:));
if meanStatic > meanMoving
    imStatic = imgPair.static.windowed;
    imMoving = min(meanStatic/meanMoving*imgPair.moving.registered,1);
else
    imStatic = min(meanMoving/meanStatic*imgPair.static.windowed,1);
    imMoving = imgPair.moving.registered;
end

imshowpair(imStatic(256:512, 256:512),imMoving(256:512, 256:512))
axis off

%get the pair of images onto two separate layers
% zerosMat = zeros(size(imgPair.static.windowed));
% 
% staticRGB = cat(3,imgPair.static.windowed,zerosMat,zerosMat);
% hs = imshow(staticRGB)
% set(hs,'AlphaData',0.5);
% hold on
% 
% registeredRGB = cat(3,zerosMat,imgPair.moving.registered,zerosMat);
% hm = imshow(registeredRGB)
% set(hm,'AlphaData',0.5);
% 
