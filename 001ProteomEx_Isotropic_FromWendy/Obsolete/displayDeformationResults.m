
spacing = 30;
mytform = deformationField;
sizeTform = size(mytform);

%decimate deformation field
Tformdec = zeros(sizeTform(1)/spacing,sizeTform(2)/spacing);
for i = 1:sizeTform(1)/spacing
    for j = 1:sizeTform(2)/spacing
        Tformdec(i,j,1) = mytform(i*spacing,j*spacing,1);
        Tformdec(i,j,2) = mytform(i*spacing,j*spacing,2);
    end
end
    
y = [spacing:spacing:sizeTform(1)];
x = [spacing:spacing:sizeTform(2)];
[x,y] = meshgrid(x,y);


%%
figure
%set(gcf,'Position',[756   533   614   247]);

subplot(1,2,1)
%set(gca,'Position',[0.1300    0.1100    0.3347    0.8150]);
imshowpair(I_SIM,I_expanded)
title('Unregistered, with deformation field')
axis on
hold on
quiver(x,y,Tformdec(:,:,2),Tformdec(:,:,1),0,'Color','white');
set(gca,'YDir','reverse');
%subplot(1,2,2)
%set(gca,'Position',[0.5703    0.1100    0.3347    0.8150]);
%imshowpair(I_SIM,I_exp_registered)
%title('Registered Images')
%xis on
