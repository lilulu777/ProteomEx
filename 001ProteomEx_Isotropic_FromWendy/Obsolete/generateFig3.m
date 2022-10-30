%%
figure
%set(gcf,'Position',[756   533   614   247]);

subplot(1,2,1)
%set(gca,'Position',[0.1300    0.1100    0.3347    0.8150]);
imshowpair(I_SIM,I_expanded)
%title('Unregistered, with deformation field')
axis on
hold on
quiver(x,y,Tformdec(:,:,2),Tformdec(:,:,1),0,'Color','white');
set(gca,'YDir','reverse');
set(gca,'FontSize',14)
%subplot(1,2,2)
%set(gca,'Position',[0.5703    0.1100    0.3347    0.8150]);
%imshowpair(I_SIM,I_exp_registered)
%title('Registered Images')
%xis on


xmin = 5;
xmax = 650;
xvals = [xmin:xmax];
notNaNindices = xvals(~isnan(measResults([xmin:xmax],3)));
xvals = notNaNindices * 0.105;
subplot(1,2,2)
[AX H1 H2] = plotyy(xvals,measResults(notNaNindices,3)./notNaNindices'*100,xvals,measResults(notNaNindices,3)*.105)
xlabel('measurement length (um)','FontSize',14)

set(get(AX(1),'Ylabel'),'String','rms error of measurement (%)','FontSize',14) 
set(get(AX(2),'Ylabel'),'String','rms error of measurement (um)','FontSize',14) 

set(get(AX(1)),'FontSize',14)
set(get(AX(2)),'FontSize',14)

set(gcf,'Color','white')
