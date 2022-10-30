function [] = measurementErrorPruned(F,mask)

mask = double(mask);
nanMask(mask==0) = NaN;
nanMask = cat(3,mask,mask);

sizeImg = size(F);

F = F.*nanMask;%Ignore transform vectors that map from parts of specimen with no features

n = 0;
Freshaped = NaN(size(F,1)*size(F,2),4);
for i = 1:size(F,1);
    for j = 1:size(F,2);
        if mask(i,j)==1
            n = n + 1;
            Freshaped(n,1) = i;
            Freshaped(n,2) = j;
            Freshaped(n,3) = F(i,j,1);
            Freshaped(n,4) = F(i,j,2);
        end
        
    end
end
Freshaped = Freshaped(1:n,:);
clear n
%Freshaped has all entries of F that are not masked out.
%col1: y pos, col2: x pos, col3: delta y, vol4: delta x

Fpadded = NaN(size(Freshaped,1)*3-2,4);
Fpadded((1:size(Freshaped,1))+size(Freshaped,1)-1,:) = Freshaped;

deltaField = NaN(size(Freshaped));
deltaMagSqField = NaN(size(Freshaped,1),2);%col1 = r, col2 = (delta r)^2
%deltaMagSqVec = NaN(sizeImg(1)*sizeImg(2),1);

%measErrorMean = NaN([sizeImg(1) sizeImg(2)]);
%measErrorSD = NaN([sizeImg(1) sizeImg(2)]);
%measLength = NaN([sizeImg(1) sizeImg(2)]);

maxR = ceil((sizeImg(1)^2 + sizeImg(2)^2)^(0.5));

measResults = zeros(maxR,3);
%row# = ceiling(measurement length)
%col1: sum of squared error values
%col2: num of measurements
%col3: rms error

r=1;

sampleInterval = 10;
size(Freshaped)
for i = sampleInterval:sampleInterval:2*size(Freshaped,1)-2
    if mod(i,1000) == 0;
        i
    end
    
    deltaField = deltaField*NaN;
    
    Fmoving = Fpadded([1:size(Freshaped,1)]+i-1,:);
    
    deltaField = Fmoving - Freshaped;
    
    deltaMagSqField(:,1) = ceil(((deltaField(:,1)).^2 + (deltaField(:,2)).^2).^(0.5));
    deltaMagSqField(:,2) = (deltaField(:,3).^2 + deltaField(:,4).^2);

    for j = 1:size(deltaField,1)
        r = deltaMagSqField(j,1);
        if r >0
            measResults(r,1) = measResults(r,1) + deltaMagSqField(j,2);
            measResults(r,2) = measResults(r,2) + 1;
        end
    end
end

measResults(:,3) = (measResults(:,1)./measResults(:,2)).^(0.5);

%% put data into buckets according to measurement length with single pixel resolution.


% figure;
% subplot(1,2,1);plot(measurementsAvg(:,1),'-k','DisplayName','mean error');hold on;plot(measurementsAvg(:,2),'DisplayName','S.D.');
% xlabel('measurement length')
% legend('show')
% subplot(1,2,2);plot(radAvgValues(:,1),radAvgValues(:,2)./radAvgValues(:,1),'-k','DisplayName','mean fractional error');hold on;plot(radAvgValues(:,1),radAvgValues(:,3)./radAvgValues(:,1),'DisplayName','S.D.');
% xlabel('measurement length')
% legend('show')

%figure;
%subplot(1,2,1);plot(radAvgValues(:,1),radAvgValues(:,2),'-k','DisplayName','mean error');hold on;plot(radAvgValues(:,1),radAvgValues(:,3),'DisplayName','S.D.');
%xlabel('measurement length')
%legend('show')
%subplot(1,2,2);plot(radAvgValues(:,1),radAvgValues(:,2)./radAvgValues(:,1),'-k','DisplayName','mean fractional error');hold on;plot(radAvgValues(:,1),radAvgValues(:,3)./radAvgValues(:,1),'DisplayName','S.D.');
%xlabel('measurement length')
%legend('show')

%% display
xmin = 4;%in pixels
xmax = size(measResults,1);
xvals = [xmin:xmax];
notNaNindices = xvals(~isnan(measResults([xmin:xmax],3)));
xvals = notNaNindices * 0.105;
figure
[AX H1 H2] = plotyy(xvals,measResults(notNaNindices,3)*.105,xvals,measResults(notNaNindices,3)./notNaNindices'*100)
xlabel('measurement length (um)','FontSize',14)

set(get(AX(1),'Ylabel'),'String','rms error of measurement (um)','FontSize',14) 
set(get(AX(2),'Ylabel'),'String','rms error of measurement (%)','FontSize',14) 


%%

sqDisp = F(:,:,1).^2 + F(:,:,2).^2;
meanSqDisp = mean(sqDisp(:));
RMSError = meanSqDisp^0.5;
expectedLongRangeError = 2^0.5*RMSError*.105






