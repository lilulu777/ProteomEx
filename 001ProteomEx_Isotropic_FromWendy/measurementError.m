function [imgPair] = measurementError(imgPair,sampleInterval,autoSave,readSaved)

computeError = 1;
displayError = 1;

F = imgPair.moving.F;
mask = imgPair.moving.mask;


if computeError


    mask = double(mask);
    nanMask = mask;
    nanMask(mask==0) = NaN;
    nanMask = cat(3,nanMask,nanMask);

    sizeImg = size(F);

    F = F.*nanMask;%Ignore transform vectors that map from parts of specimen with no features

    Fpadded = NaN(sizeImg(1)*3-2,sizeImg(2)*3-2,2);
    Xwin = [1:sizeImg(2)]+sizeImg(2)-1;
    Ywin = [1:sizeImg(1)]+sizeImg(1)-1;
    Fpadded(Ywin,Xwin,:) = F;

    deltaField = NaN(sizeImg);
    deltaMagSqField = NaN(sizeImg(1),sizeImg(2));
    deltaMagSqVec = NaN(sizeImg(1)*sizeImg(2),1);

    %measErrorMean = NaN([sizeImg(1) sizeImg(2)]);
    %measErrorSD = NaN([sizeImg(1) sizeImg(2)]);
    %measLength = NaN([sizeImg(1) sizeImg(2)]);

    maxR = ceil((sizeImg(1)^2 + sizeImg(2)^2)^(0.5));

%    r=1;

%    sampleInterval = 20;

    % 20180405 added: saving intermediate workspaces for restarting
    % For initial run, set "starting_i" to be same as "sampleInterval"
    % For later runs, set "starting_i" to be the value stored in saved file

    if readSaved
        load Workspace_During_ErrorAnalysis measResults
        load Workspace_During_ErrorAnalysis i
        starting_i = i + sampleInterval; 
    else
        measResults = zeros(maxR,6);
        %row# = ceiling(measurement length)
        %col1: sum of squared error values
        %col2: num of measurements
        %col3: rms error
        %col4: sum of y errors
        %col5: sum of x errors
        %col6: mean error magnitude
        starting_i = sampleInterval;
    end
    % Starting at the next round, to avoid double counting Round i
    
    
    for i = starting_i:sampleInterval:2*sizeImg(1)-2
        for j = 1:sampleInterval:2*sizeImg(2)-2
            deltaField = deltaField*NaN;
            Fmoving = Fpadded([1:sizeImg(1)]+i-1,[1:sizeImg(2)]+j-1,:);
            deltaField = Fmoving - F;
            deltaMagSqField = (deltaField(:,:,1).^2 + deltaField(:,:,2).^2);
            deltaMagSqVec = deltaMagSqField(:);
            r = ceil(((i-sizeImg(1))^2 + (j-sizeImg(2))^2)^(0.5));
            %disp(r)
            measResults(r,1) = measResults(r,1) + nansum(deltaMagSqVec);
            measResults(r,2) = measResults(r,2) + nansum(0*deltaMagSqVec+1);
            
            measResults(r,4) = measResults(r,4) + nansum(nansum(deltaField(:,:,1)));
            measResults(r,5) = measResults(r,5) + nansum(nansum(deltaField(:,:,2)));
        end
        
        if autoSave
            save Workspace_During_ErrorAnalysis.mat measResults i
            disp(['Round i = ' num2str(i) ' is finished and saved.']); 
        else
            disp(['Round i = ' num2str(i) ' is finished.']);
        end
    end

    measResults(:,3) = (measResults(:,1)./measResults(:,2)).^(0.5);
    measResults(:,6) = measResults(:,4).^2 + measResults(:,5).^2;
    measResults(:,6) = measResults(:,6).^0.5./measResults(:,2);

    imgPair.measResults = measResults;
end



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

if displayError
    measResults = imgPair.measResults;
    xmin = 1;%in pixels
    xmax = size(measResults,1);
    xvals = [xmin:xmax];
    notNaNindices = xvals(~isnan(measResults([xmin:xmax],3)));
    xvals = notNaNindices * imgPair.info.pixelWidth;
    
    figure1 = figure('Color',[1 1 1]);

%    [AX H1 H2] = plotyy(xvals,measResults(notNaNindices,3)*imgPair.info.pixelWidthExpanded,xvals,measResults(notNaNindices,3)./notNaNindices'*100)

%set(figure1, 'Units','inches', 'Position',[0 0 6 6]);
% figure size printed on paper
%set(figure1, 'PaperUnits','inches')
%set(figure1, 'PaperSize',[6 6])
%set(figure1, 'PaperPosition',[0  0 6 6])
%set(figure1, 'PaperOrientation','portrait')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The current version is after fixing the ExpFactor Bug  %%%
%%% (a divider is removed from both x and y display values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2018/09/20 Revised again
%%% Assumption is the measResult(3,:) has values in units of pixels
%%% of the input images (both converted to pre-ExM-pixels).
%%% Hence, to convert to pre-ExM um, 
%%% y-values need to multiply by the
%%% factor pixelWidthExpanded / ExpFactor
%%% x-values are already in post-expansion um, so divide by
%%% the factor ExpFactor
%%% THIS RETURNS TO THE EXACT SAME CODE AS IN FEI'S ORIGINAL SCRIPT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2019/11/07 Revised again 2
%%% Assumption is the measResult(3,:) has values in units of pixels
%%% of the input image (depending on pixel size of the registration
%%% target image; 

AX = axes('Parent',figure1,'YTick',[0 0.1 0.2 0.3 0.4],...
    'XTick',[0 20 40 60 80],...
    'FontSize',20);
H1 = plot(AX, ...
    xvals, ...
    measResults(notNaNindices,3)*imgPair.info.pixelWidth, 'LineWidth', 2)
    
% Create ylabel
ylabel('Measurement RMS Error (um)','FontSize',20);

% Create xlabel
xlabel('Measurement Length (um)','FontSize',20);



end

%%

% sqDisp = F(:,:,1).^2 + F(:,:,2).^2;
% meanSqDisp = mean(sqDisp(:));
% RMSError = meanSqDisp^0.5;
% expectedLongRangeError = 2^0.5*RMSError*imgPair.info.pixelWidthExpanded





