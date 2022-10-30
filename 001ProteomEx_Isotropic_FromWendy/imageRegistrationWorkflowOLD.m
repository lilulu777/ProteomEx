function [processedImagePair] = imageRegistrationWorkflow(imPair)

% Imoving and Istatic are either matrices containing images or full paths
% to images intended to be read
% processedImagePair is a structure with all relevant results from operations in
% this function


% ImgMoving.mask is the same size as ImgMoving.windowed, which is windowed
% according to the limits in window.x and window.y

ImgMoving = imPair.ImgMoving;
ImgStatic = imPair.ImgStatic;


npixX = 200;
npixY = 200;
x0 = 0;
y0 = 0;
window.x = [1 npixX] + x0;
window.y = [1 npixY] + y0;


needToRoughAlign = 0;
% ImgMoving.windowed = ImgMoving.input;
% ImgStatic.windowed = ImgStatic.input;

needToPickPoints = 0;
% ImgMoving.moving_out = imPair.sections{1,1}.moving.ctlPoints;
% ImgStatic.fixed_out = imPair.sections{1,1}.static.ctlPoints;

needToMask = 0;
% ImgMoving.mask = imPair.sections{1,1}.moving.mask;

needToRegister = 1;







%make sure pixel values are properly scaled doubles.

ImgStatic.input = ImgStatic.input/max(ImgStatic.input(:));


ImgMoving.input = ImgMoving.input/max(ImgMoving.input(:));



%% pick two points get rough translation, rotation and scaling
%Imoving is registered to Istatic and trimmed/padded to the same field of
%view as Istatic to produce roughAligned.moving
if needToRoughAlign
    [moving_out0,fixed_out0] = cpselect(ImgMoving.input,ImgStatic.input,'Wait',true);
    tform0 = cp2tform(moving_out0,fixed_out0,'nonreflective similarity');
    sizeStatic = size(ImgStatic.input);
    ImgMoving.roughAligned = imtransform(ImgMoving.input,tform0,'XData',[1 sizeStatic(2)],'YData',[1 sizeStatic(1)],'FillValues',0);
else
    ImgMoving.roughAligned = ImgMoving.input;
end

% %Once rough alignment is done, window both images
% ImgMoving.windowed = ImgMoving.roughAligned(window.y(1):window.y(2),window.x(1):window.x(2));
% ImgStatic.windowed = ImgStatic.input(window.y(1):window.y(2),window.x(1):window.x(2));


%% masks

if needToMask
    sd = 8; %Use 8 pixel s.d. blurring
    gFilter = fspecial('gaussian', sd*4, sd);
    blurred.Moving = imfilter(ImgMoving.windowed,gFilter,'replicate');
    blurred.Static = imfilter(ImgStatic.windowed,gFilter,'replicate');
    
    %initial guess for threshold value
    thresh0 = 0.3;
    
    cmapblackandred = double(zeros(64,3));
    thresh0Index = 1+floor(63*thresh0);
    cmapblackandred([1:thresh0Index],1) = double(ones(thresh0Index,1));
    greyValues = linspace(0,1,64-thresh0Index);
    cmapblackandred([thresh0Index + 1:64],1) = greyValues;
    cmapblackandred([thresh0Index + 1:64],2) = greyValues;
    cmapblackandred([thresh0Index + 1:64],3) = greyValues;

    h = figure; imagesc(horzcat(blurred.Moving,blurred.Static)); colormap(cmapblackandred);%figure out threshold by hand
    colormapeditor;
    display('Adjust colormap to mask red areas correctly and press any key')
    pause on
    pause
    cmap = colormap;
    ind = find(cmap(:,1)-1,1);%First non-zero elements are all 1. After subtracting 1, the formerly 0 element is the first non-zero element
    close(h)
    thresh1 = (ind - 1)/63;
    
    ImgMoving.mask = blurred.Moving >= thresh1;
end



%% manual pre-registration

if needToPickPoints
    [ImgMoving.moving_out,ImgStatic.fixed_out] = cpselect(ImgMoving.windowed,ImgStatic.windowed,ImgMoving.moving_out,ImgStatic.fixed_out,'Wait',true);
% else
%     ImgMoving.moving_out = [];
%     ImgStatic.fixed_out = [];
end
    

%% non-rigid registration

if needToRegister

    %renormalize images. Need to fix later to ignore masked bits. Probably
    %need two separate masks with separate thresholds.
%     hgramBuckets = 60;
%     hgram = hist(double(ImgStatic.windowed(:)),hgramBuckets);
%     ImgMoving.windowed = histeq(ImgMoving.windowed,hgram);

    Options.Similarity = 'sd';
    %sd: squared distance, gd: gradient distance, gc: gradient correlation
    
%     sd = 2; %Use 8 pixel s.d. blurring
%     gFilter = fspecial('gaussian', sd*4, sd);
%     ImgMoving.smoothed = imfilter(ImgMoving.windowed,gFilter,'replicate');
%     ImgStatic.smoothed = imfilter(ImgStatic.windowed,gFilter,'replicate');

    Options.Spacing = 4*2^3 * [1 1];%initial grid spacing
    Options.Penalty = 1e-3;%default is 1e-3
    Options.Registration = 'NonRigid';
    Options.Verbose = 1;%default is 2 (max value)
%    Options.MaskMoving = [];
    Options.MaskMoving = ImgMoving.mask;
    Options.MaskStatic = ones(size(ImgMoving.mask));%Must provide, even if it is fake
    Options.Points1 = ImgMoving.moving_out;%for Imoving
    Options.Points2 = ImgStatic.fixed_out;%for Istatic
    Options.PStrength = ones(size(Options.Points1,1),1);%might want to check and make sure this makes it a bit stronger than elastic penalty
    Options.MaxRef = 3;
    
    
    [Ireg,O_trans,Spacing,M,B,F] = image_registration(ImgMoving.windowed,ImgStatic.windowed,Options);
    ImgMoving.registered = Ireg;
    ImgMoving.O_trans = O_trans;
    ImgMoving.Spacing = Spacing;
%    ImgMoving.M = M;
    ImgMoving.B = B;
    ImgMoving.F = F;

end

processedImagePair.ImgMoving = ImgMoving;
processedImagePair.ImgStatic = ImgStatic;

transformQuiverOLD(processedImagePair)

