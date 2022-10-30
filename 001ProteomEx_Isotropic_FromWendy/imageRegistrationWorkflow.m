function [imgPair] = imageRegistrationWorkflow(imgPair)

% imgPair must contain at least the fields moving.input and static.input,
% which are the matrices containing input images

% output imgPair is a structure with all relevant results from operations in
% this function


% moving.mask is the same size as moving.windowed, which is windowed
% according to the limits in window.x and window.y

moving = imgPair.moving;
static = imgPair.static;
if isfield(imgPair,'sections'); sections = imgPair.sections; sectionInfo = imgPair.sectionInfo; end

npixX = 2000;
npixY = 2000;

% x0 = 0;
% y0 = 0;

x0 = 1;
y0 = 1;
window.x = [1 npixX] + x0;
window.y = [1 npixY] + y0;

needToRoughAlign = 1;
needToMask = 0;
needToDoSections = 0;
needToPickPoints = 0;
needToRegister = 0;

%make sure pixel values are properly scaled doubles.
static.input = double(static.input);
static.input = static.input/max(static.input(:));

moving.input = double(moving.input);
moving.input = moving.input/max(moving.input(:));


%% pick two points get rough translation, rotation and scaling
%Imoving is registered to Istatic and trimmed/padded to the same field of
%view as Istatic to produce roughAligned.moving

if needToRoughAlign
    [moving_out0,fixed_out0] = cpselect(moving.input,static.input,'Wait',true);
    tform0 = cp2tform(moving_out0,fixed_out0,'nonreflective similarity');
    sizeStatic = size(static.input);
    moving.roughAligned = imtransform(moving.input,tform0,'XData',[1 sizeStatic(2)],'YData',[1 sizeStatic(1)],'FillValues',0);


%Once rough alignment is done, window both images
moving.windowed = moving.roughAligned(window.y(1):window.y(2),window.x(1):window.x(2));
static.windowed = static.input(window.y(1):window.y(2),window.x(1):window.x(2));
end

%% masks

if needToMask
    sd = 8; %Use 8 pixel s.d. blurring
    gFilter = fspecial('gaussian', sd*4, sd);
    blurred.Moving = imfilter(moving.windowed,gFilter,'replicate');
    blurred.Static = imfilter(static.windowed,gFilter,'replicate');
    
    %initial guess for threshold value
    thresh0 = 0.111;
    
    cmapblackandred = double(zeros(64,3));
    thresh0Index = 1+floor(63*thresh0);
    cmapblackandred([1:thresh0Index],1) = double(ones(thresh0Index,1));
    greyValues = linspace(0,1,64-thresh0Index);
    cmapblackandred([thresh0Index + 1:64],1) = greyValues;
    cmapblackandred([thresh0Index + 1:64],2) = greyValues;
    cmapblackandred([thresh0Index + 1:64],3) = greyValues;

    h = figure; imagesc(blurred.Moving); colormap(cmapblackandred);%figure out threshold by hand
    colormapeditor;
    display('Adjust colormap to mask red areas correctly and press any key')
    pause on
    pause
    cmap = colormap;
    ind = find(cmap(:,1)-1,1);%First non-zero elements are all 1. After subtracting 1, the formerly 0 element is the first non-zero element
    close(h)
    thresh1 = (ind - 1)/63;
    
    moving.mask = blurred.Moving >= thresh1;
end




%% define local sections on which to run non-rigid registration

if needToDoSections
    sectionInfo.width = 100;
    sectionInfo.overlap = 20;
    if sectionInfo.overlap > sectionInfo.width
        error('Make sure section width is greater than section overlap')
    end
    if mod(size(moving.windowed,1),sectionInfo.width) ~= 0 || mod(size(moving.windowed,2),sectionInfo.width)
        error('Make sure section width is an integer multiple of windowed image size')
    end

    
    moving.padded = zeros(size(moving.windowed) + 2*sectionInfo.overlap);
    moving.padded([1:size(moving.windowed,1)] + sectionInfo.overlap,[1:size(moving.windowed,2)] + sectionInfo.overlap) = moving.windowed;
    moving.maskPadded = zeros(size(moving.windowed) + 2*sectionInfo.overlap);
    moving.maskPadded([1:size(moving.windowed,1)] + sectionInfo.overlap,[1:size(moving.windowed,2)] + sectionInfo.overlap) = moving.mask;
    static.padded = zeros(size(static.windowed) + 2*sectionInfo.overlap);
    static.padded([1:size(static.windowed,1)] + sectionInfo.overlap,[1:size(static.windowed,2)] + sectionInfo.overlap) = static.windowed;




    %Specify subfigure windows as arrays of x and y values to serve as [1 1]
    %pixel for each section. Pixels positions are referenced to original
    %windowed images, not zero-padded images, will later shift into
    %zero-padded coordinate system
    sectionTemp.origins(:,2) = [1:sectionInfo.width:size(moving.windowed,2)];
    sectionTemp.origins(:,1) = [1:sectionInfo.width:size(moving.windowed,1)];
    %Include overlaps. Images will be padded out on all sides with zeros to
    %include overlap area in all sectionInfo

    sections = cell(ceil(size(moving.windowed)/sectionInfo.width));
    for i = 1:size(sections,1)
        for j = 1:size(sections,2)
            sections{i,j}.static.windowed = static.padded([1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(i-1),[1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(j-1));
            hgramBuckets = 60;
            hgram = hist(sections{i,j}.static.windowed(:),hgramBuckets);

            %each section of moving image is hist-equalized to corresponding
            %static section
            sections{i,j}.moving.windowed = moving.padded([1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(i-1),[1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(j-1));
            sections{i,j}.moving.windowed = histeq(sections{i,j}.moving.windowed,hgram);

            sections{i,j}.moving.mask = moving.maskPadded([1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(i-1),[1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(j-1));

            sections{i,j}.moving.ctlPoints = [];
            sections{i,j}.static.ctlPoints = [];
        end
    end

end

%% manual pre-registration

if needToPickPoints
    for i = 1:size(sections,1)
        for j = 1:size(sections,2)
            if isfield(sections{i,j}.moving,'ctlPoints') && ~isempty(sections{i,j}.moving.ctlPoints)
                [sections{i,j}.moving.ctlPoints,sections{i,j}.static.ctlPoints] = cpselect(sections{i,j}.moving.windowed,sections{i,j}.static.windowed,sections{i,j}.moving.ctlPoints,sections{i,j}.static.ctlPoints,'Wait',true);
            else
                [sections{i,j}.moving.ctlPoints,sections{i,j}.static.ctlPoints] = cpselect(sections{i,j}.moving.windowed,sections{i,j}.static.windowed,'Wait',true);
            end
        end
    end
end
    

%% non-rigid registration

if needToRegister

    pool = parpool;
    
    registered = zeros(size(moving.windowed));
    B = zeros([size(moving.windowed) 2]);
    F = zeros([size(moving.windowed) 2]);

    indices = zeros(i*j,2);

    n = 1;
    for i = 1:size(sections,1)
        for j = 1:size(sections,2)
            indices(n,1) = i;
            indices(n,2) = j;
            n = n + 1;
        end
    end
    clear n

    
    parfor n = 1:size(indices,1)
        i = indices(n,1);
        j = indices(n,2);

        Options.Similarity = 'sd';
        %sd: squared distance, gd: gradient distance, gc: gradient correlation

        finalGridSpacing = 8;
        Options.MaxRef = 3;%max number of grid refinements (not incl initial)
        Options.Spacing = finalGridSpacing*2^Options.MaxRef * [1 1];%initial grid spacing
        Options.Penalty = 1e-3;%default is 1e-3
        Options.Registration = 'NonRigid';
        Options.Verbose = 1;%default is 2 (max value)
    %    Options.MaskMoving = [];
        Options.MaskMoving = sections{i,j}.moving.mask;
        Options.MaskStatic = ones(size(sections{i,j}.static.windowed));%Must provide somethign here



        Options.Points1 = fliplr(sections{i,j}.moving.ctlPoints);%for Imoving
        Options.Points2 = fliplr(sections{i,j}.static.ctlPoints);%for Istatic
        Options.PStrength = ones(size(Options.Points1,1),1);%might want to check and make sure this makes it a bit stronger than elastic penalty


        [Ireg,O_trans,Spacing,Msect,Bsect,Fsect] = image_registration(sections{i,j}.moving.windowed,sections{i,j}.static.windowed,Options);
        sections{i,j}.moving.registered = Ireg;
        sections{i,j}.moving.O_trans = O_trans;
        sections{i,j}.moving.Spacing = Spacing;
    %    sections{i,j}.moving.M = Msect;
        sections{i,j}.moving.B = Bsect;
        sections{i,j}.moving.F = Fsect;

        inputFrameYvals = [1:sectionInfo.width]+sectionInfo.width*(i-1);
        inputFrameXvals = [1:sectionInfo.width]+sectionInfo.width*(j-1);
        registered(inputFrameYvals,inputFrameXvals) = sections{i,j}.moving.registered([1:sectionInfo.width]+sectionInfo.overlap,[1:sectionInfo.width]+sectionInfo.overlap);
        B(inputFrameYvals,inputFrameXvals,:)          = sections{i,j}.moving.B([1:sectionInfo.width]+sectionInfo.overlap,[1:sectionInfo.width]+sectionInfo.overlap,:);
        F(inputFrameYvals,inputFrameXvals,:)          = sections{i,j}.moving.F([1:sectionInfo.width]+sectionInfo.overlap,[1:sectionInfo.width]+sectionInfo.overlap,:);
        %subplot(2,2,(j-1)+2*(i-1)+1);imagesc(F(:,:,1).^2+F(:,:,2).^2);

    end

    moving.registered = registered;
    moving.B = B;
    moving.F = F;

end

imgPair.moving = moving;
imgPair.static = static;

if exist('sections','var');
    imgPair.sections = sections;
    imgPair.sectionInfo = sectionInfo;
end

transformQuiver(imgPair)

