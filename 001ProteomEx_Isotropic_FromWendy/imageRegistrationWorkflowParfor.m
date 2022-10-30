function [imgPair] = imageRegistrationWorkflowParfor(imgPair,imgSize, doNotSection)

% imgPair must contain at least the fields moving.input and static.input,
% which are the matrices containing input images

% output imgPair is a structure with all relevant results from operations in
% this function


% moving.mask is the same size as moving.windowed, which is windowed
% according to the limits in window.x and window.y

moving = imgPair.moving;
static = imgPair.static;
if isfield(imgPair,'sections'); sections = imgPair.sections; sectionInfo = imgPair.sectionInfo; end

npixX = imgSize(1);
npixY = imgSize(2);

% x0 = 0;
% y0 = 0;

x0 = 0;
y0 = 0;
window.x = [1 npixX] + x0;
window.y = [1 npixY] + y0;

needToRoughAlign    = 1;
needToMask          = 1;
needToDoSections    = 1;  
singleSection = [0 1];%only applies to PickPoints and Register. elem 1 set to 0 to do all sections.
needToPickPoints    = 1;
needToRegister      = 1;
needToReconstruct   = 1;

%make sure pixel values are properly scaled doubles.
static.input = double(static.input);
static.input = static.input/max(static.input(:));

moving.input = double(moving.input);
moving.input = moving.input/max(moving.input(:));


%% pick two points get rough translation, rotation and scaling
%Imoving is registered to Istatic and trimmed/padded to the same field of
%view as Istatic to produce roug    hAligned.moving

if needToRoughAlign
    [moving_out0,fixed_out0] = cpselect(moving.input,static.input,'Wait',true);
    tform0 = cp2tform(moving_out0,fixed_out0,'nonreflective similarity');
    sizeStatic = size(static.input);
    moving.roughAligned = imtransform(moving.input,tform0,'XData',[1 sizeStatic(2)],'YData',[1 sizeStatic(1)],'FillValues',0);


end

%% masks

if needToMask
    
    %Window both images
    %moving.windowed = moving.roughAligned(window.x(1):window.x(2),window.y(1):window.y(2));
    %static.windowed = static.input(window.x(1):window.x(2),window.y(1):window.y(2));
    moving.windowed = moving.roughAligned;
    static.windowed = static.input;
    
    %make sure pixel values are properly scaled doubles.
    static.windowed = static.windowed/max(static.windowed(:));
    moving.windowed = moving.windowed/max(moving.windowed(:));

    sd = 4;
    gFilter = fspecial('gaussian', sd*4, sd);
    blurred.Moving = imfilter(moving.windowed,gFilter,'replicate');
    blurred.Static = imfilter(static.windowed,gFilter,'replicate');
    blurred.Moving = 1-(blurred.Moving - min(blurred.Moving(:)))/(max(blurred.Moving(:))-min(blurred.Moving(:)));
    blurred.Static = 1-(blurred.Static - min(blurred.Static(:)))/(max(blurred.Static(:))-min(blurred.Static(:)));
    %scale both images to [0 1]
    
    %initial guess for threshold value
    thresh0 = 3/64;
    
    h1 = figure; %subplot(1,2,1); 
    imagesc(blurred.Moving); colormap(cmapblackandred(thresh0,0,1));%figure out threshold by hand
    set(h1,'Position',[100 600 550 450]);
    colormapeditor;
    display('Adjust colormap to mask red areas correctly and press any key')
    pause on
    pause
    cmap = colormap;
    ind = find(cmap(:,1)-1,1);%Find first 0 element. First non-zero elements are all 1. After subtracting 1, the formerly 0 element is the first non-zero element
    %threshold should be the lower bound of bucket# ind
    
    thresh1 = (ind - 1)/64;
    
    moving.mask = blurred.Moving >= thresh1;
    nPixMasked = sum(~moving.mask(:));
        
    sorted = sort(blurred.Static(:));
    thresh2 = sorted(nPixMasked);
    thresh2 = floor(thresh2*64)/64;
    
    h2 = figure;
    imagesc(blurred.Static); colormap(cmapblackandred(thresh2,0,1));%figure out threshold by hand
    set(h2,'Position',[700 600 550 450]);
    colormapeditor;
    display('Adjust colormap to mask red areas correctly and press any key')
    pause on
    pause
    cmap2 = colormap;
    ind2 = find(cmap2(:,1)-1,1);%First non-zero elements are all 1. After subtracting 1, the formerly 0 element is the first non-zero element
    thresh2 = (ind2 - 1)/64;
        
    static.mask = blurred.Static >= thresh2;
    pause

    figure(h1); imagesc((moving.windowed*(63/64)+(1/64)).*moving.mask);
    colormap(cmapblackandred(1/64,0,1))
    figure(h2); imagesc((static.windowed*(63/64)+(1/64)).*static.mask);
    colormap(cmapblackandred(1/64,0,1))
    %transform images out of 0-0.2 range so mask is very clear

    pause
    %result = input('Is mask satisfactory? (y/n)','s');

    
    close(h1)
    close(h2)
    
end


%% define local sections on which to run non-rigid registratio
if needToDoSections
    sectionInfo.width = size(moving.roughAligned,1)/2;
    if doNotSection
        sectionInfo.width = size(moving.roughAligned,1);
        sectionInfo.height = size(moving.roughAligned,2);
    end
    sectionInfo.overlap = 20;
    if sectionInfo.overlap > sectionInfo.width
        error('Make sure section width is greater than section overlap')
    end
    if mod(size(moving.windowed,1),sectionInfo.width) ~= 0 || mod(size(moving.windowed,2),sectionInfo.height)
        error('Make sure section width is an integer multiple of windowed image size')
    end

    
    moving.padded = zeros(size(moving.windowed) + 2*sectionInfo.overlap);
    moving.padded([1:size(moving.windowed,1)] + sectionInfo.overlap,[1:size(moving.windowed,2)] + sectionInfo.overlap) = moving.windowed;
    moving.maskPadded = zeros(size(moving.windowed) + 2*sectionInfo.overlap);
    moving.maskPadded([1:size(moving.windowed,1)] + sectionInfo.overlap,[1:size(moving.windowed,2)] + sectionInfo.overlap) = moving.mask;

    static.padded = zeros(size(static.windowed) + 2*sectionInfo.overlap);
    static.padded([1:size(static.windowed,1)] + sectionInfo.overlap,[1:size(static.windowed,2)] + sectionInfo.overlap) = static.windowed;
    static.maskPadded = zeros(size(static.windowed) + 2*sectionInfo.overlap);
    static.maskPadded([1:size(static.windowed,1)] + sectionInfo.overlap,[1:size(static.windowed,2)] + sectionInfo.overlap) = static.mask;




    %Specify subfigure windows as arrays of x and y values to serve as [1 1]
    %pixel for each section. Pixels positions are referenced to original
    %windowed images, not zero-padded images, will later shift into
    %zero-padded coordinate system
    sectionTemp.origins(:,2) = [1:sectionInfo.width:size(moving.windowed,2)];
    sectionTemp.origins(:,1) = [1:sectionInfo.width:size(moving.windowed,1)];
    %Include overlaps. Images will be padded out on all sides with zeros to
    %include overlap area in all sectionInfo

    sections = cell(size(moving.windowed,1)/sectionInfo.width,size(moving.windowed,2)/sectionInfo.height);
    for i = 1:size(sections,1)
        for j = 1:size(sections,2)
            disp(i);
            disp(j);
            sections{i,j}.moving.mask = moving.maskPadded([1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(i-1),[1:sectionInfo.height+2*sectionInfo.overlap]+sectionInfo.height*(j-1));
            sections{i,j}.static.mask = static.maskPadded([1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(i-1),[1:sectionInfo.height+2*sectionInfo.overlap]+sectionInfo.height*(j-1));
            
            sections{i,j}.static.windowed = static.padded([1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(i-1),[1:sectionInfo.height+2*sectionInfo.overlap]+sectionInfo.height*(j-1));
            sections{i,j}.moving.windowed = moving.padded([1:sectionInfo.width+2*sectionInfo.overlap]+sectionInfo.width*(i-1),[1:sectionInfo.height+2*sectionInfo.overlap]+sectionInfo.height*(j-1));


            useMask = 0;
            sections{i,j}.static.unmasked = sections{i,j}.static.windowed;
            sections{i,j}.moving.unmasked = sections{i,j}.moving.windowed;
            if useMask
                sections{i,j}.static.unmasked = sections{i,j}.static.windowed;
                sections{i,j}.moving.unmasked = sections{i,j}.moving.windowed;

                sections{i,j}.moving.windowed = sections{i,j}.moving.windowed .* sections{i,j}.moving.mask;
                tempImg = sections{i,j}.static.windowed .* sections{i,j}.static.mask;
                sections{i,j}.static.windowed = (tempImg - min(tempImg(:)))./(max(tempImg(:))-min(tempImg(:)));%rescale to [0,1]
            end
            
            hgramBuckets = 60;
            hgram = hist(sections{i,j}.static.windowed(:),hgramBuckets);

            %each section of moving image is hist-equalized to corresponding
            %static section
            %sections{i,j}.moving.windowed = histeq(sections{i,j}.moving.windowed,hgram);

            sections{i,j}.moving.ctlPoints = [];
            sections{i,j}.static.ctlPoints = [];
        end
    end

end

%% manual pre-registration

if needToPickPoints
    for i = 1:size(sections,1)
        for j = 1:size(sections,2)
            if singleSection(1) == 0 || (i == singleSection(1) && j == singleSection(2))
                if isfield(sections{i,j}.moving,'ctlPoints') && ~isempty(sections{i,j}.moving.ctlPoints)
                    [sections{i,j}.moving.ctlPoints,sections{i,j}.static.ctlPoints] = cpselect(sections{i,j}.moving.windowed,sections{i,j}.static.windowed,sections{i,j}.moving.ctlPoints,sections{i,j}.static.ctlPoints,'Wait',true);
                else
                    [sections{i,j}.moving.ctlPoints,sections{i,j}.static.ctlPoints] = cpselect(sections{i,j}.moving.windowed,sections{i,j}.static.windowed,'Wait',true);
                end
            end
        end
    end
end
    

%% non-rigid registration

if needToRegister
    
    pool = parpool;
    
    if singleSection(1) == 0
        sectionsReshaped = reshape(sections,size(sections,1)*size(sections,2),1);
    else
        sectionsReshaped = sections(singleSection(1),singleSection(2));
    end
    
    parfor n = 1:length(sectionsReshaped)

        Options = struct;
        Options.Similarity = 'sd';
        %sd: squared distance, gd: gradient distance, gc: gradient correlation

        finalGridSpacing = 8;
        Options.MaxRef = 3;%max number of grid refinements (not incl initial)
        Options.Spacing = finalGridSpacing*2^Options.MaxRef * [1 1];%initial grid spacing
        Options.Penalty = 1e-2;%default is 1e-3
        Options.Registration = 'NonRigid';
        Options.Verbose = 1;%default is 2 (max value)
        Options.MaskMoving = [];
        %Options.MaskStatic = [];
        %Options.MaskMoving = sectionsReshaped{n}.moving.mask;
        Options.MaskStatic = ones(size(sectionsReshaped{n}.static.windowed));%Must provide somethign here



        Options.Points1 = fliplr(sectionsReshaped{n}.moving.ctlPoints);%for Imoving
        Options.Points2 = fliplr(sectionsReshaped{n}.static.ctlPoints);%for Istatic
        Options.PStrength = ones(size(Options.Points1,1),1);%might want to check and make sure this makes it a bit stronger than elastic penalty


        [Ireg,O_trans,Spacing,Msect,Bsect,Fsect] = image_registration(sectionsReshaped{n}.moving.windowed,sectionsReshaped{n}.static.windowed,Options);

        %sectionsReshaped{n}.moving.registered = Ireg;
        sectionsReshaped{n}.moving.registered = bspline_transform(O_trans,sectionsReshaped{n}.moving.unmasked,Spacing,3);
        sectionsReshaped{n}.moving.O_trans = O_trans;
        sectionsReshaped{n}.moving.Spacing = Spacing;
    %    sectionsReshaped{n}.moving.M = Msect;
        sectionsReshaped{n}.moving.B = Bsect;
        sectionsReshaped{n}.moving.F = Fsect;
    end

    delete(gcp)
    
    if singleSection(1) == 0
        sections = reshape(sectionsReshaped,size(sections,1),size(sections,2));
    else
        sections(singleSection(1),singleSection(2)) = sectionsReshaped;
    end

end


if needToReconstruct
    
    registered = zeros(size(moving.windowed));
    B = zeros([size(moving.windowed) 2]);
    F = zeros([size(moving.windowed) 2]);

    for i = 1:size(sections,1)
        for j = 1:size(sections,2)
            inputFrameYvals = [1:sectionInfo.width]+sectionInfo.width*(i-1);
            inputFrameXvals = [1:sectionInfo.height]+sectionInfo.height*(j-1);
            registered(inputFrameYvals,inputFrameXvals) = sections{i,j}.moving.registered([1:sectionInfo.width]+sectionInfo.overlap,[1:sectionInfo.height]+sectionInfo.overlap);
            B(inputFrameYvals,inputFrameXvals,:)          = sections{i,j}.moving.B([1:sectionInfo.width]+sectionInfo.overlap,[1:sectionInfo.height]+sectionInfo.overlap,:);
            F(inputFrameYvals,inputFrameXvals,:)          = sections{i,j}.moving.F([1:sectionInfo.width]+sectionInfo.overlap,[1:sectionInfo.height]+sectionInfo.overlap,:);
            %subplot(2,2,(j-1)+2*(i-1)+1);imagesc(F(:,:,1).^2+F(:,:,2).^2);
        end

    end

    moving.registered = registered;
    moving.B = B;
    moving.F = F;

end

imgPair.moving = moving;
imgPair.static = static;

if exist('sections','var')
    imgPair.sections = sections;
    imgPair.sectionInfo = sectionInfo;
end

if needToReconstruct; transformQuiver(imgPair); end

function cmap = cmapblackandred(thresh,min,max)
%returns a 64x3 array coding a colormap that shows red for pixels below
%intensity value thresh, and grey scale for pixels equal or greater than
%thresh. min and max are the minimum and maximum intensity values in the
%image.


%cmap(i,:) specifies the RGB value of pixels in bucket i.
%There are 64 buckets
%Bucket i is pixels with intensity values in the range:
%xmin + (max-min)/64*[(i-1) : i] 


threshFractional = (thresh-min)/(max-min);
%this must be between 0 and 63/64

threshIndex = 1+floor(threshFractional*64);

cmap = double(zeros(64,3));

if threshIndex > 1
    cmap([1:threshIndex-1],1) = double(ones(threshIndex-1,1));
end

greyValues = linspace(0,1,64-(threshIndex-1));
cmap([threshIndex:64],1) = greyValues;
cmap([threshIndex:64],2) = greyValues;
cmap([threshIndex:64],3) = greyValues;

