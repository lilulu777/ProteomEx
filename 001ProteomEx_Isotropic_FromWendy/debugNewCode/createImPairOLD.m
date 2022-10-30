function [imPairOLD] = createImPairOLD(imgPair)

%convert imgPair (new format) into imPair (old format), including control
%point and mask info, all for imgPair.sections{1,1}.
%(old format does not have sections implemented)


%unwrap input
section = imgPair.sections{1,1};


%extract relevant data from sections{1,1} of input
ImgMoving.input = section.moving.windowed;
ImgStatic.input = section.static.windowed;

%skip rough alignment step (or execute to overwrite these values)
ImgMoving.windowed = section.moving.windowed;
ImgStatic.windowed = section.static.windowed;

%skip control point seelction step (or execute to overwrite these values)
ImgMoving.moving_out = section.moving.ctlPoints;
ImgStatic.fixed_out = section.static.ctlPoints;

%skip mask specification step (or execute to overwrite these values)
ImgMoving.mask = section.moving.mask;


%wrap up output
imPairOLD.ImgMoving = ImgMoving;
imPairOLD.ImgStatic = ImgStatic;
