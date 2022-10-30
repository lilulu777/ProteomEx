% NonRigidReg_MasterScript
% Basic workflow for general projects in lab

% Pre-1st-script image processing
% (1) Prepare the pre and post image, run TurboReg on ImageJ (post to pre;
% "scaled rotation" mode) to rigidly register the images.
% (2) Output file is a 2-channel image.  Remove the 2nd channel.
% (3) Use "Merge channel" to overlap the pre and post images (order does 
% not matter)
% (4) Crop out a square ROI to be analyzed (i.e. an ROI that fully 
% composes of samples, i.e. leaving out the blank spaces that resulted 
% from image rotation, etc)
% (5) downsample by "Binning" (this step dramatically reduces computational
% time; we typically down-sample to < 2000 by 2000 pixels)
% (6) "Split channels" to break apart the pre and post images; save 
% seperately as TIF in a folder


%% Load Image

imgPair.moving.input = im2gray(imread('MT3 registered post.tif')); % set to post
imgPair.static.input = im2gray(imread('MT3 pre.tif')); % set to pre; this is ground truth
imgSize = size(imgPair.moving.input)

%% Run 1st Script & Save Workspace
imgPair = imageRegistrationWorkflowParfor(imgPair,imgSize,1);
currTime = clock;
timestamp = [num2str(currTime(1)*10000+ currTime(2) * 100 + currTime(3)) '-' ...
    num2str(currTime(4),'%.2d') '_' num2str(currTime(5),'%.2d')]
save(['Workspace_After_Workflow_' timestamp])

%% Set Parameters
imgPair.info.pixelWidth = 25;  % um
sampleInterval = 7; 
% this is the spacing between pixels where distortion is calculated; 
% larger values reduce the x-resolution of the output error plot, but
% significantly reduces computational time

%% Run 2nd Script & Save Workspace
autoSave = 0;
read_saved_workspace = 0;

imgPair = measurementError(imgPair,sampleInterval,...
            autoSave,read_saved_workspace);
currTime = clock;
timestamp = [num2str(currTime(1)*10000+ currTime(2) * 100 + currTime(3)) '-' ...
    num2str(currTime(4),'%.2d') '_' num2str(currTime(5),'%.2d')]
save(['Workspace_After_ErrorAnalysis_' timestamp])

