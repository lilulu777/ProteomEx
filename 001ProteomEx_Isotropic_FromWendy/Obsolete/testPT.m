  %% Read two greyscale images of Lena
  imgDir = '/Users/paultillberg/Documents/MATLAB/Image registration/SIM before after from Fei 140109/';
  Iafter=imread([imgDir 'B5_3_3_MAX_After.tif']); 
  Ibefore=imread([imgDir 'B5_3_stack_MAX_before.tif']);
  %Scale both images to max pixel value for uint16
  Ibefore = uint16(double(Ibefore)/double(max(Ibefore(:))) * (2^16-2));
  Iafter = uint16(double(Iafter)/double(max(Iafter(:))) * (2^16-2));
  
  %Shift all values up by one and use value of zero to mean no data.
  Ibefore = Ibefore + 1;
  Iafter = Iafter + 1;


  %% User selects 2 control points to roughly register before to after image
  [moving_out0,fixed_out0] = cpselect(Ibefore,Iafter,'Wait',true);
  
  
  %% Register and before image to after image and pad out so they occupy same frame
  mytform = cp2tform(moving_out0,fixed_out0,'nonreflective similarity');
  
  sizeAfter = size(Iafter);
  IbeforeRough = imtransform(Ibefore,mytform,'XData',[1 sizeAfter(2)],'YData',[1 sizeAfter(1)],'FillValues',0);
  
  
  %% User selects many control points (at least 6 pairs)
  %IMPORTANT: Export points to cpstruct
  %moving_out = []; fixed_out = [];%%RESET POINTS?
  %[moving_out,fixed_out] = cpselect(IbeforeRough,Iafter,moving_out,fixed_out,'Wait',true);
  %tformManyPoints = cp2tform(moving_out,fixed_out,'lwm');
  %IbeforeManyPoints = imtransform(IbeforeRough,tformManyPoints,'XData',[1 sizeAfter(2)],'YData',[1 sizeAfter(1)],'FillValues',0);

  
  %% Pre-register, similarity only
  
%   IbeforeRough = uint16(double(IbeforeRough)/double(max(IbeforeRough(:))) * (2^16-1));
%   Iafter = uint16(double(Iafter)/double(max(Iafter(:))) * (2^16-1));
%   
%   optimizer = registration.optimizer.RegularStepGradientDescent();
%   metric = registration.metric.MeanSquares();
%   moving_reg = imregister(IbeforeRough,Iafter,'similarity',optimizer,metric)
  
  %% All pixels in IbeforeManyPoints with value of 0 are pixels in Iafter not
  %represented in Ibefore, so zero these out in Iafter as well
%  IafterWithZeros = Iafter .* uint16(boolean(IbeforeManyPoints));
  IafterWithZeros = Iafter .* uint16(boolean(IbeforeRough));
  
  %% Register the images
  %Treat IbeforeRough as static image
  npixX = 1000;
  npixY = 800;
  window.x = [1 npixX] + 0;
  window.y = [1 npixY] + 200;

%  window.x = [1 1344];
%  window.y = [1 1024];

  sections.width = 200;%Size of each section for local registration
  sections.overlap = 20;%Overlap to include between neighboring sections (discard registration results in overlap regions)
  if sections.overlap > sections.width
      error('Make sure section width is greater than section overlap')
  end
  


  IstaticFull = IbeforeRough(window.y(1):window.y(2),window.x(1):window.x(2));
  ImovingFull = IafterWithZeros(window.y(1):window.y(2),window.x(1):window.x(2));
  
  %% mask

sd = 8; %Use 8 pixel s.d. blurring
gFilter = fspecial('gaussian', sd*4, sd);
blurred = imfilter(ImovingFull,gFilter,'replicate');

figure; imagesc(blurred); colormap(cmapblackandred);%figure out threshold by hand
colormapeditor;

%% make mask
threshold = 8200;

sizeImg = size(blurred);
maskFull = uint16(zeros(sizeImg));

for i = 1:sizeImg(1)
    for j = 1:sizeImg(2)
        if blurred(i,j) < threshold; maskFull(i,j) = 0;
        else maskFull(i,j) = 1;
        end
    end
end

%% registration
  
  IstaticFull = uint16(double(IstaticFull)/double(max(IstaticFull(:))) * (2^16-1));
  ImovingFull = uint16(double(ImovingFull)/double(max(ImovingFull(:))) * (2^16-1));
  
  imgSize = size(IstaticFull);
  
  IstaticPadded = uint16(zeros(imgSize + 2*sections.overlap));
  IstaticPadded(sections.overlap+1:sections.overlap+imgSize(1),sections.overlap+1:sections.overlap+imgSize(2)) = IstaticFull;
  ImovingPadded = uint16(zeros(imgSize + 2*sections.overlap));
  ImovingPadded(sections.overlap+1:sections.overlap+imgSize(1),sections.overlap+1:sections.overlap+imgSize(2)) = ImovingFull;
  maskPadded = uint16(zeros(imgSize + 2*sections.overlap));
  maskPadded(sections.overlap+1:sections.overlap+imgSize(1),sections.overlap+1:sections.overlap+imgSize(2)) = maskFull;
  
  
  %Specify subfigure windows as arrays of x and y values to serve as [1 1]
  %pixel for each section. Pixels positions are referenced to original
  %windowed images, not zero-padded images
  sections.xvals = [1:sections.width:npixX];
  sections.yvals = [1:sections.width:npixY];
  %Include overlaps. Images will be padded out on all sides with zeros to
  %include overlap area in all sections
  
  
  sections.xmin = sections.xvals - sections.overlap;
  sections.xmax = circshift(sections.xvals,[0 -1])-1 + sections.overlap;
  sections.xmax(end) = window.x(2) + sections.overlap;
  sections.ymin = sections.yvals - sections.overlap;
  sections.ymax = circshift(sections.yvals,[0 -1])-1 + sections.overlap;
  sections.ymax(end) = window.y(2) + sections.overlap;
  
  %Shift windows to get pixel positions in coordinate system of zero-padded
  %images
  sections.xvals = sections.xvals + sections.overlap;
  sections.yvals = sections.yvals + sections.overlap;
  sections.xmin = sections.xmin + sections.overlap;
  sections.ymin = sections.ymin + sections.overlap;
  sections.xmax = sections.xmax + sections.overlap;
  sections.ymax = sections.ymax + sections.overlap;
  
  clear Results
  
  Results{length(sections.ymin),length(sections.xmin)} = [];

  display('computing local transformations')
  display('currently on section:')
  for i = 1:length(sections.ymin)
      for j = 1:length(sections.xmin)
          display(['i=' num2str(i) ' j=' num2str(j)])
          
          secSourceYvals = [1:sections.width + 2*sections.overlap] + sections.width*(i-1);
          secSourceXvals = [1:sections.width + 2*sections.overlap] + sections.width*(j-1);

          
          Istatic = IstaticPadded(secSourceYvals,secSourceXvals);
          Imoving = ImovingPadded(secSourceYvals,secSourceXvals);
          mask = maskPadded(secSourceYvals,secSourceXvals);

          hgramBuckets = 60;
          hgram = hist(double(Istatic(:)),hgramBuckets);
          Imoving = histeq(Imoving,hgram);

          Options.Penalty = 1e-2;%default is 1e-3
        %  Options.Spacing = [32 32];
          Options.Registration = 'NonRigid';
          Options.Verbose = 1;%default is 2 (max value)
          Options.MaskMoving = [];
          %Options.MaskMoving = mask;
          [Ireg,O_trans,Spacing,M,B,F] = image_registration(Imoving,Istatic,Options);
          Results{i,j}.Istatic = Istatic;
          Results{i,j}.Imoving = Imoving;
          Results{i,j}.Ireg = Ireg;
          Results{i,j}.O_trans = O_trans;
          Results{i,j}.Spacing = Spacing;
          Results{i,j}.M = M;
          Results{i,j}.B = B;
          Results{i,j}.F = F;
      end
  end

  
  %% 
  
  imgSize = size(IstaticFull);
  
  %size of image before padding
  Ireg = zeros(imgSize);
  M = zeros([imgSize 2]);
  B = zeros([imgSize 2]);
  F = zeros([imgSize 2]);
  
  secSrceXvals = [1:sections.width]+sections.overlap;
  secSrceYvals = [1:sections.width]+sections.overlap;
    
  for i = 1:length(sections.ymin)
      for j = 1:length(sections.xmin)
          secTargXvals = [1:sections.width] + sections.width*(j-1);
          secTargYvals = [1:sections.width] + sections.width*(i-1);
          Ireg(secTargYvals,secTargXvals) = Results{i,j}.Ireg(secSrceYvals,secSrceXvals);
%          Istatic(secTargXvals,secTargYvals) = Results{i,j}.Istatic(secSrceXvals,secSrceYvals);
%          Imoving(secTargXvals,secTargYvals) = Results{i,j}.Imoving(secSrceXvals,secSrceYvals);
          B(secTargYvals,secTargXvals,:) = Results{i,j}.B(secSrceYvals,secSrceXvals,:);
          F(secTargYvals,secTargXvals,:) = Results{i,j}.F(secSrceYvals,secSrceXvals,:);

%          M(sections.yvals(i):sections.yvals(i) + sections.width-1,sections.xvals(j):sections.xvals(j) + sections.width-1,:) = Results{i,j}.M;
%Skipping affine

      end
  end
  

  % Show the registration result
%  figure,
%  subplot(2,2,1), imshow(Imoving); title('moving image');
%  subplot(2,2,2), imshow(Istatic); title('static image');
%  subplot(2,2,3), imshow(Ireg); title('registerd moving image');
  % Show also the static image transformed to the moving image
%  Ireg2=movepixels(Istatic,F);
%  subplot(2,2,4), imshow(Ireg2); title('registerd static image');

 % Show the transformation fields
%  figure,
%  subplot(2,2,1), imshow(B(:,:,1),[]), colorbar; title('Backward Transf. in y direction');
%  subplot(2,2,2), imshow(F(:,:,1),[]), colorbar; title('Forward Transf. in y direction');
%  subplot(2,2,3), imshow(B(:,:,2),[]), colorbar; title('Backward Transf. in x direction');
%  subplot(2,2,4), imshow(F(:,:,2),[]), colorbar; title('Forward Transf. in x direction');
  %PT-changed x and y, B(:,:,) stores transform in horizontal direction,
  %makes more sense to call that x.
  
  
  transformQuiver
  