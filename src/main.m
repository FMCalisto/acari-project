clear all, close all

% PERGUNTAS
%
% 2) Como e que faco a funcionalidade de vis. da trajectoria?
% 2.1) Pode ser as varias box ao longo do tempo?
% 3) O que e o num de deteccao de falhas?
% 4) Que key-frames sao estas?
%

imgbk = imread('../frames/SonofMated10/SonofMated1000262.jpg');

thr = 29; % Optimal Tested Value: 29
minArea = 7; % Optimal Tested Value: 7
baseNum = 262; % Initial Frame
seqLength = 23353;

se = strel('disk',3);

% -------------------- Backgroud -------------------- %

nFrameBKG= 2000; % 23354 Frames used to compute background image
step=20;       % Faz display de step em step frames

Bkg=zeros(size(imgbk));
BkgLast=zeros(size(imgbk));
alfa=0.05;
figure; hold on

%Exprimentar varios valores para ALPHA

for i = 0 : step : nFrameBKG
    sprintf('BKG %d',i)
    imgfr = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                   baseNum+i));
    Y = imgfr;
    Bkg = alfa * double(Y) + (1 - alfa) * double(Bkg);
    
    imgUInt8 = uint8(Bkg);
    imgUInt8Last = uint8(BkgLast);
    %imshow(imgUInt8); drawnow
    %imshow(imgUInt8Last); drawnow
    
end

% --------------------------------------------------- %

imgBkgBase = imgUInt8; % Imagem de background

% -------------------- ROI -------------------- %
% Remove object intersection
% Faz as caixinhas

stepRoi = 20;
nFrameROI = 23353;

for i = 0 : stepRoi : nFrameROI
    imgfrNew = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+i));
    
    sprintf('ROI %d',i);
    hold off
    imshow(imgfrNew); 
    
    %compare frame with background image
    imgdif = (abs(double(imgBkgBase(:,:,1))-double(imgfrNew(:,:,1)))>thr) | ...
        (abs(double(imgBkgBase(:,:,2))-double(imgfrNew(:,:,2)))>thr) | ...
        (abs(double(imgBkgBase(:,:,3))-double(imgfrNew(:,:,3)))>thr);
    
    
    bw = imclose(imgdif,se);
    str = sprintf('Frame: %d',i); title(str);
    %%%%%%%%imshow(bw); 
    [lb num]=bwlabel(bw);
    regionProps = regionprops(lb,'area','FilledImage','Centroid');
    inds = find([regionProps.Area]>minArea);
    
    regnum = length(inds);
    
    if regnum
        for j=1:regnum
            [lin col]= find(lb == inds(j));
            upLPoint = min([lin col]);
            dWindow  = max([lin col]) - upLPoint + 1;
           
            rectangle('Position', ...
                      [fliplr(upLPoint) fliplr(dWindow)], ...
                      'EdgeColor',[1 1 0], ...
                      'linewidth',2);
        end
    end
    drawnow
end

% --------------------------------------------------- %