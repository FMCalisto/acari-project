clear all

% PERGUNTAS
%
% 1) Como retiro a funcionalidade de manter uma ROI no inicio das frames?
% 2) Como e que faco a funcionalidade de vis. da trajectoria?
% 2.1) Pode ser as varias box ao longo do tempo?
% 3) O que e o num de deteccao de falhas?
% 4) Que key-frames sao estas?
%

%imgbk = imread('../frames/SonofMated2/SonofMated200001.jpg');
imgbk = imread('../frames/SonofMated10/SonofMated1000262.jpg');
%imgbk = imread('../frames/SonofMated10/SonofMated1000000.jpg');

thr = 25;
minArea = 50;
baseNum = 262;
seqLength = 23353;

se = strel('disk',3);

% --------- cBkg2 -----------%

nFrame= 40*25;
step=20;

Bkg=zeros(size(imgbk));
BkgLast=zeros(size(imgbk));
alfa=0.05;
figure; hold on

%Exprimentar varios valores para ALPHA

for i = 0 : step : nFrame
    i
    imgfr = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', baseNum+i));
    imgfr2 = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', seqLength - i));
    Y = imgfr;
    Y2 = imgfr2;
    Bkg = alfa * double(Y) + (1 - alfa) * double(Bkg);
    BkgLast = alfa * double(Y2) + (1 - alfa) * double(BkgLast);
    
    imgUInt8 = uint8(Bkg);
    imgUInt8Last = uint8(BkgLast);
    %imshow(imgUInt8); drawnow
    imshow(imgUInt8Last); drawnow
    
    hold off
    %imshow(imgfr);
    
    imgdif = (abs(double(imgUInt8(:,:,1))-double(imgfr(:,:,1)))>thr) | ...
        (abs(double(imgUInt8(:,:,2))-double(imgfr(:,:,2)))>thr) | ...
        (abs(double(imgUInt8(:,:,3))-double(imgfr(:,:,3)))>thr);
    
    bw = imclose(imgdif,se);
    %imshow(bw)
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
    %drawnow
end
    
    
