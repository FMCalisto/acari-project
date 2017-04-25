function main
clear all, close all

% PERGUNTAS
%
% 2) Como e que faco a funcionalidade de vis. da trajectoria?
% 2.1) Pode ser as varias box ao longo do tempo?
% 3) O que e o num de deteccao de falhas?
% 4) Que key-frames sao estasasd?
%

imgbk = imread('../frames/SonofMated10/SonofMated1000262.jpg');

thr = 29; % Optimal Tested Value: 29
minArea = 7; % Optimal Tested Value: 7
baseNum = 262; % Initial Frame
seqLength = 23353;

se = strel('disk',3);

% -------------------- Backgroud -------------------- %

nFrameBKG= 1000; % 23354 Frames used to compute background image
step=20;       % Faz display de step em step frames

Bkg=zeros(size(imgbk));
BkgLast=zeros(size(imgbk));
alfa=0.05;
mainfigure = figure; hold on
%movegui(mainfigure, 'northwest');
%set(mainfigure, 'Position', [100, 1000, 100, 100])

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
    if i == 500
        %touch(500);
        %figure(mainfigure);
    end
    if i == 360
        %sex(350,360);
        %figure(mainfigure);
    end
end

% --------------------------------------------------- %

imgBkgBase = imgUInt8; % Imagem de background

% -------------------- ROI -------------------- %
% Remove object intersection
% Faz as caixinhas

stepRoi = 20;
nFrameROI = 700;  % 23354 Frames used to compute background image

for i = 0 : stepRoi : nFrameROI
    imgfrNew = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+i));
    
    sprintf('ROI %d',i);
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            %
    % Tentativa de RESIZE WINDOW %
    %                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %imageAuxRoi = imresize(imgfrNew, 0.5);lixo
    %set(mainfigure, 'Position', [100, 1000, 100, 100])lixo
    %imshow(imageAuxRoi);lixo
    
    imshow(imgfrNew); %% Caminho rectangulos amarelos - Background 
    
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
            if (200 < i)&&(i < 260)    
                disp('i');
                disp(i);
                disp('upLPoint');
                disp(upLPoint);
                disp('dWindow');
                disp(dWindow);
            end
            rectangle('Position', ...
                      [fliplr(upLPoint) fliplr(dWindow)], ...
                      'EdgeColor',[1 1 0], ...
                      'linewidth',2);
            
            %shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom',...
            %    'CustomBorderColor',yellow);

            %circles = int32([upLPoint(1,1) upLPoint(1,2) 20]); %  [x1 y1 radius1]
            %J = shapeInserter(RGB, circles);
            %drawnow
            
            % --------------------------------------------------- %
            %                                                     % 
            %                    Marker                           % 
            %                                                     % 
            % --------------------------------------------------- %
            
            %center = find([regionProps.Centroid]>0.5);
            
            center = regionProps.Centroid;
            
            H = insertMarker(imgfrNew,[center(1,1) center(1,2)]);
            %J = insertMarker(H,[dWindow(1,1) dWindow(1,2)]);
            if (200 < i)&&(i < 260)
                disp('center');
                disp(center);
            end
            drawnow
            imshow(H);
            drawnow
            
        end
    end
    drawnow
end

% --------------------------------------------------- %
%                                                     % 
%                    Touch                            % 
%                                                     % 
% --------------------------------------------------- %

end

function touch(n)
    baseNum = 262; % Initial Frame
    touchFigure = figure;
    movegui(touchFigure, 'northeast');
    %subplot(2,2,1,'align');
    str = sprintf('Touch: %d',n); title(str);
    touchImageT = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+n));
    imageAux = imresize(touchImageT, 0.4);
    imshow(imageAux); title(str);
    %touch=n;
    %drawnow;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Colocar Rectangulos     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% --------------------------------------------------- %
%                                                     % 
%                      Sex                            % 
%                                                     % 
% --------------------------------------------------- %

function sex(m,l)
    baseNum = 262; % Initial Frame
    sexFigure = figure;
    movegui(sexFigure, 'east');
    
    strSex = sprintf('Before Sex: %d',m);% title(strSex);
    [S1, map1] =  imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+m));
    %imageAuxSex = imresize(sexImageT, 0.4);
    
    subplot(1,2,1,'align'), imshow(S1, map1);  title(strSex);
    
    strSex2 = sprintf('After Sex: %d',l); title(strSex2);
    [S2, map2] = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+l));
    %imageAuxSex2 = imresize(sexImageT2, 0.4);
    
    subplot(1,2,2,'align'), imshow(S2, map2),  title(strSex2);
    %imshow(imageAux); title(str);
    %drawnow;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Colocar Rectangulos     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end