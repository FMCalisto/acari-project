function main
clear
clc
close all

% PERGUNTAS
%
% 2) Como e que faco a funcionalidade de vis. da trajectoria?
% 2.1) Pode ser as varias box ao longo do tempo?
% 3) O que e o num de deteccao de falhas?
% 4) Que key-frames sao estasasd?
%

% -------------------- Timer -------------------- %

timerCopula = 0;

% -------------------- END Timer -------------------- %

% -------------------- Backgroud -------------------- %

imgbk = imread('../frames/SonofMated10/SonofMated1000262.jpg');

thr = 29; % Optimal Tested Value: 29
minArea = 50; % Optimal Tested Value: 7
alfa=0.10;
baseBkg = 262 % Initial Frame: 262

nFrameBKG = 1000; % 23354 Frames used to compute background image
step = 20;       % Faz display de step em step frames
Bkg=zeros(size(imgbk));

% -------------------- END Backgroud -------------------- %

% ----------------------- Figure ------------------------ %

touchFigure = figure(2);
mainFigure = figure(1);

%set(touchBox, 'String', 'NO TOUCH')
%set(touchBox, 'Position', [40,10,100,25])
set(touchFigure, 'Position', [630, 170, 500, 500]);
set(mainFigure, 'Position', [100, 000, 500, 1000]);
hold on

% --------------------- END Figure --------------------- %

% ---------------------- Message ----------------------- %

stringTouchAction = 'Action: TOUCH';
stringTouchSeconds = 'Touch Seconds: ';

stringCoupleAction = 'Action: COUPLE';
stringCoupleSeconds = 'Couple Seconds: ';

stringTotalLengthMale = 0;
stringTotalLengthFemale = 0;

stringMaleDistance = 'Male Distance: ';
stringFemaleDistance = 'Female Distance: ';

% --------------------- END Message -------------------- %

frameFirstCoupla = 0;
countTouch = 0;
countCoupla = 0;
counterAux = 0;
counterAuxCopula = 1;

%maleIndex = 0;
%femaleIndex = 0;

firstTouch = 1;

maleTrail = [];
femaleTrail = [];
touchDistArr = [];
allPosCopula =[];
distMaleFemale = [];

maleDistance = 0;
femaleDistance = 0;
distEuclidean = 0;

% Normal Frames %
%baseNum = 262; % Initial Frame: 262
%nTotalFrames = 5000; % Total: 23354 Frames

% Coupling Frames %
baseNum = 15500; % Couple Frame: 15500
nTotalFrames = 5000; % Total: 23354 Frames

se = strel('disk',3);

        
% Try to fix the alfa values...

for i = 0 : step : nFrameBKG
    imgfr = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                   baseNum + i));
    Y = imgfr;
    Bkg = alfa * double(Y) + (1 - alfa) * double(Bkg);
    
    imgUInt8 = uint8(Bkg);
end

% ------------------ END Backgroud ------------------ %

% --------------------------------------------------- %

imgBkgBase = imgUInt8; % Imagem de background

% -------------------- ROI -------------------------- %
% Remove object intersection
% Faz as caixinhas

stepRoi = 20;
nFrameROI = nTotalFrames;  % 23354 Frames used to compute background image

for i = 0 : stepRoi : nFrameROI

    imgfrNew = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum + i));
    
    sprintf('ROI %d', i);
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                            %
    % Tentativa de RESIZE WINDOW %
    %                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %imageResizeRoi = imresize(imgfrNew, 1);%resize mainFigure
%     set(mainFigure, 'Position', [100, 200, 500, 500])
    %imshow(imageAuxRoi);lixo
    
    % --------------------------------------------------- %
    %%%%%subplot(2,1,1,'align');  %put analized image on the left side of mainFigure
    imshow(imgfrNew); %% Caminho rectangulos amarelos - Background 
    hold on
    % --------------------------------------------------- %
    
    % Compare frame with background image:
    imgdif = (abs(double(imgBkgBase(:,:,1))-double(imgfrNew(:,:,1)))>thr) | ...
        (abs(double(imgBkgBase(:,:,2))-double(imgfrNew(:,:,2)))>thr) | ...
        (abs(double(imgBkgBase(:,:,3))-double(imgfrNew(:,:,3)))>thr);
    
    
    bw = imclose(imgdif, se);
    str = sprintf('Frame: %d', i); title(str);
    
    [lb num]=bwlabel(bw);
    regionProps = regionprops(lb, 'area', 'FilledImage', 'Centroid');
    inds = find([regionProps.Area] > minArea);
    
    regnum = length(inds);
    
    if (regnum > 1)
        regionPropsArr1 = [regionProps(inds(1), 1).Area];
        regionPropsArr2 = [regionProps(inds(2), 1).Area];
        if (regionPropsArr1 < regionPropsArr2)
            maleIndex = 1;
            areaMale = regionPropsArr1;
            femaleIndex = 2;
            areaFemale = regionPropsArr2;
        else
            maleIndex = 2;
            areaMale = regionPropsArr2;
            femaleIndex = 1;
            areaFemale = regionPropsArr1;
        end
        
        xFemale = regionProps(inds(femaleIndex), 1).Centroid(1);
        yFemale = regionProps(inds(femaleIndex), 1).Centroid(2);
        xMale = regionProps(inds(maleIndex), 1).Centroid(1);
        yMale = regionProps(inds(maleIndex), 1).Centroid(2);

        centroidMale = regionProps(inds(maleIndex), 1).Centroid;
        centroidFemale = regionProps(inds(femaleIndex), 1).Centroid;
        centroidCopola = regionProps(inds(1), 1).Centroid;
        centroidMaleArr = [centroidMale];
        centroidFemaleArr = [centroidFemale];
        centroidCopulaArr = [centroidCopola];

        plusPos = sqrt((xFemale-xMale).^2) + ((yFemale-yMale).^2);
        distMaleFemale = [distMaleFemale; plusPos];
        
        maleTrail = [maleTrail; centroidMaleArr];
        femaleTrail = [femaleTrail; centroidFemaleArr];
        
        if i > i
            plot(maleTrail(counterAux, 1), ...
                 maleTrail(counterAux, 2), ...
                 '*', 'Color', 'red', 'LineStyle', '-');
            plot(femaleTrail(counterAux, 1), ...
                 femaleTrail(counterAux, 2), ...
                 '*', 'Color', 'blue', 'LineStyle', '-');
        end
        counterAux = counterAux + 1;
    else
        % TOUCH DETECTION %
        % If the number of active regions is 1
        allPosCopula = [allPosCopula; centroidCopulaArr];
        touchDistArr = [touchDistArr; i];
        distMaleFemale = [distMaleFemale; 0];
        
        if (firstTouch)
            countTouch = countTouch + 1;
            firstTouch = false;
        end
        
        if (countTouch > 1)
            frameFirstCopula = i;
            timerCopula = timerCopula + 1;
        end
        
        if i > 1
            plot(allPosCopula(counterAuxCopula, 1), ...
                 allPosCopula(counterAuxCopula, 2), ...
                 '*', 'Color', 'magenta', 'LineStyle', '-');
        end
        
        if regnum
            for j=1:regnum
                [lin col]= find(lb == inds(j));
                upLPoint = min([lin col]);
                dWindow  = max([lin col]) - upLPoint + 1;

                rectangle('Position', ...
                          [fliplr(upLPoint) fliplr(dWindow)], ...
                          'EdgeColor',[1 1 0], ...
                          'linewidth',2);


                drawnow

            end
        end
        
        for i = 1 : size(maleTrail, 1)
            plot(maleTrail(counterAux, 1), ...
                 maleTrail(counterAux, 2), ...
                 '*', 'Color', 'red', 'LineStyle', '-');
        end
        for i = 1 : size(femaleTrail, 1)
            plot(femaleTrail(counterAux, 1), ...
                 femaleTrail(counterAux, 2), ...
                 '*', 'Color', 'blue', 'LineStyle', '-');
        end    
        for i = 1 : size(allPosCopula, 1)
            plot(allPosCopula(counterAuxCopula, 1), ...
                 allPosCopula(counterAuxCopula, 2), ...
                 '*', 'Color', 'magenta', 'LineStyle', '-');
        end
        
        subplot(2,1,2), plot(distMaleFemale), title('Distance between Male and Female Tetranychus'), xlabel('Frames (x19)'), ylabel('Distance (px)'),
        subplot(2,1,1), drawnow
            
    end

    drawnow
    clf(mainFigure, 'reset');
    
end

% --------------------------------------------------- %
%                                                     % 
%                 Output Data                         % 
%                                                     % 
% --------------------------------------------------- %

% Create the figure
mFigure = figure('Name','Output Data')


ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[.3 .1 .6 .8]);

t = 0:1000;
y = 0.25*exp(-0.005*t);
plot(ax2,t,y)

title('Resumo de informacoes do Video:')
descr = {'Distancia percorrida'
    'macho:';
    'Distancia percorrida'
    'femea: ';
    ' ';
    'N? de toques';
    'N? de copulas';
    'Frame e '
    'tempo do 1? toque'
    'Frame e '
    'tempo do 1? copula';
    };

axes(ax1) % sets ax1 to current axes
text(.025,0.6,descr)

% ------------------------------------------------------- %
%    EUCLIDEAN DISTANCE CALCULATION FOR FEMALE vs MALE    %
% ------------------------------------------------------- %

allPosCopulaLenght = length(allPosCopula) - 1;
femaleTrailLenght = length(femaleTrail) - 1;
maleTrailLenght = length(maleTrail) - 1;

for k = 1: allPosCopulaLenght
    x1 = allPosCopula(k,1);
    y1 = allPosCopula(k,2);
    x2 = allPosCopula(k + 1,1);
    y2 = allPosCopula(k + 1,2);
    
    dist = sqrt((x2-x1).^2) + ((y2-y1).^2);
    copulaDistance = copulaDistance + dist;
end

for k = 1: femaleTrailLenght
    x1 = femaleTrail(k,1);
    y1 = femaleTrail(k,2);
    x2 = femaleTrail(k + 1,1);
    y2 = femaleTrail(k + 1,2);
    
    dist = sqrt((x2-x1).^2) + ((y2-y1).^2);
    femaleDistance = femaleDistance + dist;
end

for k = 1: maleTrailLenght
    x1 = maleTrail(k,1);
    y1 = maleTrail(k,2);
    x2 = maleTrail(k + 1,1);
    y2 = maleTrail(k + 1,2);
    
    dist = sqrt((x2-x1).^2) + ((y2-y1).^2);
    maleDistance = maleDistance + dist;
end

maleDistance = maleDistance + copulaDistance;
dist('-----------> Dist: ');
dist(maleDistance);
femaleDistance = femaleDistance + copulaDistance;

% SummaryBox1 = uicontrol('style','text');
% str = {'Distancia percorrida pelo macho (px): ',num2str(maleDistance)};
% set(SummaryBox1,'String',str)
% set(SummaryBox1,'Position',[40,100,200,50])
% 
% SummaryBox2 = uicontrol('style','text');
% str = {'Distancia percorrida pela femea (px): ',num2str(femaleDistance)};
% set(SummaryBox2,'String',str)
% set(SummaryBox2,'Position',[40,150,200,50])
% 
% SummaryBox3 = uicontrol('style','text');
% str = {'Frames at? a primeira copula : ',num2str(firstCopulaFrame)};
% set(SummaryBox3,'String',str)
% set(SummaryBox3,'Position',[40,200,200,50])

% ------------------------------------------------------- %
%                          END                            %
% ------------------------------------------------------- %

end

function timerAnnotation(time, message)
    dim1 = [.2 .5 .3 .3];
    stringTouchSecondsLabel = strcat(message, ' ', time);
    annotation('textbox', dim1, 'String', stringTouchSecondsLabel, ...
                      'FitBoxToText', 'on');
end

function actionAnnotation(action)
    dim2 = [.2 .6 .3 .3];
    annotation('textbox', dim2, 'String', action, ...
                      'FitBoxToText', 'on');
end

function distanceAnnotation(msgMale, msgFemale, distMale, distFemale)
    dim3 = [.2 .4 .3 .3];
    dim4 = [.2 .3 .3 .3];
    stringMessageMaleLabel = strcat(msgMale, ' ', distMale);
    annotation('textbox', dim3, 'String', stringMessageMaleLabel, ...
                      'FitBoxToText', 'on');
    stringMessageFemaleLabel = strcat(msgFemale, ' ', distFemale);
    annotation('textbox', dim4, 'String', stringMessageFemaleLabel, ...
                      'FitBoxToText', 'on');
end

function actualTime = frameToTime(nFrame)
    frameTime = 788 / 7885;
    actualTime = (nFrame) * frameTime;
end

function touch(n,fig)
    baseNum = 262; % Initial Frame
    %touchFigure = figure(2);
%     movegui(fig, 'northeast');%mudar para set coordinates
    figure(fig); hold on
% insert touch image on figure
    subplot(2,2,1);
    str = sprintf('Touch: %d',n); title(str);
    touchImageT = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+n));
    imageAux = imresize(touchImageT, 0.4);
    imshow(imageAux); title(str);
    %touch=n;
    %drawnow;
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  !!!Colocar Rectangulos!!!  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% --------------------------------------------------- %
%                                                     % 
%                      Sex                            % 
%                                                     % 
% --------------------------------------------------- %

function sex(m,l,fig)

    baseNum = 262; % Initial Frame
    movegui(fig, 'northeast');%mudar para set coordinates
    figure(fig); hold on
% insert before sex image on figure
    subplot(2,2,1);
    strSexB = sprintf('Before Sex: %d',m);title(strSexB);    
    sexImageB = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+m));
    imageAuxB = imresize(sexImageB, 0.4);
    imshow(imageAuxB); title(strSexB);
    hold on
% insert after sex image on figure
    subplot(2,2,2);
    strSexA = sprintf('After Sex: %d',l);title(strSexA);
    sexImageA = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                      baseNum+l));
    imageAuxA = imresize(sexImageA, 0.4);
    imshow(imageAuxA); title(strSexA);
    hold on
%     drawnow;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  !!!Colocar Rectangulos!!!  % fazer o try no ROI e nao no criar BKG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
