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

%touchFigure = figure(2);
mainFigure = figure(1);

% --------------------- END Figure --------------------- %

% ---------------------- Message ----------------------- %

stringTouchAction = 'Action: TOUCH';
stringTouchSeconds = 'Touch Seconds: ';

stringCoupleAction = 'Action: COUPLE';
stringCoupleSeconds = 'Couple Seconds: ';

% --------------------- END Message -------------------- %

%set(touchFigure, 'Position', [630, 170, 500, 500]);
set(mainFigure, 'Position', [100, 000, 500, 1000]);
hold on
maleTrail = [];
femaleTrail = [];
touchDistArr = [];

% Normal Frames %
%baseNum = 262; % Initial Frame: 262
%nTotalFrames = 5000; % Total: 23354 Frames

% Coupling Frames %
baseNum = 15500; % Couple Frame: 15500
nTotalFrames = 5000; % Total: 23354 Frames

se = strel('disk',3);

        
%Exprimentar varios valores para ALPHA

for i = 0 : step : nFrameBKG
    %sprintf('BKG %d',i)
    imgfr = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                   baseNum + i));
    Y = imgfr;
    Bkg = alfa * double(Y) + (1 - alfa) * double(Bkg);
    
    imgUInt8 = uint8(Bkg);
    %imshow(imgUInt8); drawnow
    %imshow(imgUInt8Last); drawnow
%     if i == 500
%         touch(500,touchFigure);   %%%%%  TOUCH  %%%%%
%         figure(mainFigure);
%     end
%     if i == 360
%         sex(350,360, touchFigure);
%         figure(mainFigure);   %%%%%  SEX  %%%%%
%     end
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
    
    sprintf('ROI %d',i);
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
    
    %compare frame with background image
    %imgBkgBaseResize = imresize(imgBkgBase, 1);%resize back ground image
    imgdif = (abs(double(imgBkgBase(:,:,1))-double(imgfrNew(:,:,1)))>thr) | ...
        (abs(double(imgBkgBase(:,:,2))-double(imgfrNew(:,:,2)))>thr) | ...
        (abs(double(imgBkgBase(:,:,3))-double(imgfrNew(:,:,3)))>thr);
    
    
    bw = imclose(imgdif,se);
    str = sprintf('Frame: %d',i); title(str);
    % ----------------------------------------------------------- %
    %%%%%%%%imshow(bw);  %%Mete Background preto ao mesmo tempo
    % ----------------------------------------------------------- %
    [lb num]=bwlabel(bw);
    regionProps = regionprops(lb,'area','FilledImage','Centroid');
    inds = find([regionProps.Area]>minArea);
    
    regnum = length(inds);
    %%%%%%%%%%inds number explained
    %if regnum = 1 ha' 1 regiao, logo acasalamento
    %inds(1) - acaro1 e acaro2
    %if regnum = 2 ha' 2 regioes
    %inds(1) - acaro1
    %inds(2) - acaro2
    sizeTrail = size(maleTrail);
    for k = sizeTrail-20 : 1 : sizeTrail
        sizeTrail = sizeTrail(1,1);
        if i > 0
            if (sizeTrail < 21)
                plot(maleTrail(:,1),maleTrail(:,2),'*','Color', ...
                    'red', 'LineStyle','-');
                plot(femaleTrail(:,1),femaleTrail(:,2),'*','Color', ...
                    'blue','LineStyle','-');
            else
                var = k; % Var is the localization where we start using
                         % things.
                plot(maleTrail(var, 1),maleTrail(var, 2),'*', ...
                    'Color', 'red','LineStyle','-');
                plot(femaleTrail(var, 1),femaleTrail(var, 2),'*', ...
                    'Color', 'blue','LineStyle','-');         
                
            end
        end
    end
    
    trailNorm = norm(femaleTrail - maleTrail);
    %D = pdist2(femaleTrail(:, 1), maleTrail(:, 2));
    %disp('Pairwise Distance: ');
    %disp(D);
%     disp('Male Position: ');
%     disp(maleTrail(:, 1));
%     disp('regnum');
%     disp(regnum);
    if regnum
        regionProps(inds(1),1).Area
        
        %existem as 2 regioes
        
        if ( regnum > 1 )
            acariA = regionProps(inds(1));
            acariB = regionProps(inds(2));
%             disp(acariA.Centroid);
%             disp(acariB);
            if (acariA.Area) > (acariB.Area)  %se A > B, A e' femea
                femaleTrail = [femaleTrail;[acariA.Centroid(1,1) acariA.Centroid(1,2)]];
                maleTrail = [maleTrail;[acariB.Centroid(1,1) acariB.Centroid(1,2)]];
            else  %se A <= B, A e' macho
                maleTrail = [maleTrail;[acariA.Centroid(1,1) acariA.Centroid(1,2)]];
                femaleTrail = [femaleTrail;[acariB.Centroid(1,1) acariB.Centroid(1,2)]];
            end    
        else
            %existem 1 regiao
            if (regnum > 0)
                acariA = regionProps(inds(1));
                femaleTrail = [femaleTrail;[acariA.Centroid(1,1) acariA.Centroid(1,2)]];
                maleTrail = [maleTrail;[acariA.Centroid(1,1) acariA.Centroid(1,2)]];
            end
        %caso nao existao 2 regioes nao adiciona nada    
        end
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
        
        % ------------- Male/Female Position ------------- %
                
%         disp('Male Position: ');
%         disp(maleTrail(sizeTrail, 1));
%         disp('Female Position: ');
%         disp(femaleTrail(sizeTrail, 2));
        sizeMaleTrail = size(maleTrail);
        DX = [femaleTrail(sizeMaleTrail(1,1), 1), ...
              femaleTrail(sizeMaleTrail(1,1), 2); ...
              maleTrail(sizeMaleTrail(1,1), 1), ...
              maleTrail(sizeMaleTrail(1,1), 2)];
        D = pdist(DX, 'euclidean');
        disp('Male/Female Distance: ');
        disp(D);

        touchDistArr = [touchDistArr; D];

        disp(touchDistArr)

        %touchDistance = D < 10;

        couplingDistance = D < 5;
        couplingTime = 600 / i < 1;
        isCoupling = couplingDistance && couplingTime;

        sizeTouchDistArr = size(touchDistArr);

        if (touchDistArr(sizeTouchDistArr(1, 1), 1) < 10)
            if (isCoupling)
                timeStartCoupling = num2str(frameToTime(i));
                timerAnnotation(timeStartCoupling, stringCoupleSeconds);
                actionAnnotation(stringCoupleAction);
            else
                timeStartTouch = num2str(frameToTime(i));
                timerAnnotation(timeStartTouch, stringTouchSeconds);
                actionAnnotation(stringTouchAction);
            end
        end
        
    end
    
        % --------------------------------------------------- %
        %   Grafic    correr depois de criar distancia        % 
        % --------------------------------------------------- %
    
%     warning off;
% 
%     h = animatedline('Marker','v','Color','red','LineStyle','-');
%     axis([0,4*pi,-1,1])
% 
%     x = linspace(0,4*pi,1000);
%     y = sin(x);
% 
%     title('Graph of Sine and Cosine Between -2\pi and 2\pi');
%     xlabel('-2\pi < x < 2\pi'); % x-axis label
%     ylabel('sine and cosine values'); % y-axis label
%     legend('y = sin(x)','y = cos(x)','Location','southwest');
% 
%     %%get xdata is distanciaEntreAcaros and ydata is time to grafic
%     %% [xdata,ydata] = getpoints(h);
%      
%     for p = 1:length(x)
%         addpoints(h,x(p),y(p));
%         drawnow
%     end
%     %clearpoints(h) %Limpar Points

    drawnow
    clf(mainFigure, 'reset');
    
end

% --------------------------------------------------- %
%                                                     % 
%                    Touch                            % 
%                                                     % 
% --------------------------------------------------- %

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