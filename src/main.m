function main
clear
clc
close all
%     warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Se nessario volto a apagar estas linhas 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obj = VideoReader('../videos/SonofMated10.avi');
% 
% disp('info');
% info = get(obj);
% disp(info);
% 
% disp('imagebackground ');
% imagebackground = read(obj,1);
% teste = figure(100);
% imshow(imagebackground);
% 
% nFrames = obj.FrameRate;
% str = sprintf('nFrames: %d',round(nFrames));
% disp(str);
% disp('CurrentTime');
% CurrentTime = obj.CurrentTime;
% disp(CurrentTime);
%  
% disp('Duration');
% Duration = obj.Duration;
% disp(Duration);
% 
% disp('Height');
% Height = obj.Height;
% disp(Height);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        EXAMPLE 1: Normal Frames                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imgbk = imread('../frames/SonofMated10/SonofMated1000262.jpg');
%baseBkg = 262; % Initial Frame: 262
%baseNum = 262; % Initial Frame: 262
%nTotalFrames = 5000; % Total: 23354 Frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        EXAMPLE 2: Coupling Frames              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%imgbk = imread('../frames/SonofMated10/SonofMated1000262.jpg');
%baseBkg = 262; % Initial Frame: 262
%baseNum = 15500; % Couple Frame: 15500
%nTotalFrames = 5000; % Total: 23354 Frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        EXAMPLE 3: New Frames                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgbk = imread('../newframes/frame0000.jpg');
baseBkg = 0; % Initial Frame: 0
baseNum = 0;
nTotalFrames = 7885; % Total: 7885 Frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------- Timer -------------------- %

% ----------------- END Timer ------------------ %

% ----------------- Backgroud ------------------ %

thr = 29; % Optimal Tested Value: 29
minArea = 50; % Optimal Tested Value: 7
alfa=0.10;  %Exprimentar varios valores para ALPHA

nFrameBKG = 1000; % 23354 Frames used to compute background image
step = 20;       % Faz display de step em step frames
Bkg=zeros(size(imgbk));

% -------------------- END Backgroud -------------------- %

% ----------------------- Figure ------------------------ %

touchFigure = figure(2);
mainFigure = figure(1);

% --------------------- END Figure --------------------- %

% ---------------------- Message ----------------------- %

stringTouchAction = 'Action: TOUCH';
stringTouchSeconds = 'Touch Seconds: ';

stringCoupleAction = 'Action: COUPLE';
stringCoupleSeconds = 'Couple Seconds: ';

stringTotalLengthMale = 0;
stringTotalLengthFemale = 0;

% --------------------- END Message -------------------- %

frameFirstCouple = 0;
count = 0;

set(touchFigure, 'Position', [630, 170, 500, 500]);
set(mainFigure, 'Position', [100, 000, 500, 1000]);
hold on
maleTrail = [];
femaleTrail = [];
touchDistArr = [];

numKeyFrames = 0;

isCoupling = false;
isTouching = false;
isTouchingNow = false;
isCouplingNow = false;

se = strel('disk',3);

for i = 0 : step : nFrameBKG
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 1 & 2                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     imgfr = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
%                    baseNum + i));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 3: New Frames                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    imgfr = imread(sprintf('../newframes/frame%.4d.jpg', ...
                    baseNum + i));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 1 & 2                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     imgfrNew = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
%                    baseNum + i));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 3: New Frames                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    imgfrNew = imread(sprintf('../newframes/frame%.4d.jpg', ...
                    baseNum + i));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        isCouplingNow = couplingDistance && couplingTime;
        sizeTouchDistArr = size(touchDistArr);
        
        if(touchDistArr(sizeTouchDistArr(1, 1), 1) < 10)
            isTouchingNow = true;
        else
            isTouchingNow = false;
        end
        
        if(isTouching == false)
            if(isTouchingNow == true)
                touch(i, touchFigure, numKeyFrames);   %%%%%  TOUCH  %%%%%
                figure(mainFigure);
                
                numKeyFrames = numKeyFrames + 1;
                isTouching = isTouchingNow;
            end
        end
        
        if(isTouching == true)
            if(isTouchingNow == false)
                isTouching = isTouchingNow;
            end
        end
        
        if(isCoupling == false)
            if(isCouplingNow == true)
                sex(i, 'beforeSex', touchFigure, numKeyFrames);
                figure(mainFigure);   %%%%%  SEX  %%%%%
                
                numKeyFrames = numKeyFrames + 1;
                isCoupling = isCouplingNow;
            end
        end
        
        if(isCoupling == true)
            if(isCouplingNow == false)
                sex(i, 'afterSex', touchFigure, numKeyFrames);
                figure(mainFigure);   %%%%%  SEX  %%%%%
                
                numKeyFrames = numKeyFrames + 1;
                isCoupling = isCouplingNow;
            end
        end

        if (isTouchingNow)
            if (isCoupling)
                while(count == 0)
                    frameFirstCouple = i;
                    count = 1;
                end
                timeStartCoupling = num2str(frameToTime(i - frameFirstCouple));
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

function touch(n,fig, numKeyFrames)
    numKeyFrames = numKeyFrames + 1;
    baseNum = 262; % Initial Frame
    %touchFigure = figure(2);
%     movegui(fig, 'northeast');%mudar para set coordinates
    figure(fig); hold on
% insert touch image on figure
    subplot(3, 3, numKeyFrames);
    str = sprintf('Touch: %d',n); title(str);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 1 & 2                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     touchImageT = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
%                    baseNum + i));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 3: New Frames                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    touchImageT = imread(sprintf('../newframes/frame%.4d.jpg', ...
                    baseNum + i));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function sex(m, message, fig, numKeyFrames)

    numKeyFrames = numKeyFrames + 1;
    baseNum = 262; % Initial Frame
    movegui(fig, 'northeast');%mudar para set coordinates
    figure(fig); hold on
% insert  sex image on figure
    subplot(3, 3, numKeyFrames);
    strSexB = sprintf(message, ': %d', m);title(strSexB);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 1 & 2                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     sexImageB = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
%                    baseNum + m));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        EXAMPLE 3: New Frames                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sexImageB = imread(sprintf('../newframes/frame%.4d.jpg', ...
                    baseNum + m));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imageAuxB = imresize(sexImageB, 0.4);
    imshow(imageAuxB); title(strSexB);
    hold on
% % insert after sex image on figure
%     subplot(2, 2, numKeyFrames + 1);
%     strSexA = sprintf('After Sex: %d',l);title(strSexA);    
%     sexImageA = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
%                       baseNum+l));
%     imageAuxA = imresize(sexImageA, 0.4);
%     imshow(imageAuxA); title(strSexA);
%     hold on
%     drawnow;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  !!!Colocar Rectangulos!!!  % fazer o try no ROI e nao no criar BKG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
