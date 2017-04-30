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
%baseNum = 262; % Initial Frame: 262
baseNum = 15500; % Couple Frame: 15500
seqLength = 23353;

se = strel('disk',3);

% Timer Functions %
% t = timer;
% t.StartFcn = @(~,thisEvent)disp([thisEvent.Type ' executed '...
%     datestr(thisEvent.Data.time,'SS.FFF')]);
% t.TimerFcn = @(~,thisEvent)disp([thisEvent.Type ' executed '...
%      datestr(thisEvent.Data.time,'SS.FFF')]);
% t.StopFcn = @(~,thisEvent)disp([thisEvent.Type ' executed '...
%     datestr(thisEvent.Data.time,'SS.FFF')]);
% t.Period = 5;
% t.TasksToExecute = 3;
% t.ExecutionMode = 'fixedRate';
t = datetime('now');

% -------------------- Backgroud -------------------- %
nTotalFrames = 23353; %23354Frames
nFrameBKG= 1000; % 23354 Frames used to compute background image
step = 20;       % Faz display de step em step frames
Bkg=zeros(size(imgbk));
BkgLast=zeros(size(imgbk));
alfa=0.25;
touchFigure = figure(2);
mainFigure = figure(1); 
set(touchFigure, 'Position', [630, 170, 500, 500]);
set(mainFigure, 'Position', [100, 000, 500, 1000]);
hold on
maleTrail = [];
femaleTrail = [];


        
%Exprimentar varios valores para ALPHA

for i = 0 : step : nFrameBKG
    %sprintf('BKG %d',i)
    imgfr = imread(sprintf('../frames/SonofMated10/SonofMated10%.5d.jpg', ...
                   baseNum+i));
    Y = imgfr;
    Bkg = alfa * double(Y) + (1 - alfa) * double(Bkg);
    
    imgUInt8 = uint8(Bkg);
    imgUInt8Last = uint8(BkgLast);
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
                      baseNum+i));
    
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
                var = k;
                plot(maleTrail(var, 1),maleTrail(var, 2),'*', ...
                    'Color', 'red','LineStyle','-');
                plot(femaleTrail(var, 1),femaleTrail(var, 2),'*', ...
                    'Color', 'blue','LineStyle','-');
                disp('Male Position: ');
                disp(maleTrail(var, 1));
                D = pdist2(femaleTrail(var, 1), maleTrail(var, 2));
                disp('Pairwise Distance: ');
                disp(D);
                if (D < 100)
                    disp('TOUCH');
                    % Count how much time in touch
                    % if n time in touch = couple
                    Seconds = second(t);
                    disp('======= Seconds: ');
                    disp(Seconds);
                    if (Seconds > 300)
                        disp('-----------> COUPLE <------------');
                    end
                end
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
    
end

% --------------------------------------------------- %
%                                                     % 
%                    Touch                            % 
%                                                     % 
% --------------------------------------------------- %

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