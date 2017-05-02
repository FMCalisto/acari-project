clear all

obj = mmreader('Proj1/video/experiencedmale7.avi');
nFrames = obj.NumberOfFrames;
step = 19;%19

info = get(obj);

imgbk = read(obj,1);

TouchBox = uicontrol('style','text');
%MaleDistanceBox = uicontrol('style','text');
%FemaleDistanceBox = uicontrol('style','text');
%FrameFirstCopulaBox = uicontrol('style','text');
%NumberCopulaBox = uicontrol('style','text');
%NumberFramesCopulaBox = uicontrol('style','text');

%-----------------------------
%Parametros
thr = 28; 
minArea = 55;

firstTouch = false;
touchCounter = 0;
copulaTimer = 0;
firstCopulaFrame = 0;

maleDistance = 0;
femaleDistance = 0;
maleAllPositions = [];
femaleAllPositions = [];
copulaDistance = 0;
copulaAllPositions = [];
distanceBetweenMaleAndFemale = [];


touchDetectedFrames = [];


%http://stackoverflow.com/questions/3533843/how-to-draw-a-line-on-an-image-in-matlab
%-----------------------------

se = strel('disk',3);
%-------------------------------------
Bkg=zeros(size(imgbk));
alfa=0.03;
figure; hold on

for k = 1 : step : 2642
    %disp(k);
    imgfr = read(obj,k);
    Y = imgfr;
    Bkg = alfa * double(Y) + (1 - alfa) * double(Bkg);
    imgUInt8 = uint8(Bkg);
    %imshow(imgUInt8); drawnow
end
%-------------------------------------

imgBkgBase = imgUInt8;

%We will use this auxCounter, because we want to increment the positions
%on the [male|female|copula]AllPositions array in less steps than the one
%defined on 'step'
auxCounter = 1;
auxCounterCopula = 1;
for k=1: step :nFrames
    
    imgfr = read(obj,k);
    hold off
    figure(1),imshow(imgfr,[]);
    
    imgdif = (abs(double(imgBkgBase(:,:,1))-double(imgfr(:,:,1)))>thr) | ...
    (abs(double(imgBkgBase(:,:,2))-double(imgfr(:,:,2)))>thr) | ...
    (abs(double(imgBkgBase(:,:,3))-double(imgfr(:,:,3)))>thr);
    
    bw = imclose(imgdif,se);
    %imshow(bw);
    [lb num] = bwlabel(bw);
    regionProps = regionprops(lb,'area','FilledImage','Centroid');
    inds = find([regionProps.Area]>minArea);
    regnum = length(inds);
    
    set(TouchBox,'String','Not touching')
    set(TouchBox,'Position',[40,10,100,25])
    
    %Male-Female Check
    if(regnum > 1)
        firstTouch = true;
        if ([regionProps(inds(1),1).Area] < [regionProps(inds(2),1).Area]);
            maleIndex = 1;
            maleArea = [regionProps(inds(1),1).Area];
            femaleIndex = 2;
            femaleArea = [regionProps(inds(2),1).Area];
        else
            maleIndex = 2;
            maleArea = [regionProps(inds(2),1).Area];
            femaleIndex = 1;
            femaleArea = [regionProps(inds(1),1).Area];
        end
        
        xFemale = regionProps(inds(femaleIndex),1).Centroid(1);
        yFemale = regionProps(inds(femaleIndex),1).Centroid(2);
        xMale = regionProps(inds(maleIndex),1).Centroid(1);
        yMale = regionProps(inds(maleIndex),1).Centroid(2);
        distanceBetweenMaleAndFemale = [distanceBetweenMaleAndFemale; sqrt((xFemale-xMale).^2) + ((yFemale-yMale).^2)];
        
        maleAllPositions = [maleAllPositions; [regionProps(inds(maleIndex),1).Centroid]];
        femaleAllPositions = [femaleAllPositions; [regionProps(inds(femaleIndex),1).Centroid]];
        
        if k > 1
            line(maleAllPositions(auxCounter,1), maleAllPositions(auxCounter, 2),'Marker','+', 'MarkerSize', 10, 'Color', 'b');
            line(femaleAllPositions(auxCounter,1), femaleAllPositions(auxCounter,2), 'Marker', 'o', 'MarkerSize', 10, 'Color', 'm');
        end
        auxCounter = auxCounter + 1;
    else
        %Touch Detected - number of active regions = 1
        copulaAllPositions = [copulaAllPositions; [regionProps(inds(1),1).Centroid]];
        touchDetectedFrames= [touchDetectedFrames;k];
        distanceBetweenMaleAndFemale = [distanceBetweenMaleAndFemale; 0];
        
        set(TouchBox,'String','Touch detected')
        set(TouchBox,'Position',[40,10,100,25])
        
        if(firstTouch)
            touchCounter = touchCounter + 1;
            firstTouch = false;
        end
        if(touchCounter >1)
            firstCopulaFrame = k;
            copulaTimer = copulaTimer + 1;
        end
        if k > 1
            line(copulaAllPositions(auxCounterCopula,1), copulaAllPositions(auxCounterCopula, 2),'Marker','p', 'MarkerSize', 10, 'Color', 'g');
        end
        auxCounterCopula = auxCounterCopula + 1;
    end
    
    if regnum
        for j=1:regnum
                [lin col]= find(lb == inds(j));
                upLPoint = min([lin col]);
                dWindow  = max([lin col]) - upLPoint + 1;

                rectangle('Position',[fliplr(upLPoint) fliplr(dWindow)],'EdgeColor',[1 1 0],...
                    'linewidth',2);
                %center = regionProps(j).Centroid;
                %line(center(1,1), center(1,2),'Marker','+');hold on
        end
    end
    
    for i = 1 : size(maleAllPositions,1)
        line(maleAllPositions(i,1),maleAllPositions(i,2),'Marker','+', 'MarkerSize', 5, 'Color', 'b');
    end
    for i = 1 : size(femaleAllPositions,1)
        line(femaleAllPositions(i,1),femaleAllPositions(i,2), 'Marker', 'o', 'MarkerSize', 5, 'Color', 'm');
    end    
    for i = 1 : size(copulaAllPositions,1)
        line(copulaAllPositions(i,1),copulaAllPositions(i,2), 'Marker','p', 'MarkerSize', 4, 'Color', 'g');
    end
    subplot(2,1,2), plot(distanceBetweenMaleAndFemale), title('Distance between Male and Female Tetranychus'), xlabel('Frames (x19)'), ylabel('Distance (px)'),
    subplot(2,1,1), drawnow

    
end

%EUCLIDEAN DISTANCE CALCULATION FOR FEMALE; MALE; COPULA
for k = 1: length(copulaAllPositions)-1
    x1 = copulaAllPositions(k,1);
    y1 = copulaAllPositions(k,2);
    x2 = copulaAllPositions(k + 1,1);
    y2 = copulaAllPositions(k + 1,2);
    
    dist = sqrt((x2-x1).^2) + ((y2-y1).^2);
    copulaDistance = copulaDistance + dist;
end
for k = 1: length(femaleAllPositions)-1
    x1 = femaleAllPositions(k,1);
    y1 = femaleAllPositions(k,2);
    x2 = femaleAllPositions(k + 1,1);
    y2 = femaleAllPositions(k + 1,2);
    
    dist = sqrt((x2-x1).^2) + ((y2-y1).^2);
    femaleDistance = femaleDistance + dist;
end
for k = 1: length(maleAllPositions)-1
    x1 = maleAllPositions(k,1);
    y1 = maleAllPositions(k,2);
    x2 = maleAllPositions(k + 1,1);
    y2 = maleAllPositions(k + 1,2);
    
    dist = sqrt((x2-x1).^2) + ((y2-y1).^2);
    maleDistance = maleDistance + dist;
end
maleDistance = maleDistance + copulaDistance;
femaleDistance = femaleDistance + copulaDistance;

SummaryBox1 = uicontrol('style','text');
str = {'Distancia percorrida pelo macho (px): ',num2str(maleDistance)};
set(SummaryBox1,'String',str)
set(SummaryBox1,'Position',[40,100,200,50])

SummaryBox2 = uicontrol('style','text');
str = {'Distancia percorrida pela femea (px): ',num2str(femaleDistance)};
set(SummaryBox2,'String',str)
set(SummaryBox2,'Position',[40,150,200,50])

SummaryBox3 = uicontrol('style','text');
str = {'Frames até a primeira copula : ',num2str(firstCopulaFrame)};
set(SummaryBox3,'String',str)
set(SummaryBox3,'Position',[40,200,200,50])


