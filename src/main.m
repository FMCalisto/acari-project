

clear all

imgbk = imread('Pedestrian_dataBase/ped7c0000.tif');

thr = 40;
minArea = 200;

baseNum = 1350;
seqLength = 100;

% baseNum = 1374;
% seqLength = 0;
% 
%imshow(imgdif)
se = strel('disk',3);

figure;
for i=0:seqLength
    imgfr = imread(sprintf('Pedestrian_dataBase/ped7c%.4d.tif',baseNum+i));
    hold off
    imshow(imgfr);
    
    imgdif = (abs(double(imgbk(:,:,1))-double(imgfr(:,:,1)))>thr) | ...
        (abs(double(imgbk(:,:,2))-double(imgfr(:,:,2)))>thr) | ...
        (abs(double(imgbk(:,:,3))-double(imgfr(:,:,3)))>thr);
    
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
           
            rectangle('Position',[fliplr(upLPoint) fliplr(dWindow)],'EdgeColor',[1 1 0],...
                'linewidth',2);
        end
    end
    drawnow
end
    
    
