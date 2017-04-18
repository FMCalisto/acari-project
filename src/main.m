

clear all

vid = VideoReader('../SonMated/SonofMated2.avi');

nFrame= 40*25;
step=20;

img=read(vid, 1);
Bkg=zeros(size(img));
alfa=0.05;
figure; hold on

%Exprimentar varios valores para ALPHA

for i = 1 : step : nFrame
    i
    img = read(vid, i);
    Y = img;
    Bkg = alfa * double(Y) + (1-alfa) * double(Bkg);
    imshow(uint8(Bkg)); drawnow
end

% ------------ Fim ------------------- %

thr = 40;
minArea = 200;
baseNum = 1350;

vid4D = zeros([vid.Height vid.Width 3 nFrames/step]);
figure,
k = 1;

se = strel('disk',3);

for i = 1 : step : nFrames
    i
    img = read(vid, i);
    vid4D(:,:,:,k) = img;
    %imshow(img); drawnow
    k = k + 1;
    %pause
    
    Y = img;
    Bkg = alfa * double(Y) + (1-alfa) * double(Bkg);
    
    imshow(uint8(Bkg)); drawnow
    
end

bkg = median(vid4D, 4);
figure,
%imshow(uint8(bkg));

vid = VideoReader('../SonMated/SonofMated2.avi');
nFrames = 40 * 25;
step = 20;

img=read(vid,1);
Bkg=zeros (size(img));

figure; hold on
alpha=0.05; %Exprimentar varios valores para ALPHA

bkg = median(vid4D, 4);
figure,
%imshow(uint8(bkg));