% --- Main script for generating the dithered trui image
close all; clear all;

ninit    = 5e+4;        % Upper bound on potential dot positions (PDPs) to use
dotmax   = 5e+5;        % Upper bound on number of dots to place

% --- Carry out the node dropping
xy = node_drop([0 1 0 1],ninit,dotmax,@radius_trui); 

% --- Display the resulting dithered image
plot(xy(:,1),xy(:,2),'k.','MarkerSize',2.1); axis square

% print -deps trui_01.eps

% Display the resulting dithering

figure 
A = double(imread('trui.png','PNG')); 
subplot(1,2,1); % First, show the original image in regular gray scale
image(A/255)
title('(a) Original gray scale image')
axis off; axis image; axis tight;

subplot(1,2,2); % Show the full dithered image
plot(xy(:,1),xy(:,2),'k.','MarkerSize',2.1)
axis square
title('(b) Dithered version')
% print -deps trui_01.eps

figure
subplot(1,2,1)  % Show dithering detail surrounding the location (0.8, 0.4).
plot(xy(:,1),xy(:,2),'k.','MarkerSize',6); 
% axis([0.6,1.0,0.2,0.6]); axis square
axis([0.7,0.9,0.3,0.5]); axis square

subplot(1,2,2)  % Show even closer detail around (0.8, 0.4).
plot(xy(:,1),xy(:,2),'k.','MarkerSize',8); 
axis([0.75,0.85,0.35,0.45]); axis square
% print -deps trui_02.eps


figure
subplot(1,2,1)  % Show dithering detail in bottom right corner.
plot(xy(:,1),xy(:,2),'k.','MarkerSize',8); 
% axis([0.6,1.0,0.2,0.6]); axis square
axis([0.9,1.0,0.0,0.1]); axis square
box off
set(gca,'YAxisLocation','right')

subplot(1,2,2)  % Show even closer detail in top right corner.
plot(xy(:,1),xy(:,2),'k.','MarkerSize',8); 
axis([0.9,1.0,0.9,1.0]); axis square
box off
set(gca,'XAxisLocation','top')
set(gca,'YAxisLocation','right')
% print -deps trui_03.eps
