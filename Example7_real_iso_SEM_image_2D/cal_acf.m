clear all;close all;clc;

load('real_SEM_images.mat');
w = w1;
w(w1==1) = 0;
w(w1==0) = 1;
W1 = w;

w = w2;
w(w2==1) = 0;
w(w2==0) = 1;
W2 = w;

w = w3;
w(w3==1) = 0;
w(w3==0) = 1;
W3 = w;

figure;
image(w,'CDataMapping','scaled');
colormap('gray');
grid off;
set(gca,'visible','off');

% % n_X = min([size(W1,2),size(W2,2),size(W3,2)]);
% % W_X = [W1(:,1:n_X)' W2(:,1:n_X)' W3(:,1:n_X)'];
% % n_Y = min([size(W1,1),size(W2,1),size(W3,1)]);
% % W_Y = [W1(1:n_Y,:) W2(1:n_Y,:) W3(1:n_Y,:)];

n_X = size(W3,2);
W_X = [W3(:,1:n_X)'];
n_Y = size(W3,1);
W_Y = [W3(1:n_Y,:)];


% Vomlume_Frac = 1-mean([sum(sum(W1)),sum(sum(W2)),sum(sum(W3))])/n_Y/n_X;

[WW_SEM_ACF_Y,a] = corrPairs(W_Y);
figure
plot(0:length(WW_SEM_ACF_Y)-1,WW_SEM_ACF_Y);




% axis([0 300 -0.2 1 ]);
% xlabel('Vertical direction/Pixel');
% ylabel('Autocorrelation function');

[WW_SEM_ACF_X,a] = corrPairs(W_X);
figure
plot(0:length(WW_SEM_ACF_X)-1,WW_SEM_ACF_X);

