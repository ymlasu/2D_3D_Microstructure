
% The codes are freely distributed as complementary material to the article:
% Y. Gao, Y. Jiao, Y. Liu, Ultra-efficient reconstruction of 3D microstructure 
% and distribution of properties of random heterogeneous materials containing 
% multiple phases, Acta Mater. 204 (2021) 116526. 
% https://doi.org/https://doi.org/10.1016/j.actamat.2020.116526.
%--------------------------------------------------------------------------
% Copyright Notice
%    Authors:   Yi Gao             <ygao115@asu.edu>
%               Yongming Liu       <yongming.liu@asu.edu>    
%               Prognostic Analysis and Reliability Assessment Lab,
%                Arizona State University
%               https://paralab.engineering.asu.edu/
%--------------------------------------------------------------------------

% Notes:
%  2D microstructure and material property fields reconstruction for
%  random biphase materials based on the real SEM image.


clear all;clc;close all;
N_grid = 999;
x = 0:N_grid;
Vomlume_Frac =0.3205;
sigma(1,1,1) = 10;
sigma(1,1,2) = 10;
load('acf_w4_sample_autocorrelation.mat');


%%
%********GMP Margin ( Zero mean, sampling variance)********
mu = [212;230]-Vomlume_Frac*(212-230)-230;
GMP_Marginal = gmdistribution(mu,sigma,[Vomlume_Frac,1-Vomlume_Frac]);
Sam_GMP_Marginal = random(GMP_Marginal,1e5);
GMP_Marginal_Var = var(Sam_GMP_Marginal);
% % % % % % GM_sigma = Vomlume_Frac*(mu(1).^2+sigma(1,1,1).^2-GM_mu)+(1-Vomlume_Frac)*(mu(2).^2+sigma(1,1,2).^2-GM_mu);

GMP_ygrid = (-25:0.1:25)';
GMP_quantile = cdf(GMP_Marginal,GMP_ygrid);
GMP_Marginal_icdf = @(y)interp1(GMP_quantile,GMP_ygrid,y,'pchip');

%%
%********Binary Margin (Not zero mean, sampling variance)********
Binary_mu = [212;230];
% % % % % mean_Trans = Vomlume_Frac*(212-230)+230;
Threshold_binary = Vomlume_Frac*(212-230)+230+GMP_Marginal_icdf(Vomlume_Frac); %%%icdf(Vomlume_Frac);
Binary_Marginal = gmdistribution(Binary_mu,sigma,[Vomlume_Frac,1-Vomlume_Frac]);
Sam_Marginal = random(Binary_Marginal,1e5);
Binary_Marginal_Var = var(Sam_Marginal);

Binary_ygrid = (190:0.1:250)';
Binary_quantile = cdf(Binary_Marginal,Binary_ygrid);
Binary_Marginal_icdf = @(y) (interp1(Binary_quantile,Binary_ygrid,y,'pchip')>=Threshold_binary);


%%
%%%%% ACF function 
N_grid_Y = length(acf_y)*2;
% N_grid_Y = 247;
% N_grid_Y = 120;
[R_Y,R_Y_target] = acf_updating_y(N_grid_Y,Vomlume_Frac);
[eigenVec_RY,eigenVal_RY]= eig(R_Y);
eigenVal_RY = diag(eigenVal_RY);
eigenVal_RY(eigenVal_RY<0) = 0;

N_grid_X = length(acf_x)*2;
% N_grid_X = 374;
% N_grid_X = 130;
[R_X,R_X_target] = acf_updating_x(N_grid_X,Vomlume_Frac);
[eigenVec_RX,eigenVal_RX]= eig(R_X);
eigenVal_RX = diag(eigenVal_RX);
eigenVal_RX(eigenVal_RX<0) = 0;

%%
%%%%% S2 function 
% N_grid_Y = length(WW_SEM_ACF_Y)-1;
% [R_Y,R_Y_target] = acf_updating_S2(N_grid_Y,Vomlume_Frac);
% [eigenVec_RY,eigenVal_RY]= eig(R_Y);
% eigenVal_RY = diag(eigenVal_RY);
% 
% N_grid_X = length(WW_SEM_ACF_X)-1;
% [R_X,R_X_target] = acf_updating_S2(N_grid_X,Vomlume_Frac);
% [eigenVec_RX,eigenVal_RX]= eig(R_X);
% eigenVal_RX = diag(eigenVal_RX);


%%

tic
eta = normrnd(0,1,N_grid_Y+1,size(eigenVec_RX,2));
ksi_Inter = eta.*(sqrt(eigenVal_RX))'*eigenVec_RX';
W = eigenVec_RY*(ksi_Inter.*(repmat(sqrt(eigenVal_RY),1,N_grid_X+1)));
GMP_W = GMP_Marginal_icdf(normcdf(W)) + Vomlume_Frac*(212-230)+230;
GMP_W_Binary = zeros(size(GMP_W));
GMP_W_Binary(GMP_W>=(GMP_Marginal_icdf(Vomlume_Frac)+Vomlume_Frac*(212-230)+230)) = 1;
toc
% 

GMP_W_Binary1 = GMP_W_Binary;
GMP_W_Binary1(GMP_W_Binary==1)=0;
GMP_W_Binary1(GMP_W_Binary==0)=1;


figure(1);
image(GMP_W_Binary1,'CDataMapping','scaled');
colormap('gray');
grid off;
set(gca,'visible','off');
axis equal;
axis([1 N_grid_X 1 N_grid_Y ]);


figure;
image(GMP_W,'CDataMapping','scaled');
% colormap('gray');
grid off;
set(gca,'visible','off');
colorbar;
axis equal;
axis([1 N_grid_X 1 N_grid_Y ]);

tic
for j= 1:1e4
eta = normrnd(0,1,N_grid_Y+1,size(eigenVec_RX,2));
ksi_Inter = eta.*(sqrt(eigenVal_RX))'*eigenVec_RX';
W = eigenVec_RY*(ksi_Inter.*(repmat(sqrt(eigenVal_RY),1,N_grid_X+1)));

WW_Y(:,j) =  W(:,1);
WW_X(:,j) =  W(2,:)';
end
toc


%%

GMP_W_Y = GMP_Marginal_icdf(normcdf(WW_Y));
GMP_W_Y = GMP_W_Y+Vomlume_Frac*(212-230)+230;
GMP_W_Binary_Y = zeros(size(GMP_W_Y));
GMP_W_Binary_Y(GMP_W_Y>=(GMP_Marginal_icdf(Vomlume_Frac)+Vomlume_Frac*(212-230)+230)) = 1;


GMP_W_X = GMP_Marginal_icdf(normcdf(WW_X));
GMP_W_X = GMP_W_X+Vomlume_Frac*(212-230)+230;
GMP_W_Binary_X = zeros(size(GMP_W_X));
GMP_W_Binary_X(GMP_W_X>=(GMP_Marginal_icdf(Vomlume_Frac)+Vomlume_Frac*(212-230)+230)) = 1;





%%

% figure;
% hist(WW_Y(1,:),100);
% mu = mean(WW_Y,2);
% variance = var(WW_Y,0,2);
% 
% 
[WW_samples_ACF_Y,a] = corrPairs(GMP_W_Binary_Y);
figure
plot(0:N_grid_Y,WW_samples_ACF_Y);
hold on
plot(0:N_grid_Y,R_Y_target(1,1:N_grid_Y+1),'-.');
axis([0 300 -0.2 1 ]);
xlabel('Autocorrelation along vertical direction/Pixel');
ylabel('Autocorrelation function');
legend("Samples' ACF", "Reference ACF");

[WW_samples_ACF_X,~] = corrPairs(GMP_W_Binary_X);
figure
plot(0:N_grid_X,WW_samples_ACF_X);
hold on
plot(0:N_grid_X,R_X_target(1,1:N_grid_X+1),'-.');
axis([0 300 -0.2 1 ]);
xlabel('Autocorrelation along horizontal direction/Pixel');
ylabel('Autocorrelation function');
legend("Samples' ACF", "Reference ACF");


[CDF, Cooordinate] = ecdf(GMP_W_Y(1,:));
figure
plot(Cooordinate,CDF);
hold on
plot(GMP_Marginal_icdf(CDF)+Vomlume_Frac*(212-230)+230,CDF,'-.');
ylabel('Probability');
xlabel("Young's modulus/GPa");
legend("Samples' marginal CDF", "Reference marginal CDF");



load('real_SEM_images_example2.mat', 'w4');
tic
volume_fraction_fit = sum(sum(w4==1))/(size(w4,1)*size(w4,2));
GMModel = fitgmdist(GMP_W_Y(1,:)',2);
toc

tic
A=autocorr2d(w4);
toc
