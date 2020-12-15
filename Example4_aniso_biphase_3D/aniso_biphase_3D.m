
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
%  3D microstructure and material property fields reconstruction for
%  random biphase materials with anisotropic ACF.

clear all;clc;close all;
N_grid = 499;
x = 0:N_grid;

sigma(1,1,1) = 10;
sigma(1,1,2) = 10;
Vomlume_Frac = 0.2675; %%% Anisotropic distance

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


[R_Y,R_Y_target] = acf_updating_y(N_grid,Vomlume_Frac);
[eigenVec_RY,eigenVal_RY]= eig(R_Y);
eigenVal_RY = diag(eigenVal_RY);

[R_X,R_X_target] = acf_updating_x(N_grid,Vomlume_Frac);
[eigenVec_RX,eigenVal_RX]= eig(R_X);
eigenVal_RX = diag(eigenVal_RX);

[R_Z,R_Z_target] = acf_updating_z(N_grid,Vomlume_Frac);
[eigenVec_RZ,eigenVal_RZ]= eig(R_Z);
eigenVal_RZ = diag(eigenVal_RZ);




N_length = N_grid;

tic
eta = normrnd(0,1,N_grid+1,N_grid+1,N_grid+1);
%%% 1 and 2 dimension
for i = 1:N_grid+1
    
    ksi_Inter11(:,:,i) = eta(:,:,i).*(sqrt(eigenVal_RX))'*eigenVec_RX';
 ksi_Inter22(:,:,i) = eigenVec_RY*(ksi_Inter11(:,:,i).*sqrt(eigenVal_RY));
    
    
end

%%% 3rd dimension 
ksi_Inter22 = permute(ksi_Inter22,[3,2,1]);
for i = 1:N_grid+1
    
    ksi_Inter33(:,:,i) = eigenVec_RZ*(ksi_Inter22(:,:,i).*sqrt(eigenVal_RZ));

end
W = permute(ksi_Inter33,[3,2,1]);


GMP_W = GMP_Marginal_icdf(normcdf(W)) + Vomlume_Frac*(212-230)+230;
 GMP_W_Binary = zeros(size(GMP_W));
GMP_W_Binary(GMP_W>=(GMP_Marginal_icdf(Vomlume_Frac)+Vomlume_Frac*(212-230)+230)) = 1;
toc
plot_3D(GMP_W,GMP_W_Binary,N_grid);

%%
%%%%large number of generations
tic
for k = 1:1e4
eta = normrnd(0,1,N_grid+1,N_grid+1,N_grid+1);
%%% 1 and 2 dimension
for i = 1:N_grid+1
    
    ksi_Inter11(:,:,i) = eta(:,:,i).*(sqrt(eigenVal_RX))'*eigenVec_RX';
 ksi_Inter22(:,:,i) = eigenVec_RY*(ksi_Inter11(:,:,i).*sqrt(eigenVal_RY));
    
    
end

%%% 3rd dimension 
ksi_Inter22 = permute(ksi_Inter22,[3,2,1]);
for i = 1:N_grid+1
    
    ksi_Inter33(:,:,i) = eigenVec_RZ*(ksi_Inter22(:,:,i).*sqrt(eigenVal_RZ));

end
ksi_Inter33 = permute(ksi_Inter33,[3,2,1]);


W_Z(:,k) = ksi_Inter33(1,1,1:N_length+1);
W_X(:,k) = ksi_Inter33(1,1:N_length+1,1);
W_Y(:,k) = ksi_Inter33(1:N_length+1,1,1);


end
toc



GMP_W_X = GMP_Marginal_icdf(normcdf(W_X));
GMP_W_X = GMP_W_X+Vomlume_Frac*(212-230)+230;
GMP_W_Binary_X = zeros(size(GMP_W_X));
GMP_W_Binary_X(GMP_W_X>=(GMP_Marginal_icdf(Vomlume_Frac)+Vomlume_Frac*(212-230)+230)) = 1;


GMP_W_Y = GMP_Marginal_icdf(normcdf(W_Y));
GMP_W_Y = GMP_W_Y+Vomlume_Frac*(212-230)+230;
GMP_W_Binary_Y = zeros(size(GMP_W_Y));
GMP_W_Binary_Y(GMP_W_Y>=(GMP_Marginal_icdf(Vomlume_Frac)+Vomlume_Frac*(212-230)+230)) = 1;

GMP_W_Z = GMP_Marginal_icdf(normcdf(W_Z));
GMP_W_Z = GMP_W_Z+Vomlume_Frac*(212-230)+230;
GMP_W_Binary_Z = zeros(size(GMP_W_Z));
GMP_W_Binary_Z(GMP_W_Z>=(GMP_Marginal_icdf(Vomlume_Frac)+Vomlume_Frac*(212-230)+230)) = 1;










% 
% figure;
% image(GMP_W_Binary,'CDataMapping','scaled');
% colormap('gray');
% grid off;
% set(gca,'visible','off');
% 
% 
% figure;
% image(GMP_W,'CDataMapping','scaled');
% % colormap('gray');
% grid off;
% set(gca,'visible','off');






%%

% figure;
% hist(WW_Y(1,:),100);
% mu = mean(WW_Y,2);
% variance = var(WW_Y,0,2);
% 
% 


%%
%%%%%  x y z coordinates to x y z ACF: (x,y,z) = acf(z,x,y);
%%%%%  x y z coordinates to 3-dimensional matrix: (z,y,x) = dim(1,2,3);
%%%%% 3-dimensional matrix to x y z ACF: (1,2,3) = acf(y,x,z);

[WW_samples_ACF_X,a] = corrPairs(GMP_W_Binary_X);
figure
plot(0:N_length,WW_samples_ACF_X);
hold on
plot(0:N_length,R_X_target(1,1:N_length+1),'-.');
axis([0 300 -0.2 1 ]);
xlabel('Autocorrelation along Y direction/Pixel');
ylabel('Autocorrelation function');




[WW_samples_ACF_Y,~] = corrPairs(GMP_W_Binary_Y);
figure
plot(0:N_length,WW_samples_ACF_Y);
hold on
plot(0:N_length,R_Y_target(1,1:N_length+1),'-.');
axis([0 300 -0.2 1 ]);
xlabel('Autocorrelation along Z direction/Pixel');
ylabel('Autocorrelation function');


[WW_samples_ACF_Z,~] = corrPairs(GMP_W_Binary_Z);
figure
plot(0:N_length,WW_samples_ACF_Z);
hold on
plot(0:N_length,R_Z_target(1,1:N_length+1),'-.');
axis([0 300 -0.2 1 ]);
xlabel('Autocorrelation along X direction/Pixel');
ylabel('Autocorrelation function');


[CDF, Cooordinate] = ecdf(GMP_W_X(1,:));
figure
plot(Cooordinate,CDF);
hold on
plot(GMP_Marginal_icdf(CDF)+Vomlume_Frac*(212-230)+230,CDF,'-.');



% 

% [WW_samples_ACF_diag,~] = corrPairs(WW_diag,0);
% figure
% plot(0:N_grid,WW_samples_ACF_diag);
% hold on
% R_diag = diag(R_Gaussian);
% plot(0:N_grid,R_diag(1:N_grid+1),'-.');

% R_Gaussian_reshape = reshape(R_Gaussian(1:10,:)',[],1);
% [WW_samples_ACF_ALL,~] = corrPairs(WW);
% figure
% plot(0:length(WW_samples_ACF_ALL)-1,WW_samples_ACF_ALL);
% hold on
% plot(0:length(WW_samples_ACF_ALL)-1,R_Gaussian_reshape);



