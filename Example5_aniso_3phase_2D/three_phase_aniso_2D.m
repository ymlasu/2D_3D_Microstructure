
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
%  random 3-phase materials with anisotropic ACF.

clear all;clc;close all;


sigma(1,1,1) = 10;
sigma(1,1,2) = 10;

sympref('HeavisideAtOrigin', 0);
II = [0 0.5 1];   %%% phase indictor representation
Vomlume_Frac = [0.3 0.4 0.3];
Mean_Marginal = sum(II.*Vomlume_Frac);
Var_Marginal = sum(((II-Mean_Marginal).^2).*Vomlume_Frac);

%%
%********GMP Margin ( Zero mean, sampling variance)********
% mu1 = [212;230;250];
a_rand = lhsdesign(1e5,1);
loc1 = find(a_rand <= Vomlume_Frac(1));
loc2 = find(a_rand > Vomlume_Frac(1) & a_rand <= (Vomlume_Frac(1)+Vomlume_Frac(2)));
loc3 = find(a_rand > (Vomlume_Frac(1)+Vomlume_Frac(2)));
aa_rand(loc1,1) = normrnd(212,sqrt(10),length(loc1),1);
aa_rand(loc2,1) = random('Lognormal', 5.4380,0.013748,length(loc2),1);
aa_rand(loc3,1) = unifrnd(240, 250, length(loc3),1);
Mean_GMP = mean(aa_rand);
GMP_Marginal_Var = var(aa_rand);
[CDF_Ref, Cooordinate_Ref] = ecdf(aa_rand);
% % % figure
% % % plot(Cooordinate_Ref,CDF_Ref);
Cooordinate_Ref = Cooordinate_Ref- Mean_GMP;
GMP_Marginal_icdf = @(y)interp1(CDF_Ref,Cooordinate_Ref,y,'pchip');





N_grid = 159;

[R_Y,R_Y_target] = acf_updating_x(N_grid,Vomlume_Frac);
[eigenVec_RY,eigenVal_RY]= eig(R_Y);
eigenVal_RY = diag(eigenVal_RY);

N_grid_y = 299;
[R_X,R_X_target] = acf_updating_y(N_grid_y,Vomlume_Frac);
[eigenVec_RX,eigenVal_RX]= eig(R_X);
eigenVal_RX = diag(eigenVal_RX);


%%
%%%%large number of generations

N_length = N_grid;

tic
for j= 1:1e3
eta = normrnd(0,1,N_grid+1,size(eigenVec_RX,2));
ksi_Inter = eta.*(sqrt(eigenVal_RX))'*eigenVec_RX';

W = eigenVec_RY*(ksi_Inter.*(repmat(sqrt(eigenVal_RY),1,N_grid+1)));


WW_Y(:,j) =  W(1:N_length+1,1);
WW_Y1(:,j) = W(1:N_length+1,10);
WW_X(:,j) =  W(1,1:N_length+1)';

end
toc







GMP_W_Y = GMP_Marginal_icdf(normcdf(WW_Y));
GMP_W_Y = GMP_W_Y+Mean_GMP;
GMP_W_Discrete_Y = zeros(size(GMP_W_Y));
GMP_W_Discrete_Y(GMP_W_Y<=(GMP_Marginal_icdf(Vomlume_Frac(1)) + Mean_GMP)) = II(1);
GMP_W_Discrete_Y(GMP_W_Y>(GMP_Marginal_icdf(Vomlume_Frac(1)) + Mean_GMP) & GMP_W_Y<=(GMP_Marginal_icdf(Vomlume_Frac(1)+Vomlume_Frac(2)) + Mean_GMP)) = II(2);
GMP_W_Discrete_Y(GMP_W_Y>(GMP_Marginal_icdf(Vomlume_Frac(1)+Vomlume_Frac(2)) + Mean_GMP)) = II(3);


GMP_W_X = GMP_Marginal_icdf(normcdf(WW_X));
GMP_W_X = GMP_W_X+Mean_GMP;
GMP_W_Discrete_X = zeros(size(GMP_W_X));
GMP_W_Discrete_X(GMP_W_X<=(GMP_Marginal_icdf(Vomlume_Frac(1)) + Mean_GMP)) = II(1);
GMP_W_Discrete_X(GMP_W_X>(GMP_Marginal_icdf(Vomlume_Frac(1)) + Mean_GMP) & GMP_W_X<=(GMP_Marginal_icdf(Vomlume_Frac(1)+Vomlume_Frac(2)) + Mean_GMP)) = II(2);
GMP_W_Discrete_X(GMP_W_X>(GMP_Marginal_icdf(Vomlume_Frac(1)+Vomlume_Frac(2)) + Mean_GMP)) = II(3);



tic
eta = normrnd(0,1,size(eigenVec_RY,2),size(eigenVec_RX,2));
ksi_Inter = eta.*(sqrt(eigenVal_RX))'*eigenVec_RX';
W = eigenVec_RY*(ksi_Inter.*(repmat(sqrt(eigenVal_RY),1,size(eigenVec_RX,2))));
GMP_W = GMP_Marginal_icdf(normcdf(W))+Mean_GMP;
GMP_W_Discrete = zeros(size(GMP_W));
GMP_W_Discrete(GMP_W<=(GMP_Marginal_icdf(Vomlume_Frac(1)) + Mean_GMP)) = II(1);
GMP_W_Discrete(GMP_W>(GMP_Marginal_icdf(Vomlume_Frac(1)) + Mean_GMP) & GMP_W<=(GMP_Marginal_icdf(Vomlume_Frac(1)+Vomlume_Frac(2)) + Mean_GMP)) = II(2);
GMP_W_Discrete(GMP_W>(GMP_Marginal_icdf(Vomlume_Frac(1)+Vomlume_Frac(2)) + Mean_GMP)) = II(3);
toc


figure;
image(GMP_W_Discrete,'CDataMapping','scaled');
% colormap([gray;winter;hsv]);
% colormap([gray;winter;autumn]);
colormap(gray);
grid off;
set(gca,'visible','off');


A=autocorr2d(GMP_W_Discrete);
figure
plot(A(501:end,501));
acf_y_reference = A(501:end,501);
save('acf_y_reference.mat','acf_y_reference');
figure
plot(A(501,501:end));
acf_x_reference = A(501,501:end)';
save('acf_x_reference.mat','acf_x_reference');

figure;
image(GMP_W,'CDataMapping','scaled');
% colormap('gray');
grid off;
set(gca,'visible','off');
colorbar;






%%

% figure;
% hist(WW_Y(1,:),100);
% mu = mean(WW_Y,2);
% variance = var(WW_Y,0,2);
% 
% 
[WW_samples_ACF_Y,a] = corrPairs(GMP_W_Discrete_Y);
figure
plot(0:N_length,WW_samples_ACF_Y);
hold on
plot(0:N_length,R_Y_target(1,1:N_length+1),'-.');
axis([0 300 -0.2 1 ]);
xlabel('Along vertical direction/Pixel');
ylabel('Autocorrelation function');
% 
[WW_samples_ACF_X,~] = corrPairs(GMP_W_Discrete_X);
figure
plot(0:N_length,WW_samples_ACF_X);
hold on
plot(0:N_length,R_X_target(1,1:N_length+1),'-.');
axis([0 300 -0.2 1 ]);
xlabel('Along horizontal direction/Pixel');
ylabel('Autocorrelation function');



[CDF, Cooordinate] = ecdf(GMP_W_Y(1,:));
figure
plot(Cooordinate,CDF);
hold on
plot(Cooordinate_Ref+Mean_GMP,CDF_Ref);
% plot(GMP_Marginal_icdf(CDF)+Vomlume_Frac*(212-230)+230,CDF,'-.');



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



