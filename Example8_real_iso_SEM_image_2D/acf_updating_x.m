function [R_Gaussian_nearestCorr,R_NonGaussianT] = acf_updating_x(N_grid,Vomlume_Frac)


% load('acf_W3.mat','WW_SEM_ACF_X');
load('acf_w4_sample_autocorrelation.mat');
acf_x = acf_x/acf_x(1);
x = 0:N_grid;

alpha = 1;
dis = alpha*abs(x-x'); 

acf_x(length(acf_x)+1:N_grid+1) = acf_x(end);



sigma(1,1,1) = 10;
sigma(1,1,2) = 10;




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

R_NonGaussianT = acf_x';

Integral_limit_low = -6;
Integral_limit_upper = 6;
R_Gaussian = R_NonGaussianT;
% % R_NonGaussianT = R_NonGaussianT(1:N_grid+1);
% % R_Gaussian = R_NonGaussianT;
% warning off;
trapz_grid = linspace(Integral_limit_low,Integral_limit_upper,128);
[Trapz_grid_x,Trapz_grid_y] = meshgrid(trapz_grid,trapz_grid,128);
% parpool(6);
options = optimoptions('fmincon','Display','off','Algorithm','interior-point');

tic
for i = 1:length(x)
    if i == 1
        %             R_Gaussian_est = 0.9999;
        %             optimal_R_NonGaussian1 = @(R) abs(-1+1/Vomlume_Frac*mvncdf([norminv(Vomlume_Frac),norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian_est; R_Gaussian_est 1])+1/(1-Vomlume_Frac)*mvncdf([-norminv(Vomlume_Frac),-norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian_est; R_Gaussian_est 1]) - R_NonGaussianT(1,i));
        %             [solu,fval,exitflag] = fmincon(@(R) optimal_R_NonGaussian1(R),R_Gaussian(1,i),[],[],[],[],0,1,[],options);
        R_Gaussian_opti(1,i) = 1;
    else
        optimal_R_NonGaussian1 = @(R) abs(-1+1/Vomlume_Frac*mvncdf([norminv(Vomlume_Frac),norminv(Vomlume_Frac)],[0, 0],[1 R; R 1])+1/(1-Vomlume_Frac)*mvncdf([-norminv(Vomlume_Frac),-norminv(Vomlume_Frac)],[0, 0],[1 R; R 1]) - R_NonGaussianT(1,i));
        [solu,fval,exitflag] = fmincon(@(R) optimal_R_NonGaussian1(R),R_Gaussian(1,i),[],[],[],[],0,1,[],options);
        R_Gaussian_opti(1,i) = solu;
    end
end
toc


R_Gaussian = R_Gaussian_opti;
for i = 1:length(x)
    if i == 1
        R_NonGaussian1(1,i) = 1;
        %             R_NonGaussian1(1,i) = -1+1/Vomlume_Frac*mvncdf([norminv(Vomlume_Frac),norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian_est; R_Gaussian_est 1])+1/(1-Vomlume_Frac)*mvncdf([-norminv(Vomlume_Frac),-norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian_est; R_Gaussian_est 1]);
    else
        R_NonGaussian1(1,i) = -1+1/Vomlume_Frac*mvncdf([norminv(Vomlume_Frac),norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian(1,i); R_Gaussian(1,i) 1])+1/(1-Vomlume_Frac)*mvncdf([-norminv(Vomlume_Frac),-norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian(1,i); R_Gaussian(1,i) 1]);
    end
end




R_NonGaussian_after_iter = R_NonGaussian1;
R_Gaussian_after_iter = R_Gaussian;


[R_NOnGaussian_unique, index_R_unique] = unique(dis(1,:));


% %     R_interp = @(x) interp1(R_NOnGaussian_unique,R_NonGaussian1(1,index_R_unique),x,'pchip');
% %     R_NonGaussian(1,:) = R_NonGaussian1;
% %     R_NonGaussian(2:length(x),:) = R_interp(dis(2:end,:));

R_interp = @(x) interp1(R_NOnGaussian_unique,R_Gaussian(1,index_R_unique),x,'pchip');
R_Gaussian(1,:) = R_Gaussian;
R_Gaussian(2:length(x),:) = R_interp(dis(2:end,:));

R_Gaussian(:,1) =  R_Gaussian(1,:);

% figure
% plot(0:N_grid,R_NonGaussian_after_iter(1,1:N_grid+1));
% hold on;
% plot(0:N_grid,R_NonGaussianT(1,1:N_grid+1));



% [R_Gaussian_nearestCorr] = nearcorr(R_Gaussian);
[R_Gaussian_nearestCorr] = nearcorr(R_Gaussian,'Tolerance',1e-8,'MaxIterations',500);
% figure
% hold on
% plot(0:N_grid,R_Gaussian(1,1:N_grid+1));
% plot(0:N_grid,R_Gaussian_nearestCorr(1,1:N_grid+1));

%%
%********Binary Autocorrelation********
R_Gaussian = R_Gaussian_nearestCorr(1,:);
R_NonGaussian1 = [];

for i = 1:length(x)
    if i == 1
        R_NonGaussian1(1,i) = 1;
        %             R_NonGaussian1(1,i) = -1+1/Vomlume_Frac*mvncdf([norminv(Vomlume_Frac),norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian_est; R_Gaussian_est 1])+1/(1-Vomlume_Frac)*mvncdf([-norminv(Vomlume_Frac),-norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian_est; R_Gaussian_est 1]);
    else
        R_NonGaussian1(1,i) = -1+1/Vomlume_Frac*mvncdf([norminv(Vomlume_Frac),norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian(1,i); R_Gaussian(1,i) 1])+1/(1-Vomlume_Frac)*mvncdf([-norminv(Vomlume_Frac),-norminv(Vomlume_Frac)],[0, 0],[1 R_Gaussian(1,i); R_Gaussian(1,i) 1]);
    end
end

R_NonGaussian1 = R_NonGaussian1/R_NonGaussian1(1,1);

figure
plot(0:N_grid,R_NonGaussian1(1:N_grid+1));
hold on;
plot(0:N_grid,R_NonGaussianT(1:N_grid+1),'-.');
% axis([0 300 -0.2 1 ]);

% [R_NOnGaussian_unique, index_R_unique] = unique(dis(1,:));
%     R_interp = @(x) interp1(R_NOnGaussian_unique,R_NonGaussian1(1,index_R_unique),x,'pchip');
%     R_NonGaussian(1,:) = R_NonGaussian1;
%     R_NonGaussian(2:length(x),:) = R_interp(dis(2:end,:));
% k = k+1;
% Error(k,1) = 100*sqrt(sum((R_NonGaussian1(1,:)-R_NonGaussianT(1,:)).^2)./sum(R_NonGaussianT(1,:).^2));

%%
%********GMP Autocorrelation********
% R_NonGaussian_GMP = [];
% for i = 1:length(x)
%     if i == 1
%         R_Gaussian_est = 0.999;
%         %         gg = @(x,y) 1./GMP_Marginal_Var*GMP_Marginal_icdf(normcdf(x)).*GMP_Marginal_icdf(normcdf(y))./(2*pi*sqrt(1-R_Gaussian_est.^2)).*exp(-(x.^2+y.^2-2*R_Gaussian_est.*x.*y)./(2*(1-R_Gaussian_est.^2)));
%         %         R_NonGaussian_GMP(1,i) = trapz(trapz_grid,trapz(trapz_grid,gg(Trapz_grid_x,Trapz_grid_y),2));
%     else
%         gg = @(x,y) 1./GMP_Marginal_Var*GMP_Marginal_icdf(normcdf(x)).*GMP_Marginal_icdf(normcdf(y))./(2*pi*sqrt(1-R_Gaussian(1,i).^2)).*exp(-(x.^2+y.^2-2*R_Gaussian(1,i).*x.*y)./(2*(1-R_Gaussian(1,i).^2)));
%         R_NonGaussian_GMP(1,i) = trapz(trapz_grid,trapz(trapz_grid,gg(Trapz_grid_x,Trapz_grid_y),2));
%         
%     end
% end
% R_NonGaussian_GMP(1,1) = 1;
% R_NonGaussian_GMP = R_NonGaussian_GMP/R_NonGaussian_GMP(1,1);
% delete(gcp('nocreate'));

%%



% Binary_R_NonGaussian_final = R_NonGaussian1;
% Binary_R_Gaussian_after_final = R_Gaussian;
% 
% %%%%% Horizontal direction
% figure
% plot(0:N_grid,Binary_R_NonGaussian_final(1,1:N_grid+1),'b');
% hold on
% plot(0:N_grid,R_NonGaussianT(1,1:N_grid+1),'-.');
% % plot(0:N_grid,Binary_R_Gaussian_after_final(1,1:N_grid+1));
% plot(0:N_grid,R_NonGaussian_GMP(1,1:N_grid+1),'k:');
% 
% %%%%% Vertical direction
% figure
% plot(0:N_grid,Binary_R_NonGaussian_final(1,1:N_grid+1:end),'b');
% hold on
% plot(0:N_grid,R_NonGaussianT(1,1:N_grid+1:end),'-.');
% % plot(0:N_grid,Binary_R_Gaussian_after_final(1,1:N_grid+1));
% plot(0:N_grid,R_NonGaussian_GMP(1,1:N_grid+1:end),'k:');

% [eigenVec,eigenVal]=eig(R_Gaussian_nearestCorr);
% eigenVal =  diag(eigenVal);

end