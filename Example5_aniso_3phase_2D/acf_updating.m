function [R_Gaussian_nearestCorr,R_NonGaussianT] = acf_updating(N_grid)

x = 0:N_grid;

alpha = 1;
dis = alpha*abs(x-x'); 


% % Cov_func = @(x) exp(-(2e-3)*abs(x).^2); % covariance function

S2_r = textread('S2_216.txt');
S2(:,1) = S2_r(:,1);
S2(:,2) = (S2_r(:,2)-S2_r(1,2).^2)./(S2_r(1,2).*(1-S2_r(1,2)));
nn = length(S2);
S2(nn+1:N_grid+100,2) = S2(75,2);
S2(nn+1:N_grid+100,1) = nn:(N_grid+99);
% 
% 
S2(:,1) = S2(:,1)*2.5;
Cov_func = @(x) interp1(S2(:,1),S2(:,2),x,'pchip');


% % % % Cov_func = @(x1,x2) exp(-abs(x1-x2)).*cos(4*pi*(x1-x2)); % covariance function
% % % % Cov_func = @(x) exp(-abs(sqrt(sum(x.^2,2)))); % covariance function
% % % % Marginal_icdf =  @(y) exp(norminv(y)-0.7707)-0.7628;
% % % % Marginal_icdf =  @(y) exp(norminv(y)-0.7707)-0.7628;

sigma(1,1,1) = 10;
sigma(1,1,2) = 10;
% Vomlume_Frac = S2_r(1,2);
Vomlume_Frac = 0.3; %%% Anisotropic distance
% mu = [212;230]-Vomlume_Frac*(212-230)-230;




%%
%********GMP Margin ( Zero mean, sampling variance)********
% mu = [212;230]-Vomlume_Frac*(212-230)-230;
% GMP_Marginal = gmdistribution(mu,sigma,[Vomlume_Frac,1-Vomlume_Frac]);
% Sam_GMP_Marginal = random(GMP_Marginal,1e5);
% GMP_Marginal_Var = var(Sam_GMP_Marginal);
% % % % % % % GM_sigma = Vomlume_Frac*(mu(1).^2+sigma(1,1,1).^2-GM_mu)+(1-Vomlume_Frac)*(mu(2).^2+sigma(1,1,2).^2-GM_mu);
% 
% GMP_ygrid = (-25:0.1:25)';
% GMP_quantile = cdf(GMP_Marginal,GMP_ygrid);
% GMP_Marginal_icdf = @(y)interp1(GMP_quantile,GMP_ygrid,y,'pchip');
% 
% %%
% %********Binary Margin (Not zero mean, sampling variance)********
% Binary_mu = [212;230];
% % % % % % mean_Trans = Vomlume_Frac*(212-230)+230;
% Threshold_binary = Vomlume_Frac*(212-230)+230+GMP_Marginal_icdf(Vomlume_Frac); %%%icdf(Vomlume_Frac);
% Binary_Marginal = gmdistribution(Binary_mu,sigma,[Vomlume_Frac,1-Vomlume_Frac]);
% Sam_Marginal = random(Binary_Marginal,1e5);
% Binary_Marginal_Var = var(Sam_Marginal);
% 
% Binary_ygrid = (190:0.1:250)';
% Binary_quantile = cdf(Binary_Marginal,Binary_ygrid);
% Binary_Marginal_icdf = @(y) (interp1(Binary_quantile,Binary_ygrid,y,'pchip')>=Threshold_binary);





%%
%********discrete phase indicator random field ( Zero mean, sampling variance)********
sympref('HeavisideAtOrigin', 0);
II = [0 0.5 1];   %%% phase indictor representation
Vomlume_Frac = [0.2 0.6 0.2];
Mean_Marginal = sum(II.*Vomlume_Frac);
Var_Marginal = sum(((II-Mean_Marginal).^2).*Vomlume_Frac);
% % % % % Translation_Func = @(y) II(1)+(II(2)-II(1))*heaviside(y-norminv(Vomlume_Frac(1)))+(II(3)-II(2))*heaviside(y-norminv(Vomlume_Frac(1)+Vomlume_Frac(2)));

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



%%
R_NonGaussianT = Cov_func(dis(1,:));
% R_NonGaussianT = sparse(length(x),length(x));
% for i = 1:length(x)
%     R_NonGaussianT(i,i:length(x)) = R_NonGaussianT(1,1:length(i:length(x)));
% end
% R_NonGaussianT = R_NonGaussianT+R_NonGaussianT'-diag(ones(length(x),1)).*R_NonGaussianT(1,1);



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


%% 2D numerical trapz integral
tic
for i = 1:length(x)
    %     g = @(x,y) Marginal_icdf(normcdf(x)).*Marginal_icdf(normcdf(y)).*mvnpdf([x y],[0 0],[1 R_Gaussian(1,i);R_Gaussian(1,i) 1]);
    if i == 1
        R_Gaussian_opti(1,i) = 1;
    else
        gg=@(x,y,R) 1/Var_Marginal*(Translation_Func(x)-Mean_Marginal).*(Translation_Func(y)-Mean_Marginal)./(2*pi*sqrt(1-R.^2)).*exp(-(x.^2+y.^2-2*R.*x.*y)./(2*(1-R.^2)));            %             R_NonGaussian1(1,i) = integral2(g,-Integral_limit,Integral_limit,-Integral_limit,Integral_limit);
        %
        %             gg = @(x,y) Marginal_icdf(normcdf(x)).*Marginal_icdf(normcdf(y))./(2*pi*sqrt(1-R_Gaussian(1,i).^2)).*exp(-(x.^2+y.^2-2*R_Gaussian(1,i).*x.*y)./(2*(1-R_Gaussian(1,i).^2)));
        optimal_R_NonGaussian1 = @(R) abs(trapz(trapz_grid,trapz(trapz_grid,gg(Trapz_grid_x,Trapz_grid_y,R),2)) - R_NonGaussianT(1,i));
        [solu,fval,exitflag] = fmincon(@(R) optimal_R_NonGaussian1(R),R_Gaussian(1,i),[],[],[],[],0,1,[],options);
        R_Gaussian_opti(1,i) = solu;
    end
end
toc

R_Gaussian = R_Gaussian_opti;

for i = 1:length(x)
    if i == 1
        R_NonGaussian1(1,i) = 1;
    else
        gg=@(x,y,R) 1/Var_Marginal*(Translation_Func(x)-Mean_Marginal).*(Translation_Func(y)-Mean_Marginal)./(2*pi*sqrt(1-R.^2)).*exp(-(x.^2+y.^2-2*R.*x.*y)./(2*(1-R.^2)));            %             R_NonGaussian1(1,i) = integral2(g,-Integral_limit,Integral_limit,-Integral_limit,Integral_limit);
        optimal_R_NonGaussian1 = @(R) trapz(trapz_grid,trapz(trapz_grid,gg(Trapz_grid_x,Trapz_grid_y,R),2));
        R_NonGaussian1(1,i) = optimal_R_NonGaussian1(R_Gaussian(1,i));
    end
end



R_NonGaussian_after_iter = R_NonGaussian1;
R_Gaussian_after_iter = R_Gaussian;


[R_NOnGaussian_unique, index_R_unique] = unique(dis(1,:));


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
    else
        gg=@(x,y,R) 1/Var_Marginal*(Translation_Func(x)-Mean_Marginal).*(Translation_Func(y)-Mean_Marginal)./(2*pi*sqrt(1-R.^2)).*exp(-(x.^2+y.^2-2*R.*x.*y)./(2*(1-R.^2)));            %             R_NonGaussian1(1,i) = integral2(g,-Integral_limit,Integral_limit,-Integral_limit,Integral_limit);
        optimal_R_NonGaussian1 = @(R) trapz(trapz_grid,trapz(trapz_grid,gg(Trapz_grid_x,Trapz_grid_y,R),2));
        R_NonGaussian1(1,i) = optimal_R_NonGaussian1(R_Gaussian(1,i));
    end
end


R_NonGaussian1 = R_NonGaussian1/R_NonGaussian1(1,1);

figure
plot(0:N_grid,R_NonGaussian1(1:N_grid+1));
hold on;
plot(0:N_grid,R_NonGaussianT(1:N_grid+1),'-.');
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