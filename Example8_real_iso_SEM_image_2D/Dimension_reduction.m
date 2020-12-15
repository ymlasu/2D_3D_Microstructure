
eigenVal_RX = flip(eigenVal_RX);
eigenVal_RY = flip(eigenVal_RY);
eigenVec_RX= flip(eigenVec_RX,2);
eigenVec_RY= flip(eigenVec_RY,2);

eigenVal = kron(eigenVal_RX,eigenVal_RY);
[eigenVal, index] = sort(eigenVal,'descend');
GMP_quantile = CDF_Ref;
GMP_ygrid = Cooordinate_Ref + Mean_GMP;
GMP_Marginal_mean = Mean_GMP;

%%% selected number of kl terms, i.e. the number of reduced dimensionality
M_terms = 821;
eigenVal = eigenVal(1:M_terms);
index = index(1:M_terms);
index = int16(index);

% resolution
N_grid_Y = 999;
N_grid_X = 999;
tic
eigenVec = [];
i = 1:length(eigenVal_RX)*length(eigenVal_RY);
I = rem(i,N_grid_Y+1);
I(I==0) = N_grid_Y+1;
for j = 1:M_terms
JJ = index(j);
J = rem(JJ,N_grid_Y+1);
if J == 0
J = N_grid_Y+1;
else
end
eigenVec(:,j) = eigenVec_RX(idivide(int16(i),int16(N_grid_Y+1),'ceil'),idivide(JJ,int16(N_grid_Y+1),'ceil')).*eigenVec_RY(I,J);
end
toc
eigenVal_GMP = eigenVal;
GMP_quantile = CDF_Ref;
GMP_ygrid = Cooordinate_Ref + Mean_GMP;
GMP_Marginal_mean = Mean_GMP;
% % save('microstructure_data_RBTO.mat','GMP_quantile','GMP_ygrid','eigenVec','eigenVal_GMP','GMP_Marginal_Var','GMP_Marginal_mean')
% % save('microstructure_data_RBTO_small.mat','GMP_quantile','GMP_ygrid','eigenVec','eigenVal_GMP','GMP_Marginal_Var','GMP_Marginal_mean')