function f = plot_3D(GMP_W,GMP_W_Binary,N_grid)


GMP_W_1(:,:)=GMP_W(:,:,1);
GMP_W_2(:,:)=GMP_W(:,end,:);
GMP_W_3(:,:)=GMP_W(1,:,:);
[X Y] = meshgrid(0:N_grid,0:N_grid);
figure
h_GMP = surf(zeros(size(X)),X,Y,flip(flip(GMP_W_1,1),2));
xlabel('x(Pixel)');
ylabel('y(Pixel)');
zlabel('z(Pixel)');
set(h_GMP,'edgecolor','none');
hold on;
h_GMP = surf(X,zeros(size(X))*N_grid,Y,flip(GMP_W_2,1));
set(h_GMP,'edgecolor','none');
h_GMP = surf(X,Y,ones(size(X))*N_grid,flip(GMP_W_3,1));
set(h_GMP,'edgecolor','none');
plot3(zeros(N_grid+1,1),zeros(N_grid+1,1),0:N_grid,'k:');
plot3(zeros(N_grid+1,1),0:N_grid,ones(N_grid+1,1)*N_grid,'k:');
plot3(0:N_grid,zeros(N_grid+1,1),ones(N_grid+1,1)*N_grid,'k:');
% set(gca,'visible','off');
axis([0 N_grid+1 0 N_grid+1 0 N_grid+1]);
AZ = -43.9805;
EL = 27.7508;
view(AZ,EL);
colorbar;
% colormap(jet);
% colormap(autumn);
% colormap([bone;spring;jet;winter;]);
% colormap([winter;cool;parula;jet]);











%%
%%% x=0 in plot corresponding toflip(flip(W(:,:,1),1),2) 
Binary_W1(:,:)=GMP_W_Binary(:,:,1);
%%
%%% y=0 in plot corresponding to flip(W(:,end,:),1)
Binary_W2(:,:)=GMP_W_Binary(:,end,:);
%%
%%% Z= z_top in plot corresponding to flip(W(1,:,:),1)
Binary_W3(:,:)=GMP_W_Binary(1,:,:);
[X Y] = meshgrid(0:N_grid,0:N_grid);
figure
h = surf(zeros(size(X)),X,Y,flip(flip(Binary_W1,1),2));
xlabel('x(Pixel)');
ylabel('y(Pixel)');
zlabel('z(Pixel)');
set(h,'edgecolor','none');
hold on;
h = surf(X,zeros(size(X))*N_grid,Y,flip(Binary_W2,1));
set(h,'edgecolor','none');
h = surf(X,Y,ones(size(X))*N_grid,flip(Binary_W3,1));
set(h,'edgecolor','none');
plot3(zeros(N_grid+1,1),zeros(N_grid+1,1),0:N_grid,'k:');
plot3(zeros(N_grid+1,1),0:N_grid,ones(N_grid+1,1)*N_grid,'k:');
plot3(0:N_grid,zeros(N_grid+1,1),ones(N_grid+1,1)*N_grid,'k:');
% set(gca,'visible','off');
axis([0 N_grid+1 0 N_grid+1 0 N_grid+1]);
colormap([gray;autumn;summer]);
AZ = -43.9805;
EL = 27.7508;
view(AZ,EL);
end




