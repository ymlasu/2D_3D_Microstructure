
clear all;
S2_r_w1 = textread('TS2_w1.txt');
S2_r_w1(S2_r_w1(:,2)==0,:)=[];
S2_w1(:,1) = S2_r_w1(:,1);
S2_w1(:,2) = (S2_r_w1(:,2)-S2_r_w1(1,2).^2)./(S2_r_w1(1,2).*(1-S2_r_w1(1,2)));

S2_r_w2 = textread('TS2_w2.txt');
S2_r_w2(S2_r_w2(:,2)==0,:)=[];
S2_w2(:,1) = S2_r_w2(:,1);
S2_w2(:,2) = (S2_r_w2(:,2)-S2_r_w2(1,2).^2)./(S2_r_w2(1,2).*(1-S2_r_w2(1,2)));


S2_r_w3 = textread('TS2_w3.txt');
S2_r_w3(S2_r_w3(:,2)==0,:)=[];
S2_w3(:,1) = S2_r_w3(:,1);
S2_w3(:,2) = (S2_r_w3(:,2)-S2_r_w3(1,2).^2)./(S2_r_w3(1,2).*(1-S2_r_w3(1,2)));


figure;
hold on;
plot(S2_w1(:,1),S2_w1(:,2));
plot(S2_w2(:,1),S2_w2(:,2));
plot(S2_w3(:,1),S2_w3(:,2));


S2(:,1) = S2_w3(1:120,1);
S2(:,2) =  S2_w3(1:120,2);


% S2(:,1) = S2_w1(1:120,1);
% S2(:,2) = mean([S2_w1(1:120,2),S2_w2(1:120,2),S2_w3(1:120,2)],2);
% Vomlume_Frac = mean([S2_r_w1(1,2),S2_r_w2(1,2),S2_r_w3(1,2)]);
figure;
plot(S2(:,1),S2(:,2));
