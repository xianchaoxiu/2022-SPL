close all; clear; clc;
addpath(genpath(pwd));

% -----------------------------------
% Author: ...
% Date: 12-July-2020
% -----------------------------------


%% 数据生成
u   = normrnd(0,1,[500 1]);
v11 = ones(1,10);
v12 = -ones(1,10);
v13 = -zeros(1,80);
v1  = [v11 v12 v13];
v21 = zeros(1,80);
v22 = ones(1,10);
v23 = -ones(1,10);
v2  = [v21 v22 v23];
e1  = normrnd(0,0.01,[1,100]);
e2  = normrnd(0,0.01,[1,100]);
X1 = u*(v1+e1);
X2 = u*(v2+e2);

% 验证X列稀�?
% X1 = u*v1;
% X2 = u*v2;
% 另外�?种加入噪�?
% X1 = u*v1+e1;
% X2 = u*v2+e2;

%% 投影维数
d = 10;

%% 第一种方法：GCCA
tic;
[U_GCCA,P1_GCCA,P2_GCCA] = GCCA(X1,X2,d);
Time_GCCA = toc;

%% 第二种方法：SGCCA
r1 = 1;
r2 = 10; 
tic;
[U_SGCCA,P1_SGCCA,P2_SGCCA] = SGCCA(X1,X2,r1,r2,d);
Time_SGCCA = toc;

%% 第三种方法：JSGCCA
r1 = 1;
r2 = 10; 
tic;
[U_JSGCCA,P1_JSGCCA,P2_JSGCCA] = JSGCCA(X1,X2,r1,r2,d);
Time_JSGCCA = toc;

%% 第四种发方法：JJSCGCCA
s1 = 10;
s2 = 10;
tic;
[U_JSCGCCA,P1_JSCGCCA,P2_JSCGCCA] = JSCGCCA(X1,X2,s1,s2,d);
Time_JSCGCCA = toc;


