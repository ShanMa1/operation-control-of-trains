function [p,A,E,w_max,w_min,K,F1,Q,TTT1 ] = firstT( ~ )
%UNTITLED14 此处显示有关此函数的摘要
%   此处显示详细说明
u_min = -110000.*ones(180,1);  %牵引力下限值,构造相应的180维行向量 值为-110000
a1 = zeros(64,1); %1-2站台间[x,v]变量的下限
t = [135000 0]';
a2 = repmat(t,32,1);%2-3站台间[x,v]变量的下限
t = [270000 0]';
a3 = repmat(t,32,1);%3-4站台间[x,v]变量的下限
t = [405000 0]';
a4 = repmat(t,32,1);%4-5站台间[x,v]变量的下限
t = [540000 0]';
a5 = repmat(t,32,1);%5-6站台间[x,v]变量的下限
t = [675000 0]';
a6 = repmat(t,20,1);
w_min = [u_min;a1;a2;a3;a4;a5;a6]; %540维的列向量
%==========================================================================
u_max = 500000.*ones(180,1);  %牵引力下限值,构造相应的180维行向量 值为-110000
t = [135000 250/3]';
a1 = repmat(t,30,1); %1-2站台间[x,v]变量的下限
a1 = [a1;135000;0;135000;0];
t = [270000 250/3]';
a2 = repmat(t,30,1);%3-4站台间[x,v]变量的下限
a2 = [a2;270000;0;270000;0];
t = [405000 250/3]';
a3 = repmat(t,30,1);%4-5站台间[x,v]变量的下限
a3 = [a3;405000;0;405000;0];
t = [540000 250/3]';
a4 = repmat(t,30,1);%5-6站台间[x,v]变量的下限
a4 = [a4;540000;0;540000;0];
t = [675000 250/3]';
a5 = repmat(t,30,1);
a5 = [a5;675000;0;675000;0];
t = [675000 0]';
a6 = repmat(t,20,1);
w_max = [u_max;a1;a2;a3;a4;a5;a6]; %540维的列向量
%==========================================================================
%构建稀疏等式矩阵中的A， 构建了360*540 维的矩阵
D = [-0.004;-0.0004/3];
z = repmat(D,160,1); %构建一个 320*1 维的矩阵  repmat 复制 D向量160次，构成320*1的矩阵
b = mat2cell(z,2*ones(1,160),1); %构建以 2*1 维矩阵的 原包组，并以此为块 构建块对角矩阵
k = blkdiag(b{:});
j = zeros(40,180);
i = zeros(320,20);
m = [k i;j];
%构建以G为对角线的 159个块 对角矩阵
G = [-1 -60;0 -1];
z = repmat(G,159,1);
b = mat2cell(z,2*ones(1,159),2);
n = blkdiag(b{:});
p = eye(320);
z = zeros(2,318);
b = zeros(2,2);
i = zeros(318,2);
k = p + [z b;n i];
i = zeros(40,360);
j = zeros(320,40);
x = [k j;i];
A = [m x];
%得到相应的解析解 E值
e1 = [-1000;-100/3];
e2 = zeros(40,1);
e3 = repmat(e1,30,1);
e3 = [e3;0;0;0;0];
e3 = repmat(e3,5,1);
E = [e3;e2];
% 对目标函数的系数矩阵初始化 
R = 0.00001*eye(180);
%构建 位置系数 p 和速度系数 q  
%h = [100000000 0;0 100000000];
%构建H1矩阵
z1 = zeros(58,58);
b = [100000 0;0 1000000];
z2 = zeros(4,4);
H1 = blkdiag(z1,b,z2);
r = repmat(H1,5,1);
r1 = mat2cell(r,64*ones(1,5),64);
r2 = blkdiag(r1{:});
r = zeros(40,360);
r1 = zeros(320,40);
H = [r2 r1;r];
r1 = zeros(180,360);
r = zeros(360,180);
%得到相应的Q矩阵
Q = [R r1;r H];
%得到相应的k^T * k 矩阵
i = zeros(180,180);
j = [1 0];
z = repmat(j,180,1);
j = mat2cell(z,1*ones(180,1),2);
k = blkdiag(j{:});
K = [i k];
K1 = K'*K;
%得到相应的Q1矩阵
Q1 = Q+0.5*K1;
zz1 = zeros(180,1);
zz2 = zeros(58,1);
zz3 = [zz2;135000;0;0;0;0;0];
zz4 = [zz2;270000;0;0;0;0;0];
zz5 = [zz2;405000;0;0;0;0;0];
zz6 = [zz2;540000;0;0;0;0;0];
zz7 = [zz2;675000;0;0;0;0;0];
zz2 = zeros(40,1);
c = [zz1;zz3;zz4;zz5;zz6;zz7;zz2];
n = 2*c'*Q;
p = -n';
F1 = [eye(180) zeros(180,720)];
TTT1 = c;
end

