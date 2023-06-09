function [p,A2,O,A,E,w_max3,w_min3,K,F3,Q,TTT3 ] = thirdT( ~ )
%UNTITLED17 此处显示有关此函数的摘要
u_min = -110000.*ones(180,1); %牵引力下限值,构造相应的180维行向量 值为-110000
%先构造x1的约束
% a1 = zeros(111,1); %1-2站台间[x1,x2,v]变量的下限 x1范围
% t = [130000 0 0]';
% a2 = repmat(t,32,1);%2-3站台间[x,v]变量的下限
% t = [260000 0 0]';
% a3 = repmat(t,32,1);%3-4站台间[x,v]变量的下限
% t = [390000 0 0]';
% a4 = repmat(t,32,1);%4-5站台间[x,v]变量的下限
% t = [520000 0 0]';
% a5 = repmat(t,32,1);%5-6站台间[x,v]变量的下限
% t = [650000 0 0]';
% a6 = repmat(t,15,1);
% w_x1 = [u_min;a1;a2;a3;a4;a5;a6]; %720维的列向量
% %构造x2 的约束
% u_min2 = zeros(180,1);
ax1 = zeros(126,1); %1-2站台间[x1,x2,v]变量的下限 x1范围
t = [0 135000 0]';
ax2 = repmat(t,32,1);%2-3站台间[x,v]变量的下限
t = [0 270000 0]';
ax3 = repmat(t,32,1);%3-4站台间[x,v]变量的下限
t = [0 405000 0]';
ax4 = repmat(t,32,1);%4-5站台间[x,v]变量的下限
t = [0 540000 0]';
ax5 = repmat(t,32,1);%5-6站台间[x,v]变量的下限
t = [0 675000 0]';
ax6 = repmat(t,10,1);
w_x2 = [u_min;ax1;ax2;ax3;ax4;ax5;ax6];
w_min3 = w_x2; %+ w_x1;
%==========================================================================
u_max = 500000.*ones(180,1);  %牵引力上限值,构造相应的180维行向量 值为500000
% ax1 = zeros(15,1); %1-2站台间[x1,x2,v]变量的下限 x1范围
% t = [130000 0 0]';
% ax2 = repmat(t,32,1);%2-3站台间[x,v]变量的下限
% t = [260000 0 0]';
% ax3 = repmat(t,32,1);%3-4站台间[x,v]变量的下限
% t = [390000 0 0]';
% ax4 = repmat(t,32,1);%4-5站台间[x,v]变量的下限
% t = [520000 0 0]';
% ax5 = repmat(t,32,1);%5-6站台间[x,v]变量的下限
% t = [650000 0 0]';
% ax6 = repmat(t,32,1);
% t = [650000 0 0]';
% ax7 = repmat(t,15,1);
% w_x1 = [u_max;ax1;ax2;ax3;ax4;ax5;ax6;ax7];
% %构造x2 的约束
% u_min2 = zeros(180,1);
 %1-2站台间[x1,x2,v]变量的上限 x1范围
ax1 = [675000 0 0]';
ax1 = repmat(ax1,10,1);

t = [675000 135000 250/3]';
ax2 = repmat(t,30,1);%2-3站台间[x,v]变量的下限
ax2 = [ax2;675000;135000;0;675000;135000;0];
t = [675000 270000 250/3]';
ax3 = repmat(t,30,1);%3-4站台间[x,v]变量的下限
ax3 = [ax3;675000;270000;0;675000;270000;0];
t = [675000 405000 250/3]';
ax4 = repmat(t,30,1);%4-5站台间[x,v]变量的下限
ax4 = [ax4;675000;405000;0;675000;405000;0];
t = [675000 540000 250/3]';
ax5 = repmat(t,30,1);%5-6站台间[x,v]变量的下限
ax5 = [ax5;675000;540000;0;675000;540000;0];
t = [675000 675000 250/3]';
ax6 = repmat(t,30,1);
ax6 = [ax6;675000;675000;0;675000;675000;0];
t = [675000 675000 0]';
ax7 = repmat(t,10,1);
w_x2 = [u_max;ax1;ax2;ax3;ax4;ax5;ax6;ax7];
w_max3 = w_x2; %+ w_x1;
%==========================================================================
%构建稀疏等式矩阵中的A， 构建了360*720 维的矩阵
D = [-0.004;-0.0004/3];
z = repmat(D,160,1); %构建一个 320*1 维的矩阵  repmat 复制 D向量160次，构成320*1的矩阵
b = mat2cell(z,2*ones(1,160),1); %构建以 2*1 维矩阵的 原包组，并以此为块 构建块对角矩阵
k = blkdiag(b{:});
block = zeros(20,10);
block2 = zeros(20,10);
m = blkdiag(block,k,block2);

%构建以G为对角线的 159个块 对角矩阵
G = [0 -1 -60;0 0 -1];
z = repmat(G,159,1);
b = mat2cell(z,2*ones(1,159),3);
n = blkdiag(b{:});

C = [0 1 0;0 0 1]; % 构建文章中的以C矩阵为对角块的矩阵
zp = repmat(C,160,1);
bp = mat2cell(zp,2*ones(1,160),3);
np = blkdiag(bp{:});
%构建文章中的 A1矩阵的右半部分

z = zeros(2,477);
b = zeros(2,3);
i = zeros(318,3);
k = np + [z b;n i];

%合并成 A1矩阵
block4 = zeros(20,30);
block3 = zeros(20,30);
x = blkdiag(block3,k,block4);
A = [m x];
%得到相应的解析解 E值
e1 = [-1000;-100/3];
e2 = zeros(20,1);
e4 = zeros(20,1);
e3 = repmat(e1,30,1);
e3 = [e3;0;0;0;0];
e3 = repmat(e3,5,1);
E = [e2;e3;e4];

%构建不等式矩阵 A2
blockA2 = zeros(180,180);
Zb = [1,-1,-625/12];
Zx = repmat(Zb,160,1);
Z1 = mat2cell(Zx,ones(1,160),3);
Z2 = blkdiag(Z1{:});
blockA2z1 = zeros(10,30);
blockA2z2 = zeros(10,30);
Z = blkdiag(blockA2z1,Z2,blockA2z2);
A2 = [blockA2 Z];
A2 = -A2;
O = zeros(180,1);


% 对目标函数的系数矩阵初始化 
R = 0.00001*eye(180);
%构建 位置系数 p 和速度系数 q  
%h = [100000000 0;0 100000000];
%构建H1矩阵
blockH1 = zeros(87,87);
b = [0 0 0;0 100000 0;0 0 1000000];
z2 = zeros(6,6);
H1 = blkdiag(blockH1,b,z2);
r = repmat(H1,5,1);
r1 = mat2cell(r,96*ones(1,5),96);
H = blkdiag(r1{:});

blockH = zeros(30,30);
blockH2 = zeros(30,30);
%得到相应的Q矩阵
Q = blkdiag(R,blockH,H,blockH2);

%得到相应的k^T * k 矩阵
i = zeros(360,180);
j = [1 0 0;0 1 0];
z = repmat(j,180,1);
j = mat2cell(z,2*ones(180,1),3);
k = blkdiag(j{:});
K = [i k];
K1 = K'*K;
%得到相应的Q1矩阵
Q1 = Q+0.5*K1;
c1 = zeros(210,1);
c2 = zeros(88,1);
ze1 = zeros(6,1);
c3 = [c2;135000;0;ze1];
c4 = [c2;270000;0;ze1];
c5 = [c2;405000;0;ze1];
c6 = [c2;540000;0;ze1];
c7 = [c2;675000;0;ze1];
c2 = zeros(30,1);
c = [c1;c3;c4;c5;c6;c7;c2];
n = 2*c'*Q;
p = -n';
f1 = [1;0];
f1 = repmat(f1,180,1);
f1 = mat2cell(f1,2*ones(180,1),1);
f1 = blkdiag(f1{:});
f2 = [0;1];
f2 = repmat(f2,180,1);
f2 = mat2cell(f2,2*ones(180,1),1);
f2 = blkdiag(f2{:});
F3 = [zeros(360,180) f1 f2 zeros(360,360)];
TTT3 = c;

end

