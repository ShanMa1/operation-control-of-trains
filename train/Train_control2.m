t_start = tic;
rho = 9;
[p1,A1_1,E1,w_max1,w_min1,k1,F1,QQ,cc1] = firstT( );
[p2,A2_2,O2,A1_2,E2,w_max2,w_min2,k2,F2,QQ2,cc2] = secondT( );
[p3,A2_3,O3,A1_3,E3,w_max3,w_min3,k3,F3,QQ3,cc3] = thirdT( );
[p4,A2_4,O4,A1_4,E4,w_max4,w_min4,k4,F4,QQ4,cc4] = fourthT( );
[p5,A2_5,O5,A1_5,E5,w_max5,w_min5,k5,F5,QQ5,cc5] = fifthT( );
Q1 = QQ+(rho/2)*k1'*k1; Q2 = QQ2+(rho/2)*k2'*k2; Q3 = QQ3+(rho/2)*k3'*k3; Q4 = QQ4+(rho/2)*k4'*k4; Q5 = QQ5+(rho/2)*k5'*k5;
Lambda1 = zeros(180,1);
Lambda2 = zeros(360,1);
Lambda3 = zeros(360,1);
Lambda4 = zeros(360,1);
Lambda5 = zeros(360,1);
z = zeros(900,1);
z1 = zeros(180,1);
z2 = zeros(180,1);
z3 = zeros(180,1);
z4 = zeros(180,1);
z5 = zeros(180,1);
ADMM_Max_Iter = 100;
v = [1 0];
j = repmat(v,180,1);
v = mat2cell(j,ones(180,1),2);
v1 = blkdiag(v{:});
v2 = [0 1];
j2 = repmat(v2,180,1);
v2 = mat2cell(j2,ones(180,1),2);
v2 = blkdiag(v2{:});
g1 = zeros(180,1);
g2 = zeros(360,1);
g3 = zeros(360,1);
g4 = zeros(360,1);
g5 = zeros(360,1);
w1 = zeros(540,1); w2 = zeros(720,1); w3 = zeros(720,1); w4 = zeros(720,1); w5 = zeros(720,1);
ABSTOL = 3;
RELTOL = 1e-5;
QUIET = 0;
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','MaxIterations',500,'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8,'Display','off');
for step = 1 : ADMM_Max_Iter  %此时迭代最大值为1000

            % 五辆车，平行运算 f(w_k+1) = argmin (f(w_k)+lambda(y_k-z_k)+rho/2*norm(y_k-z_k)^2 ) 其中 精简为
            % f(w_k+1) = f(w_k)+ c*w_k + d     quadprog()二次规划算法
            % ========================================================
            f1 = p1+(Lambda1'*k1-rho*g1'*k1)';
            f2 = p2+(Lambda2'*k2-rho*g2'*k2)';
            f3 = p3+(Lambda3'*k3-rho*g3'*k3)';
            f4 = p4+(Lambda4'*k4-rho*g4'*k4)';
            f5 = p5+(Lambda5'*k5-rho*g5'*k5)';
            w1q = w1; w2q = w2; w3q = w3; w4q = w4; w5q = w5;
            [w1,fval,exitflag,output] = quadprog(2.*Q1,f1,[],[],A1_1,E1,w_min1,w_max1,[],options);
            [w2,fval2,exitflag2,output2] = quadprog(2.*Q2,f2,A2_2,O2,A1_2,E2,w_min2,w_max2,[],options);
            [w3,fval3,exitflag3,output3] = quadprog(2.*Q3,f3,A2_3,O3,A1_3,E3,w_min3,w_max3,[],options);
            [w4,fval4,exitflag4,output4] = quadprog(2.*Q4,f4,A2_4,O4,A1_4,E4,w_min4,w_max4,[],options);
            [w5,fval5,exitflag5,output5] = quadprog(2.*Q5,f5,A2_5,O5,A1_5,E5,w_min5,w_max5,[],options);
            % 更新对偶变量===============================================
            
                 c1 = k1*w1; c2 = k2*w2; c3 = k3*w3; c4 = k4*w4; c5 = k5*w5;  % 代表y向量值
               if step > 1
                 Lambda1 = Lambda1 + rho * ( c1 - g1 );  % F2*z 表示从z变量中提取相应的在w中的相应y向量值
                 Lambda2 = Lambda2 + rho * ( c2 - g2 );
                 Lambda3 = Lambda3 + rho * ( c3 - g3 );
                 Lambda4 = Lambda4 + rho * ( c4 - g4 );
                 Lambda5 = Lambda5 + rho * ( c5 - g5 );
              end
            % G(z) 
            % =======================================================
%             c1 = k1*w1; c2 = k2*w2; c3 = k3*w3; c4 = k4*w4; c5 = k5*w5;  % 代表y向量值
            tmp1 = z1; tmp2 = z2; tmp3 = z3; tmp4 = z4; tmp5 = z5;
            z1 = 0.5*(c1+v1*c2);
            z2 = 0.5*(v2*c2+v1*c3);
            z3 = 0.5*(v1*c4+v2*c3);
            z4 = 0.5*(v1*c5+v2*c4);
            z5 = v2*c5;
            z = [z1;z2;z3;z4;z5];
            % 更新对偶变量===============================================
            g1 = F1*z;  %代表全局变量中相应关于y向量的对应值
            g2 = F2*z;
            g3 = F3*z;
            g4 = F4*z;
            g5 = F5*z;
            
            l1 = Lambda1; l2 = Lambda2;l3 = Lambda3;l4 = Lambda4;l5 = Lambda5;
            
                Lambda1 = Lambda1 + rho * ( c1 - g1 );  % F2*z 表示从z变量中提取相应的在w中的相应y向量值
                Lambda2 = Lambda2 + rho * ( c2 - g2 );
                Lambda3 = Lambda3 + rho * ( c3 - g3 );
                Lambda4 = Lambda4 + rho * ( c4 - g4 );
                Lambda5 = Lambda5 + rho * ( c5 - g5 );
            
            %===========================================================                          
            
            % stoppin criterion for ADMM
            % primal误差
            history.fval1(step) = fval; %+Lambda1'*(-g1)+cc1'*QQ*cc1+0.25*g1'*g1;%(w1-cc1)'*QQ*(w1-cc1);
            history.fval2(step) = fval2;%+Lambda2'*(-g2)+cc2'*QQ2*cc2+0.25*g2'*g2;   %(w2-cc2)'*QQ2*(w2-cc2);
            history.fval3(step) = fval3;%+Lambda3'*(-g3)+cc3'*QQ3*cc3+0.25*g3'*g3;       %(w3-cc3)'*QQ3*(w3-cc3);
            history.fval4(step) = fval4;%+Lambda4'*(-g4)+cc4'*QQ4*cc4+0.25*g4'*g4;  %(w4-cc4)'*QQ4*(w4-cc4);
            history.fval5(step) = fval5;%+Lambda5'*(-g5)+cc5'*QQ5*cc5+0.25*g5'*g5;  %(w5-cc5)'*QQ5*(w5-cc5);
            history.fval(step) = (w2-cc2)'*QQ2*(w2-cc2)+Lambda2'*(c2-g2)+rho/2*norm(c2-g2).^2;
            
            history.r_norm1(step)  = norm(Lambda1-l1)/rho;  history.r_norm2(step)  = norm(Lambda2-l2)/rho;
            history.r_norm3(step)  = norm(Lambda3-l3)/rho;  history.r_norm4(step)  = norm(Lambda4-l4)/rho;
            history.r_norm5(step)  = norm(Lambda5-l5)/rho;
            
             % dual 误差  
             if step == 1
                 history.s_norm1(step) = sqrt(1.25)*norm(z1); history.s_norm2(step) = sqrt(1.25)*norm(z2);
                 history.s_norm3(step) = sqrt(1.25)*norm(z3); history.s_norm4(step) = sqrt(1.25)*norm(z4);
                 history.s_norm5(step) = sqrt(1.25)*norm(z5);
             else 
                 history.s_norm1(step) = sqrt(1.25)*(norm(z1-tmp1));  history.s_norm2(step) = sqrt(1.25)*(norm(z2-tmp2)); % 5 * 0.25
                 history.s_norm3(step) = sqrt(1.25)*(norm(z3-tmp3));  history.s_norm4(step) = sqrt(1.25)*(norm(z4-tmp4));
                 history.s_norm5(step) = sqrt(1.25)*(norm(z5-tmp5));
             end
             history.eps_pri1(step) = 30*ABSTOL + RELTOL*max([norm(c1),norm(g1)]);         %sqrt();
             history.eps_pri2(step) = 30*ABSTOL + RELTOL*max([norm(c2),norm(g2)]);
             history.eps_pri3(step) = 30*ABSTOL + RELTOL*max([norm(c3),norm(g3)]);
             history.eps_pri4(step) = 30*ABSTOL + RELTOL*max([norm(c4),norm(g4)]);
             history.eps_pri5(step) = 30*ABSTOL + RELTOL*max([norm(c5),norm(g5)]);
             
             history.eps_dual1(step)= 30*ABSTOL + RELTOL*(norm(rho*Lambda1)); history.eps_dual2(step)= 30*ABSTOL + RELTOL*(norm(rho*Lambda2));
             history.eps_dual3(step)= 30*ABSTOL + RELTOL*(norm(rho*Lambda3)); history.eps_dual4(step)= 30*ABSTOL + RELTOL*(norm(rho*Lambda4));
             history.eps_dual5(step)= 30*ABSTOL + RELTOL*(norm(rho*Lambda5));
             
            if ~QUIET
                fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', step, ...
                history.r_norm3(step), history.eps_pri3(step), ...
                history.s_norm3(step), history.eps_dual3(step),norm(w1q-w1),norm(w2q-w2),norm(w3q-w3));
            end

            if (history.r_norm1(step) < history.eps_pri1(step) && history.r_norm2(step) < history.eps_pri2(step) && history.r_norm3(step) < history.eps_pri3(step) &&...
                    history.r_norm4(step) < history.eps_pri4(step) && history.r_norm5(step) < history.eps_pri5(step) && ...
                    history.s_norm1(step) < history.eps_dual1(step) && history.s_norm2(step) < history.eps_dual2(step) && ...
                    history.s_norm3(step) < history.eps_dual3(step) && history.s_norm4(step) < history.eps_dual4(step) && ...
                    history.s_norm5(step) < history.eps_dual5(step))
                    break;
            end
end 
toc(t_start);
T = 1:180;
L1 = [0 1];
L1 = repmat(L1,180,1);
L1 = mat2cell(L1,ones(180,1),2);
L1 = blkdiag(L1{:});
L1 = [zeros(180,180) L1];
L2 = [0 0 1];
L2 = repmat(L2,180,1);
L2 = mat2cell(L2,ones(180,1),3);
L2 = blkdiag(L2{:});
L2 = [zeros(180,180) L2];


figure(9);
plot3(c1,T,L1*w1);
set(gca, 'YDir','reverse');
xlabel('diatance(m)'),ylabel('time(min)'),zlabel('velocity(m/s)');
grid on;
figure(1);

plot3(v2*c2,T,L2*w2);
set(gca, 'YDir','reverse');
xlabel('diatance(m)'),ylabel('time(min)'),zlabel('velocity(m/s)');
grid on;
figure(2);


plot3(v2*c3,T,L2*w3);
set(gca, 'YDir','reverse');
xlabel('diatance(m)'),ylabel('time(min)'),zlabel('velocity(m/s)');
grid on;
figure(3);


plot3(v2*c4,T,L2*w4);
set(gca, 'YDir','reverse');
xlabel('diatance(m)'),ylabel('time(min)'),zlabel('velocity(m/s)');
set(get(gca,'XLabel'),'FontSize',16);%图上文字为8 point或小5号
set(get(gca,'YLabel'),'FontSize',16);
set(get(gca,'ZLabel'),'FontSize',16);
grid on;
figure(4);


plot3(v2*c5,T,L2*w5);
set(gca, 'YDir','reverse');
xlabel('diatance(m)'),ylabel('time(min)'),zlabel('velocity(m/s)');
grid on;
figure(5);


plot3(c1,T,L1*w1,v2*c2,T,L2*w2,v2*c3,T,L2*w3,v2*c4,T,L2*w4,v2*c5,T,L2*w5);
set(gca, 'YDir','reverse');
xlabel('diatance(m)'),ylabel('time(min)'),zlabel('velocity(m/s)');
legend('train1','train2','train3','train4','train5','Location','northeast');
grid on;
figure(6);



plot(T,c1,T,v2*c2,T,v2*c3,T,v2*c4,T,v2*c5);
xlabel('time'),ylabel('diatance');
grid on;
figure(7);

blk = zeros(40,60);
blk2 = eye(40);
blk3 = zeros(40,80);
T1 = 61:100;
zzz = [blk blk2 blk3];
plot(T1,zzz*c1,T1,zzz*v2*c2,T1,zzz*v2*c3,T1,zzz*v2*c4,T1,zzz*v2*c5);
grid on;


figure(8);
s = length(history.fval);
plot(1:s,history.fval);
xlabel('iteration number'),ylabel('\Psi'),
grid on;



