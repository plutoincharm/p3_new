clc; close all

clearvars -except GRFx GRFy GRFxr GRFyr GRFxl GRFyl; 
global A B Nx Nu pert MI L m  nx ny tx ty g r lam vars misc alp alpval indic kc lamall xdata lamx lamy val af acal fx fy Mmat2 invM phi

g = 9.81; % gravity
Nx = 12;
Nu  = 4;
Tf = 0.075;
%Tf = 1;
dt = 0.001;
Nt = round(Tf/dt)+1;
A = zeros(Nx,Nx);
B = zeros(Nx,Nu);
pert = 0.0001;
nx = 0;tx = 1;
ny = 1;ty = 0;
%s_frame = 20;      
s_frame = 1; 
e_frame = 76; 


dynamics_midpoint = @(x,u,dt) x + fun_xdot((x + fun_xdot(x,u,dt)*dt/2),u,dt)*dt;

vars = fun_data(); 
                                             
% initial conditions

BC      =  readmatrix("BC.xlsx");
omg     =  readmatrix("omg.xlsx");

%%%% extract the heel data

tht1=BC(3,s_frame);tht2=BC(4,s_frame);tht3=BC(5,s_frame);tht4=BC(6,s_frame);tht5=BC(7,s_frame);tht6=BC(8,s_frame);tht7=BC(9,s_frame);hx=BC(1,s_frame) ;hy=BC(2,s_frame);
omg1 = omg(3,s_frame); omg2 =omg(4,s_frame); omg3 = omg(5,s_frame);  omg5 = omg(7,s_frame); omg6 = omg(8,s_frame); omg7 = omg(9,s_frame);vhx =  omg(1,s_frame); vhy = omg(2,s_frame);
omg4 = omg(6,s_frame);
%x0 = [hx;hy;tht1;tht2;tht3;tht4;tht5;tht6;tht7;vhx;vhy;omg1;omg2;omg3;omg4;omg5;omg6;omg7];
x0 = [hx;hy;tht1;tht2;tht3;tht4;vhx;vhy;omg1;omg2;omg3;omg4];


% goal
thtf1=BC(3,e_frame);thtf2=BC(4,e_frame);thtf3=BC(5,e_frame);thtf5=BC(7,e_frame);thtf6=BC(8,e_frame);thtf7=BC(9,e_frame);hfx=BC(1,e_frame);hfy=BC(2,e_frame);
omgf1 = omg(3,e_frame);omgf2 = omg(4,e_frame); omgf3 =  omg(5,e_frame); omgf5 = omg(7,e_frame);
omgf6 =   omg(8,e_frame);  omgf7 =  omg(9,e_frame); vhfx=  omg(1,e_frame); vhfy =  omg(2,e_frame);
omgf4 = omg(6,e_frame); 
thtf4=BC(6,e_frame);
%xf = [hfx;hfy;thtf1;thtf2;thtf3;thtf4;thtf5;thtf6;thtf7;vhfx;vhfy;omgf1;omgf2;omgf3;omgf4;omgf5;omgf6;omgf7];
xf = [hfx;hfy;thtf1;thtf2;thtf3;thtf4;vhfx;vhfy;omgf1;omgf2;omgf3;omgf4];

% costs
%Q = 1e-5*eye(4);
Q =  1e-5*eye(Nx);
Qf = 25*eye(Nx);
R = 5*1e-7*eye(Nu);
I = eye(Nu);
e_dJ = 1e-12;



% initialization
u = rand(Nu,Nt-1)*20;
u = ones(Nu,Nt-1)*20;
%u = zeros(Nu,Nt-1);
uf  = ones(Nu,Nt-1);
x = zeros(Nx,Nt);
x_prev = zeros(Nx,Nt);
x(:,1) = x0;

alp = zeros(Nx/2);
alp =alpval(1:6,s_frame);
lamx = fx(s_frame);
lamy =  fy(s_frame);
afall = [alp];

% first roll-out
for k = 2:Nt-1
           x(:,k) = dynamics_midpoint(x(:,k-1),u(:,k-1),dt);
           alp =  alpval(:,s_frame+k-1);

      %{
        alp(1) =af(1);   
        alp(2) =af(2);     
        alp(3) =af(3);   
        alp(4) =af(4);
        alp(5) =af(5);   
        alp(6) =af(6);     
        alp(7) =af(7);   
        alp(8) =af(8);
        alp(9) =af(9);
      %}

       lamx = fx(s_frame+k-1);
        lamy = fy(s_frame+k-1);
        afall = [afall,af];
        k;
        pause()           
        
end




val = x;
tdt = 0.2:0.01:0.74;
tt = 0.23:dt:1.58;
fr = 1:1:125;
cr = s_frame:1:76;
hr = s_frame:1:75;


%{

%%%  lamx
figure;
%plot(tt,GRFxl,'b-','LineWidth',1);
%grid on;
%hold on;
%plot(tt,GRFx,'r-','LineWidth',1);
%grid on;
%hold on;
plot(tdt,lamallx,'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_x (N) \rightarrow');
legend('2D Lagrangian','Experimental');

 
%%%  lamy
figure;
%plot(tt,GRFyl,'b-','LineWidth',1);
%grid on;
%hold on;
%plot(tt,GRFy,'r-','LineWidth',1);
%grid on;
%hold on;
plot(tdt,lamally,'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('GRF_y (N) \rightarrow');
legend('2D Lagrangian','Experimental');


%}

%{
figure;
%hold on;
plot(fr,alpval(4,:),'g-','LineWidth',1);
%grid on;
hold on;
plot(cr,afall,'b-','LineWidth',1);
%}


%%% alpha x
figure;
plot(hr,afall(1,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(1,s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('ax (N) \rightarrow');
legend('calc','dataset');



%%% alpha y
figure;
plot(hr,afall(2,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(2, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('ay (N) \rightarrow');
legend('calc','dataset');

%%% alpha 1
figure;
plot(hr,afall(3,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(3,s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a1 (N) \rightarrow');
legend('calc','dataset');

%%% alpha 2
figure;
plot(hr,afall(4,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(4, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a2 (N) \rightarrow');
legend('calc','dataset');


%%% alpha 3
figure;
plot(hr,afall(5,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(5, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a3 (N) \rightarrow');
legend('calc','dataset');



%%% alpha 4
figure;
plot(hr,afall(6,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(6, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a4 (N) \rightarrow');
legend('calc','dataset');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%left leg


%%% theta3
figure;
plot(cr,val(5,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(5, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta3 (N) \rightarrow');
legend('calc','dataset');


%%% theta4
figure;
plot(cr,val(6,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(6, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta4 (N) \rightarrow');
legend('calc','dataset');


%%% theta3
figure;
plot(cr,val(5,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(5, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta3 (N) \rightarrow');
legend('calc','dataset');


%%% theta2
figure;
plot(cr,val(4,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(4, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta2 (N) \rightarrow');
legend('calc','dataset');


%%% theta1
figure;
plot(cr,val(3,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(3, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta1 (N) \rightarrow');
legend('calc','dataset');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Right leg


%%% theta5
figure;
plot(cr,val(7,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(7,s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta5 (N) \rightarrow');
legend('calc','dataset');


%%% theta6
figure;
plot(cr,val(8,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(8, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta6 (N) \rightarrow');
legend('calc','dataset');


%{
%%% theta7
figure;
plot(hr,val(9,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(9,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta7 (N) \rightarrow');
legend('calc','dataset');

%}




%%% hx
figure;
plot(cr,val(1,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(1,s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hx (N) \rightarrow');
legend('calc','dataset');


%%% hy
figure;
plot(cr,val(2,1:e_frame),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(2, s_frame:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hy (N) \rightarrow');
legend('calc','dataset');

%{
%%% hax
figure;
plot(hr,val(3,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(10,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hax (N) \rightarrow');
legend('calc','dataset');


%%% hay
figure;
plot(hr,val(4,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(11,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hay (N) \rightarrow');
legend('calc','dataset');
%}



 
% original cost
J = 0;                                              
for k = 1:Nt-1
    J = J +  (x(:,k)-xf)'*Q*( x(:,k)-xf) + (u(:,k))'*R*(u(:,k)) 
end 
disp('Original cost:')
J = 0.5*(J + (x(:,Nt)-xf)'*Qf*(x(:,Nt)-xf)) 
val = x;

%{

pause()
disp(' ILQR starts--------------------------------------- ')

%%%%%%%%%%%%%%%% ILQR Algorithm  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   
p = ones(Nx,Nt);
P = zeros(Nx,Nx,Nt);
%d = ones(Nu,Nu,Nt-1);
d = ones(Nu,Nt-1);
K = zeros(Nu,Nx,Nt-1);
%pdim = ones(Nx,Nu,Nt);
dJ = 0.0;  % change in cost

xn = zeros(Nx,Nt);
un = zeros(Nu,Nt-1);
% func g(dx,du) is perturbation of val func
% grad- g/ hessian-G of change in value fun
gx = zeros(Nx);
gu = zeros(Nu);
Gxx = zeros(Nx,Nx);
Guu = zeros(Nu,Nu);
Gxu = zeros(Nx,Nu);
Gux = zeros(Nu,Nx);

iter = 0;
while max(abs(d(:))) >  1e-3
    
    iter = iter +  1 

 %%%%% Backward Pass %%%%%
    dJ = 0.0;
    p(:,Nt) = Qf*(x(:,Nt)-xf);     %%% P is vx
    P(:,:,Nt) = Qf;                %%% P is vxx
    mu_reg = 0;
    for k = (Nt-1):-1:1
   
          %Calculate derivatives of stage cost
           q = Q*( x(:,k)-xf);     % lx
           r = R*u(:,k);           % lu
            
           
           % select alpha value for finding  A and B for each config at K

            alp =  alpval(:,s_frame+k-1);
            lamx =fx(s_frame+k-1)
            lamy = fy(s_frame+k-1)


            A = fun_amat(x(:,k),u(:,k),dt);
            B = fun_bmat(x(:,k),u(:,k),dt);

           %gradient of change in val fn
            gx = q + A'*p(:,k+1);% gx = dg/dx  
            gu = r + B'*p(:,k+1)   % gu = dg/du
    
          %iLQR (Gauss-Newton) version
          %Hessian
             Gxx = Q + A'*(P(:,:,k+1))*A;
             Guu = R + B'*(P(:,:,k+1)+ mu_reg*eye(Nx))*B
             Gxu = A'*P(:,:,k+1)*B;
             Gux = B'*(P(:,:,k+1) + mu_reg*eye(Nx))*A;     
             
             %beta = 0.1;
             log = issymmetric([Guu]);
             eigv = eig([Guu]);

          if any(eig(Guu)<0)
            mu_reg = mu_reg + 1;
            k = Nt-1;
            disp('regularized')
          end
          %% 
        %{
              while (log==0) || all(eigv < 0) 
                    Gxx = Gxx + A'*beta*I*A
                    Guu = Guu + B'*beta*I*B
                    Gxu = Gxu + A'*beta*I*B
                    Gux = Gux + B'*beta*I*A
                    beta = 2*beta
                    %display("regularizing G")
                    display(beta)
                    log = issymmetric([Gxx Gxu; Gux Guu]);
                    eigv = eig([Gxx Gxu; Gux Guu]);
              end
         %}
            d(:,k) = Guu\gu;  % feedforward term
            K(:,:,k) = Guu\Gux; % feedback gain term
    
             p(:,k) = gx - K(:,:,k)'*gu + K(:,:,k)'*Guu*d(:,k) - Gxu*d(:,k);
             P(:,:,k) = Gxx + K(:,:,k)'*Guu*K(:,:,k) - Gxu*K(:,:,k) - K(:,:,k)'*Gux;
             dJ = dJ +  gu'*d(:,k);
 disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS')
       k
       iter
    
      %pause()
    end
    disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS- cOMPLETED')
         % pause()
    
  
%%%% End of Backward Pass %%%%%
     %alp = zeros(Nx/2);%
   % alp =alpval(:,frame) 
    %alp(4) = 0; 
    %Forward rollout with line search
    xn(:,1) = x(:,1);
    alpha = 1.0;
    indic = 15;
    xdata = x0;   
     
         %  alp(1) =  alpval(3,frame);
          % alp(2) =  alpval(4,frame);
          % alp(3) =  alpval(10,frame);
           %alp(4) =  alpval(11,frame);  
            lamx = [];
            lamallx = [lamx];
            lamy = [];
            lamally = [lamy];
            lamx =fx(s_frame)
            lamy = fy(s_frame)
   for kc = 1:(Nt-1)
        un(:,kc) = u(:,kc) - alpha*d(:,kc) - (K(:,:,kc)*(xn(:,kc)-x(:,kc)));
        xn(:,kc+1) = dynamics_midpoint(xn(:,kc),un(:,kc),dt);
        lamx =fx(s_frame+kc-1)
        lamy = fy(s_frame+kc-1)

        %alp(1) =  alpval(3,frame+kc);
        %alp(2) =  alpval(4,frame+kc);
        %alp(3) =  alpval(10,frame+kc);
        %alp(4) =  alpval(11,frame+kc);
        lamallx = [lamallx,lamx];
        lamally = [lamally,lamy];
        afall = [afall,af];
        %thetval = BC(:,frame+kc);
        %omgval = omg(:,frame+kc);
        %xdata = [thetval;omgval]; 
        %alp = zeros(Nx/2);
       % alp =  alpval(:,frame+kc);
       % alp(4) = 0;
        disp('ITERATIPONSSSSSSSSSSSSSSSSSSSSS')
       kc
     %pause()
    end
    disp('EOFP')
     %pause() 
    Jn = 0;
    for k = 1:Nt-1
        Jn = Jn + (xn(:,k)-xf)'*Q*(xn(:,k)-xf) + un(:,k)'*R*un(:,k);
    end
   Jn = 0.5*(Jn + (xn(:,Nt)-xf)'*Qf*(xn(:,Nt)-xf))
    
     disp('line search')
     %pause() 
 liter =1;
    while isnan(Jn) || Jn > (J - 1e-2*alpha*dJ)
         disp('Inside search')
        alpha = 0.5*alpha
         %xdata = x0;   
         %alp =alpval(:,frame)
       
           %alp(1) =  alpval(3,frame);
           %alp(2) =  alpval(4,frame);
           %alp(3) =  alpval(10,frame);
          % alp(4) =  alpval(11,frame); 
           lamx = [];
           lamallx = [lamx];
           lamy = [];
           lamally = [lamy];
           lamx =fx(s_frame)
           lamy = fy(s_frame)

        for k = 1:(Nt-1)
            un(:,k) = u(:,k) - alpha*d(:,k) - (K(:,:,k)*(xn(:,k)-x(:,k)));
            xn(:,k+1) = dynamics_midpoint(xn(:,k),un(:,k),dt);
            lamx = fx(s_frame+k-1);
            lamy = fy(s_frame+k-1);
            %alp(1) =  alpval(3,frame+k);
           % alp(2) =  alpval(4,frame+k);
           % alp(3) =  alpval(10,frame+k);
            %alp(4) =  alpval(11,frame+k);
            lamallx = [lamallx,lamx];
            lamally = [lamally,lamy];
            afall = [afall,af];
            %thetval = BC(:,frame+k);
            %omgval = omg(:,frame+k);
            %xdata = [thetval;omgval];
           % alp =  alpval(:,frame+k);
            %xn(:,k+1) = dynamics_rk4(xn(:,k),un(k)
            liter = liter + 1;
        end
        %Jn = cost(xn,un);
        Jn = 0;
        for k = 1:Nt-1
            Jn = Jn + (xn(:,k) - xf)'*Q*(xn(:,k) - xf) + un(:,k)'*R*un(:,k);
        end
     Jn = 0.5*(Jn + (xn(:,Nt) - xf)'*Qf*(xn(:,Nt) - xf))
     % pause()
    end
     %  pause() 
 
   % J = Jn;
    x = xn;
    u = un;
   %if iter > 5
    %   break
    %end
  end

val = x;

close all


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%



%%% alpha x
figure;
plot(hr,afall(1,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(1,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('ax (N) \rightarrow');
legend('calc','dataset');



%%% alpha y
figure;
plot(hr,afall(2,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(2,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('ay (N) \rightarrow');
legend('calc','dataset');

%%% alpha 1
figure;
plot(hr,afall(3,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(3,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a1 (N) \rightarrow');
legend('calc','dataset');

%%% alpha 2
figure;
plot(hr,afall(4,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(4,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a2 (N) \rightarrow');
legend('calc','dataset');


%%% alpha 3
figure;
plot(hr,afall(5,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(5,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a3 (N) \rightarrow');
legend('calc','dataset');



%%% alpha 4
figure;
plot(hr,afall(6,:),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,alpval(6,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('a4 (N) \rightarrow');
legend('calc','dataset');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%left leg


%%% theta3
figure;
plot(hr,val(5,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(5,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta3 (N) \rightarrow');
legend('calc','dataset');


%%% theta4
figure;
plot(hr,val(6,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(6,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta4 (N) \rightarrow');
legend('calc','dataset');


%%% theta3
figure;
plot(hr,val(5,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(5,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta3 (N) \rightarrow');
legend('calc','dataset');


%%% theta2
figure;
plot(hr,val(4,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(4,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta2 (N) \rightarrow');
legend('calc','dataset');


%%% theta1
figure;
plot(hr,val(3,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(3,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta1 (N) \rightarrow');
legend('calc','dataset');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Right leg


%%% theta5
figure;
plot(hr,val(7,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(7,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta5 (N) \rightarrow');
legend('calc','dataset');


%%% theta6
figure;
plot(hr,val(8,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(8,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta6 (N) \rightarrow');
legend('calc','dataset');


%%% theta7
figure;
plot(hr,val(9,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(9,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('theta7 (N) \rightarrow');
legend('calc','dataset');



%%% hx
figure;
plot(hr,val(1,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(1,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hx (N) \rightarrow');
legend('calc','dataset');


%%% hy
figure;
plot(hr,val(2,1:56),'b-','LineWidth',1);
grid on;
hold on;
plot(cr,BC(2,20:76),'g-','LineWidth',1);
xlabel('time (sec) \rightarrow');
ylabel('hy (N) \rightarrow');
legend('calc','dataset');




%%% Torque %%%%%%%%%%%%
%%%%%%%%%%

%%% T1
figure;
plot(hr,un(1,1:56),'b-','LineWidth',1);
grid on;
xlabel('time (sec) \rightarrow');
ylabel('u1 (N) \rightarrow');
legend('calc');
%%% T2
figure;
plot(hr,un(2,1:56),'b-','LineWidth',1);
grid on;
xlabel('time (sec) \rightarrow');
ylabel('u2 (N) \rightarrow');
legend('calc');
%%% T3
figure;
plot(hr,un(3,1:56),'b-','LineWidth',1);
grid on;
xlabel('time (sec) \rightarrow');
ylabel('u3 (N) \rightarrow');
legend('calc');
%%% T4
figure;
plot(hr,un(4,1:56),'b-','LineWidth',1);
grid on;
xlabel('time (sec) \rightarrow');
ylabel('u4(N) \rightarrow');
legend('calc');
%%% T5
figure;
plot(hr,un(5,1:56),'b-','LineWidth',1);
grid on;
xlabel('time (sec) \rightarrow');
ylabel('u5 (N) \rightarrow');
legend('calc');
%%% T6
figure;
plot(hr,un(6,1:56),'b-','LineWidth',1);
grid on;
xlabel('time (sec) \rightarrow');
ylabel('u6(N) \rightarrow');
legend('calc');

%%% T7
figure;
plot(hr,un(7,1:56),'b-','LineWidth',1);
grid on;
xlabel('time (sec) \rightarrow');
ylabel('u7(N) \rightarrow');
legend('calc');

%}

