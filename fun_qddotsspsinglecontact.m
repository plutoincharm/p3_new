function qddot = fun_qddotsspsinglecontact(x,u,dt)
%%%% SSP phase- 1   %%%%%%%%%%%%%%


global A B Nx Nu pert MI L m  nx ny tx ty g r lam vars misc alp alpval indic kc lamall xdata lamx lamy val af acal fx fy Mmat2 invM phi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = m(1);
m2 = m(2);
m3 = m(3);
m4 = m(4);
m5 = m(5);
m6 = m(6);
m7 = m(7);

%L1 = L(1);
L2 = L(1);
L3 = L(2);
%L4 = L(4);
L5 = L(3);
L6 = L(4);

r1 = r(1);% trunk
r2 = r(2);
r3 = r(3);
%r4 = ra(4);
%r5 = r(4);
%r6 = r(5);
%r7 = ra(7);

MI1 = MI(1);
MI2 = MI(2);
MI3 = MI(3);
MI4 = MI(4);
MI5 = MI(5);
MI6 = MI(6);
MI7 = MI(7);




%M= 0;
d1 = 0.1081;
d3 = 0.0682;
h4 = 0.0437;
h2 = 0.0722;
r4c =  vars(5)   %0.1638;
% 0.0872 new value considering lycop = 0
gamma43 =  vars(6);% 2.3752 ;  
% 2.6481 new value considering lycop = 0



r14 = vars(1);
L4  = vars(2);
gamma1 = vars(3);
gamma2 = vars(4);
beta1 = vars(5);
beta2 = vars(6);
r74 = vars(7);
L7  = vars(8);

%{
tht1 = xdata(1);
tht2 = xdata(2);
tht3 = xdata(3);
tht4 = xdata(4);
tht5 = xdata(5);
tht6 = xdata(6);
tht7 = xdata(7);
x1   = xdata(8);
y1   = xdata(9);
omg1 = xdata(10);
omg2 = xdata(11);
omg3 = xdata(12);
omg4 = xdata(13);
omg5 = xdata(14);
omg6 = xdata(15);
omg7 = xdata(16);
vhx = xdata(17);
vhy = xdata(18);

%}
%%%%%%%%%%%%%%%%%%%%%%

%{
tht1 = x(3);
tht2 = x(4);
tht3 = x(5);
tht4 = x(6);
tht5 = x(7);
tht6 = x(8);
tht7 = x(9);
x1   = x(1);
y1   = x(2);
omg1 = x(12);
omg2 = x(13);
omg3 = x(14);
omg4 = x(15);
omg5 = x(16);
omg6 = x(17);
omg7 = x(18);
vhx = x(10);
vhy = x(11);


%}

tht1 = x(3);
tht2 = x(4);
tht3 = x(5);
tht4 = x(6);
x1   = x(1);
y1   = x(2);
omg1 = x(9);
omg2 = x(10);
omg3 = x(11);
omg4 = x(12);
vhx = x(7);
vhy = x(8);




%{
alp1 = alp(3);
alp2 = alp(4);
alp3 = alp(5);
alp4 = alp(6);
alp5 = alp(7);
alp6 = alp(8);
alp7 = alp(9);
ax1 =  alp(1);
ay1 =  alp(2);
%}
T1  =  u(1);
T2  =  u(2);
T3  =  u(3);
T4  =  u(4);
%T5  =  u(5);
%T6  =  u(6);
%T7  =  u(7);
c = -1.6;
%T1  = u(1);
%{
T2  = c*omg2*r2 +u(1);
T3  = c*omg3*r3 +u(2);
T4  = c*omg4*r4 +u(3);
T5  = c*omg5*r5 +u(4);
T6  = c*omg6*r6 +u(5);
T7  = c*omg7*r7 +u(6);
%}
lam1 = lamy;
lam2 =lamx;
%{
T2  = -c*omg2*r2 +u(1);
T3  = -c*omg3*r3 +u(2);
T4  = -MI4*alp4  -c*omg4*r4 +u(3);
T5  = -c*omg5*r5 +u(4);
T6  = -c*omg6*r6 +u(5);
T7  = -c*omg7*r7 +u(6);

%}

%T4  =  c*omg4 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding external forces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat=[m1+m2+m3+m4,0,(-1).*m1.*r1.*sin(tht1),(-1).*L2.*m1.*sin(tht2)+( ...
  -1).*m2.*r2.*sin(tht2),(-1).*L3.*m1.*sin(tht3)+(-1).*L3.*m2.*sin( ...
  tht3)+(-1).*m3.*r3.*sin(tht3),(-1).*m4.*r14.*sin(gamma1+tht4)+(-1) ...
  .*L4.*m1.*sin(gamma2+tht4)+(-1).*L4.*m2.*sin(gamma2+tht4)+(-1).* ...
  L4.*m3.*sin(gamma2+tht4);0,m1+m2+m3+m4,m1.*r1.*cos(tht1),L2.*m1.* ...
  cos(tht2)+m2.*r2.*cos(tht2),L3.*m1.*cos(tht3)+L3.*m2.*cos(tht3)+ ...
  m3.*r3.*cos(tht3),m4.*r14.*cos(gamma1+tht4)+L4.*m1.*cos(gamma2+ ...
  tht4)+L4.*m2.*cos(gamma2+tht4)+L4.*m3.*cos(gamma2+tht4);(-1).*m1.* ...
  r1.*sin(tht1),m1.*r1.*cos(tht1),MI1+m1.*r1.^2,L2.*m1.*r1.*cos( ...
  tht1+(-1).*tht2),L3.*m1.*r1.*cos(tht1+(-1).*tht3),L4.*m1.*r1.*cos( ...
  gamma2+(-1).*tht1+tht4);(-1).*(L2.*m1+m2.*r2).*sin(tht2),(L2.*m1+ ...
  m2.*r2).*cos(tht2),L2.*m1.*r1.*cos(tht1+(-1).*tht2),L2.^2.*m1+MI2+ ...
  m2.*r2.^2,L2.*L3.*m1.*cos(tht2+(-1).*tht3)+L3.*m2.*r2.*cos(tht2+( ...
  -1).*tht3),L2.*L4.*m1.*cos(gamma2+(-1).*tht2+tht4)+L4.*m2.*r2.* ...
  cos(gamma2+(-1).*tht2+tht4);(-1).*(L3.*(m1+m2)+m3.*r3).*sin(tht3), ...
  (L3.*(m1+m2)+m3.*r3).*cos(tht3),L3.*m1.*r1.*cos(tht1+(-1).*tht3), ...
  L2.*L3.*m1.*cos(tht2+(-1).*tht3)+L3.*m2.*r2.*cos(tht2+(-1).*tht3), ...
  L3.^2.*m1+L3.^2.*m2+MI3+m3.*r3.^2,L3.*L4.*m1.*cos(gamma2+(-1).* ...
  tht3+tht4)+L3.*L4.*m2.*cos(gamma2+(-1).*tht3+tht4)+L4.*m3.*r3.* ...
  cos(gamma2+(-1).*tht3+tht4);(-1).*m4.*r14.*sin(gamma1+tht4)+(-1).* ...
  L4.*(m1+m2+m3).*sin(gamma2+tht4),m4.*r14.*cos(gamma1+tht4)+L4.*( ...
  m1+m2+m3).*cos(gamma2+tht4),L4.*m1.*r1.*cos(gamma2+(-1).*tht1+ ...
  tht4),L2.*L4.*m1.*cos(gamma2+(-1).*tht2+tht4)+L4.*m2.*r2.*cos( ...
  gamma2+(-1).*tht2+tht4),L3.*L4.*m1.*cos(gamma2+(-1).*tht3+tht4)+ ...
  L3.*L4.*m2.*cos(gamma2+(-1).*tht3+tht4)+L4.*m3.*r3.*cos(gamma2+( ...
  -1).*tht3+tht4),L4.^2.*m1+L4.^2.*m2+L4.^2.*m3+MI4+m4.*r14.^2]





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mmat'
Mmat - Mmat'
issymmetric(Mmat)
eig(Mmat)



%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% phi is for finding qddot with ext lambda calculation 
phi=[lam2+m1.*omg1.^2.*r1.*cos(tht1)+L2.*m1.*omg2.^2.*cos(tht2)+m2.* ...
  omg2.^2.*r2.*cos(tht2)+L3.*m1.*omg3.^2.*cos(tht3)+L3.*m2.* ...
  omg3.^2.*cos(tht3)+m3.*omg3.^2.*r3.*cos(tht3)+m4.*omg4.^2.*r14.* ...
  cos(gamma1+tht4)+L4.*m1.*omg4.^2.*cos(gamma2+tht4)+L4.*m2.* ...
  omg4.^2.*cos(gamma2+tht4)+L4.*m3.*omg4.^2.*cos(gamma2+tht4),lam1+ ...
  m1.*omg1.^2.*r1.*sin(tht1)+L2.*m1.*omg2.^2.*sin(tht2)+m2.* ...
  omg2.^2.*r2.*sin(tht2)+L3.*m1.*omg3.^2.*sin(tht3)+L3.*m2.* ...
  omg3.^2.*sin(tht3)+m3.*omg3.^2.*r3.*sin(tht3)+m4.*omg4.^2.*r14.* ...
  sin(gamma1+tht4)+L4.*m1.*omg4.^2.*sin(gamma2+tht4)+L4.*m2.* ...
  omg4.^2.*sin(gamma2+tht4)+L4.*m3.*omg4.^2.*sin(gamma2+tht4),T1+( ...
  -1).*L2.*m1.*omg2.^2.*r1.*sin(tht1+(-1).*tht2)+(-1).*L3.*m1.* ...
  omg3.^2.*r1.*sin(tht1+(-1).*tht3)+L4.*m1.*omg4.^2.*r1.*sin(gamma2+ ...
  (-1).*tht1+tht4),(-1).*T1+T2+L2.*m1.*omg1.^2.*r1.*sin(tht1+(-1).* ...
  tht2)+(-1).*L2.*L3.*m1.*omg3.^2.*sin(tht2+(-1).*tht3)+(-1).*L3.* ...
  m2.*omg3.^2.*r2.*sin(tht2+(-1).*tht3)+L2.*L4.*m1.*omg4.^2.*sin( ...
  gamma2+(-1).*tht2+tht4)+L4.*m2.*omg4.^2.*r2.*sin(gamma2+(-1).* ...
  tht2+tht4),(-1).*T2+T3+L3.*m1.*omg1.^2.*r1.*sin(tht1+(-1).*tht3)+ ...
  L2.*L3.*m1.*omg2.^2.*sin(tht2+(-1).*tht3)+L3.*m2.*omg2.^2.*r2.* ...
  sin(tht2+(-1).*tht3)+L3.*L4.*m1.*omg4.^2.*sin(gamma2+(-1).*tht3+ ...
  tht4)+L3.*L4.*m2.*omg4.^2.*sin(gamma2+(-1).*tht3+tht4)+L4.*m3.* ...
  omg4.^2.*r3.*sin(gamma2+(-1).*tht3+tht4),(-1).*T3+T4+(-1).*L4.* ...
  m1.*omg1.^2.*r1.*sin(gamma2+(-1).*tht1+tht4)+(-1).*L2.*L4.*m1.* ...
  omg2.^2.*sin(gamma2+(-1).*tht2+tht4)+(-1).*L4.*m2.*omg2.^2.*r2.* ...
  sin(gamma2+(-1).*tht2+tht4)+(-1).*L3.*L4.*m1.*omg3.^2.*sin(gamma2+ ...
  (-1).*tht3+tht4)+(-1).*L3.*L4.*m2.*omg3.^2.*sin(gamma2+(-1).*tht3+ ...
  tht4)+(-1).*L4.*m3.*omg3.^2.*r3.*sin(gamma2+(-1).*tht3+tht4)]



%{
%%%%%%%%%%%%%%%%%%%%%%%
inv(Mmat);
qddot = Mmat\phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rhs1 = rhs;
%Cn  = Cmat;
%rhs2 = Gd';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%rhs2 -gn
%%rhs1 - phi
%Mmat
invM = inv(Mmat)

%Gmat = Cn*invM*Cn';
%invG = inv(Gmat);

%nr  = Gd' - Cn*invM*rhs1';
%P1 = Gd' ;
%P2 = Cn*invM*rhs1';
%lam = Gmat\( Gd'- Cn*invM*rhs1')

%lam = [lam1;lam2;lam3;lam4]
%(Cn'*lam + rhs1')
%qddot_invdy = invM*(Cn'*lam + rhs1');
%}

% qddot = qddot_invdy

phi;
qddot_invdy = invM*phi'
%qddot_invdy(4) = 0
af = qddot_invdy ;
%qddot = [omg1;omg2;omg3;omg4;omg5;omg6;omg7;vhx;vhy;qddot_invdy]
 qddot = qddot_invdy;
 %lamx = lam2 + lam4;
 %lamy = lam1 + lam3;
% lamx = lamt;
% lamy = lamn;

end






 






