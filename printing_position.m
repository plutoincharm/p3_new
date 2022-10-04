close all;
clc;
clearvars -except val;

frame  = 20; 
% frame 20 taken for start of SSP1
% frame 55 midswing phase
% frame 75 is taken as end of ssp
%

BC      =  readmatrix("BC.xlsx"); 
omg     =  readmatrix("omg.xlsx");
%alpval  =  readmatrix("alp.xlsx");
%tht1=BC(1,frame);tht2=BC(2,frame);tht3=BC(3,frame);tht4=BC(4,frame);tht5=BC(5,frame);tht6=BC(6,frame);tht7=BC(7,frame);hx=BC(8,frame) ;hy=BC(9,frame);
%omg1 = 0.1;omg2 = 0;omg3 = 0.2;omg4 = 0.0;omg5 = 0.5;omg6 = 0.5;omg7 = 0.5;vhx =-0.8;vhy = 0.2;
%omg1 = omg(1,frame); omg2 =omg(2,frame); omg3 = omg(3,frame);  omg5 = omg(5,frame); omg6 = omg(6,frame); omg7 = omg(7,frame);vhx =  omg(8,frame); vhy = omg(9,frame);
%omg4 = 0;
% omg4 = omg(4,frame);
%hx=BC(10,frame) ;hy=BC(11,frame);
%vhx =  omg(10,frame); vhy = omg(11,frame);






%omg1 = vel(1,5);omg2 = vel(2,5); omg3 =  vel(3,5); omg4 = vel(4,5);  omg5 =  vel(5,5);
%omg6 =   vel(6,5);  omg7 =   vel(7,5); vhx=  vel(8,5); vhy =  vel(9,5);
% x0 = [tht1;tht2;tht3;tht4;tht5;tht6;tht7;hx;hy;omg1;omg2;omg3;omg4;omg5;omg6;omg7;vhx;vhy];





%omg1 = 0.1;omg2 = 0;omg3 = 0.2;omg4 = 0.0;omg5 = 0.5;omg6 = 0.5;omg7 = 0.5;vhx =-0.8;vhy = 0.2;
%omg1 = omg(1,55); omg2 =omg(2,55); omg3 = omg(3,55);  omg5 = omg(5,55); omg6 = omg(6,55); omg7 = omg(7,55);vhx =  omg(8,55); vhy = omg(9,55);
%omg4 = 0;
 %omg4 = omg(4,55);


%omg1 = vel(1,5);omg2 = vel(2,5); omg3 =  vel(3,5); omg4 = vel(4,5);  omg5 =  vel(5,5);
%omg6 =   vel(6,5);  omg7 =   vel(7,5); vhx=  vel(8,5); vhy =  vel(9,5);
% x0 = [tht1;tht2;tht3;tht4;tht5;tht6;tht7;hx;hy];





Lf = 0.1968;
L4  = 0.0841;
L1=  0.4041;
L2 = 0.4418;
L5 = 0.4418;
L3 = 0.4033;
L6 = 0.4033;
r4 = 0.0432;
r7 = 0.0442;
gamma71 =  1.2147;
gamma72 = 3.9843; 
gamma61 = 1.2104;
gamma62 = 3.9711;
r7t = 0.1115;
r7h = 0.0887;
r4t = 0.1114 ;
r4h = 0.0877 ;
r4c =   0.0873;%    0.0992;
gamma43 =  3.9381 %3.6666;
gamma2 = 0.7110;

%i = 20;
i = 1;
%while i < 76
while i < 56 
  
  %tht1=val(1,i);tht2=val(2,i);tht3=val(3,i);tht4=val(4,i);tht5=val(5,i);tht6=val(6,i);tht7=val(7,i);hx=val(8,i);hy=val(9,i);
 %tht1=BC(3,i);tht2=BC(4,i);tht3=BC(5,i);tht4=BC(6,i);tht5=BC(7,i);tht6=BC(8,i);tht7=BC(9,i);hx=BC(1,i);hy=BC(2,i); copx=BC(11,i);
 %copy=BC(10,i)
%  x0 = [tht1;tht2;tht3;tht4;tht5;tht6;tht7;hx;hy]
%tht3=val(1,i);tht4=val(2,i);hx=val(3,i);hy=val(4,i);
%hx=BC(10,frame+i);
%hy=BC(11,frame+i); 
tht1=val(3,i);tht2=val(4,i);tht3=val(5,i);tht4=val(6,i);tht5=val(7,i);tht6=val(8,i);%tht7=val(9,i);
hx=val(1,i);hy=val(2,i);    % copx=BC(10,i);copy=BC(11,i);




xm = hx  + Lf*cos(tht4); 
ym = hy  + Lf*sin(tht4); 
xal = hx + L4*cos(tht4 + gamma2);
yal = hy + L4*sin(tht4 + gamma2);
xt = hx  + L4*cos(tht4 + gamma2) + L3*cos(tht3)+ L2*cos(tht2) + L1*cos(tht1);
yt = hy  + L4*sin(tht4 + gamma2) + L3*sin( tht3)+ L2*sin( tht2) + L1*sin(tht1);
xkl = hx + L4*cos(tht4 + gamma2) + L3*cos(tht3);
ykl = hy + L4*sin(tht4 + gamma2) + L3*sin( tht3);
xhl = hx + L4*cos(tht4 + gamma2) + L3*cos(tht3)  + L2*cos(tht2);
yhl = hy + L4*sin(tht4 + gamma2) + L3*sin( tht3) + L2*sin( tht2);


%%%%
xkr = hx + L4*cos(tht4 + gamma2) + L3*cos(tht3) + L2*cos(tht2) + L5*cos(pi + tht5);
ykr = hy + L4*sin(tht4 + gamma2) + L3*sin( tht3) + L2*sin( tht2) + L5*sin(pi + tht5);
xar = hx + L4*cos(tht4 + gamma2) + L3*cos(tht3) + L2*cos(tht2) + L5*cos(pi + tht5) + L6*cos(pi + tht6);
yar = hy + L4*sin(tht4 + gamma2) + L3*sin( tht3) + L2*sin( tht2) + L5*sin(pi + tht5) + L6*sin(pi + tht6);
%xhr = hx + L4*cos(tht4 + gamma2) + L3*cos(tht3) + L2*cos(tht2) + L5*cos(pi + tht5) + L6*cos(pi + tht6);
%yhr = hy + L4*sin(tht4 + gamma2) + L3*sin( tht3) + L2*sin( tht2) + L5*sin(pi + tht5) + L6*sin(pi + tht6);





xcop = hx + L1*cos(pi + tht1) + L2*cos(pi + tht2) + L3*cos(pi + tht3) + r4*cos(pi + tht4) + r4c*cos(gamma43 + tht4);
ycop = hy + L1*sin(pi + tht1) + L2*sin(pi + tht2) + L3*sin(pi + tht3) + r4*sin(pi + tht4) +  r4c*sin(gamma43 + tht4);





axis([-0.5 2 -0.5 2]) 
%axis([-0.5 1.5 -0.5 0.5]) 
base =line([-1 2],[0.02 0.02],'LineWidth',1,'Color','black');
%pelvic =line([lhip_x(i) rhip_x(i)],[lhip_z(i) rhip_z(i)],'LineWidth',1,'Color','black');
  T=line([xhl  xt],[yhl  yt],'LineWidth',1,'Color','black');
    
   u=line([xhl xkl],[yhl ykl],'LineWidth',1,'Color','red');
   v=line([xkl xal],[ykl yal ],'LineWidth',1,'Color','red');
   w=line([xal hx],[yal hy],'LineWidth',1,'Color','red'); 
  x=line([xm hx],[ym hy],'LineWidth',1,'Color','red');


   %y=line([xfl xtl],[yfl ytl],'LineWidth',1,'Color','red');
   %lo =line([xal xhl],[yal yhl],'LineWidth',1,'Color','red');
%   ko=line([xal xtl],[yal ytl],'LineWidth',1,'Color','red');
 %  yhj =line([xhl  xtl],[yhl  ytl],'LineWidth',1,'Color','red'); 
%   yhj =line([xfl  xcop],[yfl  ycop],'LineWidth',1,'Color','red'); 
   hold on
   plot(xcop,ycop,'g*')
   % legend('Left leg')

   u1=line([xhl xkr],[yhl ykr],'LineWidth',1,'Color','blue'); 
   v1=line([xkr xar],[ykr yar],'LineWidth',1,'Color','blue'); 

    

   %w1=line([xar xfr],[yar yfr],'LineWidth',1,'Color','blue');
   %x1=line([xfr xhr],[yfr yhr],'LineWidth',1,'Color','blue');
   %y1=line([xfr xtr],[yfr ytr],'LineWidth',1,'Color','blue');
  %yh1 =line([xhr xtr],[yhr ytr],'LineWidth',1,'Color','blue');
  % lo1 =line([xar xhr],[yar yhr],'LineWidth',1,'Color','blue');
  % ko1 =line([xar xtr],[yar ytr],'LineWidth',1,'Color','blue');
   %legend('','Torso','Left leg','','','','','Right leg','','','','')
  % t1= text(1000,1500,txt,'Color','red','FontSize',14);

 i = i +1;
 pause
end  



