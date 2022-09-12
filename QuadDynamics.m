function [sys,x0,str,ts] = QuadDynamicsC(t,x,u,flag,quad)
%CSFUNC An example MATLAB file S-function for defining a continuous system.  
%   Example MATLAB file S-function implementing continuous equations: 
%      x' = Ax + Bu
%      sys  = Cx + Du
%   
%   See sfuntmpl.m for a general S-function template.
%
%   See also SFUNTMPL.
  
%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $

Pos0 = [0 0 0];          
Att0 = [0 0 0];             
Vel0 = [0 0 0];              
Ang0 = [0 0 0];             
InitSys = [Pos0 Att0 Vel0 Ang0];

ro=1.184; % air density
r=0.165;%0.127;  %prop diameter
Ct=0.0048; %0.001396; 
Cp=Ct*sqrt(Ct/2);
S=r^2*pi;
b=Ct*ro*S*r.^2;
k=Cp*ro*S*r.^3;

m=1.284;%1.282;
g=9.81;

b=0.9937e-8;
k=b/26;
%b=1.3234e-005;%9.7575e-6; 
%k=1.0697e-007%(9.7575e-6)/30;

l=0.254;

Ixx = 8.283e-2;
Iyy = 8.283e-2;
Izz = 8.283e-2;
I = [Ixx 0 0; 
    0 Iyy 0; 
    0 0 Izz];

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(InitSys,quad);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(x,u,quad);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(x);

  %%%%%%%%%%%%%%%%%%%
  % Unhandled flags %
  %%%%%%%%%%%%%%%%%%%
  case { 2, 4, 9 },
    sys = [];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
end
% end csfunc

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(InitSys,quad)

sizes = simsizes;
sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 4;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);
x0  = InitSys;
str = [];
ts  = [0 0];

end
% end mdlInitializeSizes
%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(x,u,quad)

Pos=x(1:3);
Att=x(4:6);
Vel=x(7:9);
Ang=x(10:12);

phi=Att(1);
the=Att(2);
psi=Att(3);


R = [cos(the)*cos(psi)  sin(phi)*sin(the)*cos(psi)-cos(phi)*sin(psi)  cos(phi)*sin(the)*cos(psi)+sin(phi)*sin(psi);
     cos(the)*sin(psi)  sin(phi)*sin(the)*sin(psi)+cos(phi)*cos(psi)  cos(phi)*sin(the)*sin(psi)-sin(phi)*cos(psi);
     -sin(the)  sin(phi)*cos(the)  cos(phi)*cos(the)];

W =[-sin(the) 0 1;
    cos(the)*sin(phi) cos(phi) 0;
    cos(the)*cos(phi) -sin(phi) 0];
      
%----------------------------------------                                
%for i=1:4
    %w2(i)=u(i)^2;
%end;
F=quad.b*sum(u);
tau=[quad.l*quad.b*(-u(1)+u(2)+u(3)-u(4))
     quad.l*quad.b*(-u(1)-u(2)+u(3)+u(4))
     quad.k*(-u(1)+u(2)-u(3)+u(4))];

%----------------------------------------

dPos=Vel; 
dAtt=W*Ang; %(W*Ang???)
dVel=R/quad.m*[0;0;F]-[0;0;quad.g]; 
%[cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi); 
%sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi); 
%cos(the)*cos(phi)]
dAng=inv(quad.J)*inv(W)*(-cross(Ang,quad.J*Ang)+tau); %(-I*dW*dAtt)

sys = [dPos' dAtt' dVel' dAng'];
end

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(x)

if x(3)<0 
    x(3)=0;
 %   if x(6) > 0
  %      x(6) = 0;
  % end
end

Att=x(4:6);
phi=Att(1);
the=Att(2);
psi=Att(3);

R = [cos(the)*cos(psi)  sin(phi)*sin(the)*cos(psi)-cos(phi)*sin(psi)  cos(phi)*sin(the)*cos(psi)+sin(phi)*sin(psi);
     cos(the)*sin(psi)  sin(phi)*sin(the)*sin(psi)+cos(phi)*cos(psi)  cos(phi)*sin(the)*sin(psi)-sin(phi)*cos(psi);
     -sin(the)  sin(phi)*cos(the)  cos(phi)*cos(the)];

W =[-sin(the) 0 1;
    cos(the)*sin(phi) cos(phi) 0;
    cos(the)*cos(phi) -sin(phi) 0];

sys = [x(1:6); inv(R)*x(7:9);  W*x(10:12)];

end
% end mdlOutputs