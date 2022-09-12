function [sys,x0,str,ts] = ForcesTest(t,x,u,flag)
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

T=zeros(12,1);

m=1.282;
g=9.81;

ro=1.2041; % air density
r=0.127;  %prop diameter
Ct=0.001396; Cp=0.000581; %Cp=Ct*sqrt(Ct/2);
A=r^2*pi;

%b=Ct*ro*A*r.^2;
%k=Cp*ro*A*r.^3;
b=9.7575e-6; 
k=(9.7575e-6)/30;

l=0.3;

Ixx = 0.08283;
Iyy = 0.08283;
Izz = 0.149;
J = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(InitSys);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u,l,b,k,m,g,J);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u,b);

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
function [sys,x0,str,ts]=mdlInitializeSizes(InitSys)

sizes = simsizes;
sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 27;
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
function sys=mdlDerivatives(~,x,u,l,b,k,m,g,J)

Pos=x(1:3);
Att=x(4:6);
Vel=x(7:9);
Ang=x(10:12);

w2=([u(1);u(2);u(3);u(4)]).^2;

Phi=Att(1); Theta=Att(2); Psi=Att(3);

A=[0 -l*b 0 l*b; l*b 0 -l*b 0; k -k k -k];

Rx=[1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
Ry=[cos(Theta) 0 sin(Theta); 0 1 0; -sin(Theta) 0 cos(Theta)];
Rz=[cos(Psi) -sin(Psi) 0; sin(Psi) cos(Psi) 0; 0 0 1];
R=Rx*Ry*Rz;

Rr = [1 0 -sin(Theta); 0 cos(Phi) sin(Phi)*cos(Theta); 0 -sin(Phi) cos(Phi)*cos(Theta)];     
%----------------------------------------      

F=b*sum(w2);
tau=A*w2;
%----------------------------------------

dPos=Vel;
dAtt=Rr*Ang;
dVel=R*[0;0;F/m]-[0; 0; g];
dAng=J\(cross(-Ang,J*Ang)+tau);

 sys = [dPos; dAtt; dVel; dAng];
end

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u,b)

w2=([u(1);u(2);u(3);u(4)]).^2;
Att=x(4:6);
Phi=Att(1); Theta=Att(2); Psi=Att(3);

Rx=[1 0 0; 0 cos(Phi) -sin(Phi); 0 sin(Phi) cos(Phi)];
Ry=[cos(Theta) 0 sin(Theta); 0 1 0; -sin(Theta) 0 cos(Theta)];
Rz=[cos(Psi) -sin(Psi) 0; sin(Psi) cos(Psi) 0; 0 0 1];
R=Rx*Ry*Rz;

Rr = [1 0 -sin(Theta); 0 cos(Phi) sin(Phi)*cos(Theta); 0 -sin(Phi) cos(Phi)*cos(Theta)]; 

k=0;
for i=1:4
        f=R*[0;0;b*w2(i)];
        for j=1:3
           k=k+1;
           T(k)=f(j);
        end;
end;

F=R*[0;0;b*sum(w2)];
%F=[F(1) F(2) F(3)];

%sys = x;
sys = [x(1:6); R\x(7:9);  Rr*x(10:12); T; F];

end
% end mdlOutputs