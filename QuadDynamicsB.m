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
    sys=mdlDerivatives(x,t,u,quad);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(x,quad);

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
function sys=mdlDerivatives(x,t,u,quad)

Pos=x(1:3);
Att=x(4:6);
Vel=x(7:9);
Ang=x(10:12);

phi=Att(1);
the=Att(2);
psi=Att(3);


R = [cos(psi)*cos(the)  cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi)  cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);
     sin(psi)*cos(the)  sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi)  sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi);
     -sin(the)  cos(the)*sin(phi)  cos(the)*cos(phi)];

Rr_inv =[1  sin(phi)*tan(the)  cos(phi)*tan(the);
     0  cos(phi)  -sin(phi);
     0  sin(phi)/cos(the) cos(phi)/cos(the)]; 

w2=u.^2;
 
%for i=1:4
%w2(i)=abs(u(i))*u(i);
%end; 

Wg=u(1)-u(2)+u(3)-u(4);

%----------------------------------------                                
F=quad.b*sum(w2);
tau=[quad.l*quad.b*(-w2(2)+w2(4));
     quad.l*quad.b*(-w2(1)+w2(3));
     quad.k*(w2(1)-w2(2)+w2(3)-w2(4))];

%----------------------------------------
dPos=Vel; 
dAtt=Rr_inv*Ang;
dVel=R/quad.m*[0;0;F]-[0;0;quad.g]-0*1/quad.m*quad.A*Vel; 
%Drag force in the corresponding directions: -1/quad.m*quad.A*Vel
dAng=quad.J\(-cross(Ang,quad.J*Ang)-0*quad.Jr*cross(Ang,[0;0;1])*Wg+tau);  
%Gyroscopic forces: -quad.Jr*cross(Ang,[0;0;1])*Wg

sys = [dPos; dAtt; dVel; dAng];
end

% end mdlDerivatives
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(x,quad)

if x(3)<0 
    x(3)=0;
 %   if x(6) > 0
  %      x(6) = 0;
  % end
end;

Att=x(4:6);
phi=Att(1);
the=Att(2);
psi=Att(3);

R = [cos(psi)*cos(the)  cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi)  cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);
     sin(psi)*cos(the)  sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi)  sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi);
     -sin(the)  cos(the)*sin(phi)  cos(the)*cos(phi)];

Rr =[1  0  -sin(the);
     0  cos(phi)  cos(the)*sin(phi);
     0  -sin(phi) cos(the)*cos(phi)];

sys = [x(1:6); R\x(7:9);  Rr*x(10:12)];

end
% end mdlOutputs