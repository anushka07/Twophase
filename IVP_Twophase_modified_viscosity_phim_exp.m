clearvars
%close all

% This is a code to solve for evolution of phi along a pipe for a Newtonian solid phase and a Darcy fluid phase.
% The critical particle fraction 'phim' controls the flow of solid phase. 
% phim is considered to be a linearly decreasing function phim = 1+c*z, where c is -ve


%% Define Parameters

Kr=0.001;
%ar=linspace(-4,0,500);
%Kr = 2.^ar;

%k1r=linspace(0.01,.1,10);                       %Ps/L scales with k1
%k1r = [1.0 0.75 0.5 0.25];
k1r = 1e-5;
Kalpha=0.00;                   %1/alpha = 0, No slip condition
Qar=[1.0];%0.35 0.4 0.45];
etaf= 1.0;

%phin= 0.7;
phinr = linspace(.1, .45, 25);
%phinr = [0.2];


b = 0.5;
%b =b/2;
d = 1000;


Varar = phinr;
%Varar = cr;
%Varar = ar;
%Varar = k1r;
Es = zeros(length(Varar),1);
Ps = zeros(length(Varar),1);
Pf = zeros(length(Varar),1);

%% Initial guess 

Q = Qar(1);


for j=1:length(Varar)
   
    Q = Qar(1);phin=phinr(j);K=Kr(1);k1 = k1r(1);

%% ode15s solver
Opt    = odeset('RelTol',1e-4,'MaxStep',1e-5);
I = bcfun(K,k1,Q,Kalpha,etaf,phin,b,d);
[sol] = ode23s(@(z,y)bvpfunc(z,y,K,k1,Q,Kalpha,etaf,b,d),[0 1],I,Opt);
%solinit=sol;

%% Finding Qf, Qs, phim, phi G, Psz
A = K/(2*etaf);
%zmesh_new = z;

% G = y(:,2);                 %-dEs/dz
% phi = y(:,1);
zmesh_new = linspace(0,sol.x(end),1000);

G = deval(sol,zmesh_new,2);                 %-dEs/dz
phi = deval(sol,zmesh_new,1);
phi_z = gradient(phi, zmesh_new);
%phi_z = gradient(y(:,1), z);

%phim = 1 + c*zmesh_new;

% c = 0.5;
% m=5;
% b = -log(1-c)/m;
% phim = c + exp(-m.*(zmesh_new+b));
% phimz = -m*exp(-m.*(zmesh_new+b));

a = 1-b;
phim = a + b*abs(erf((zmesh_new-0.5).^3*d));
phimz = 6*b*d*exp(-(d*(zmesh_new-0.5).^3).^2).*(zmesh_new-0.5).^2.*erf((zmesh_new-0.5).^3*d)./(sqrt(pi).*abs(erf((zmesh_new-0.5).^3*d)));
% phim = a - b*(erf(d*(zmesh_new-0.5)));
% phimz = -(2*b*d*exp(-d^2*(zmesh_new - 1/2).^2))/pi^(1/2);

f = ((phi).^2./(phim-phi).^2);

Qs = (2*pi*G.*(1./(f*16) + Kalpha/(4)));

dfm = ((-2*(phi).^2)./(phim - phi).^3).*phimz;
df = 2*phim.*phi./(phim - phi).^3;

Qf = (A*2*pi*(G + k1*dfm + k1*df.*phi_z)./(1-phi) + Qs);
Psz = k1*dfm + k1*df.*phi_z;                % dPs/dz
Pfz = Psz + G;                              % -dPf/dz

Es(j) = trapz(zmesh_new,G);                    % DEs
Ps (j)= trapz(zmesh_new, Psz);
Pf(j) = trapz(zmesh_new, Pfz);

%% Plot solution
%if(mod(j,10)==0)

g=figure(12);
plot(zmesh_new,phi,'LineWidth',1);%zmesh_new,phim,':',
str='$\phi $ vs z';
title([str],'Interpreter','latex');
g.Position = [10 40 280 210];
xlabel('z','Interpreter', 'latex') ;
str='$\phi$';
ylabel(str,'Interpreter','latex') ;
hold on

% figure(1);
% clf;
% plot(zmesh_new,phi,'k-',zmesh_new,phim,'r--')
% ylabel('\phi and \phi_m (red dashed)')
% xlabel('z')

% g=figure(20);
% plot(zmesh_new,Qs.*phi,'-','LineWidth',1);
% str='$Qs*\phi$';
% title(str,'Interpreter','latex');
% g.Position = [600 40 280 210];
% xlabel('z','Interpreter', 'latex') ;
% str='$Qs*\phi$';
% ylabel(str,'Interpreter','latex') ;
% hold on

g=figure(20);
plot(zmesh_new,G,'-','LineWidth',1);
str='$-\Sigma_z $ vs z';
title(str,'Interpreter','latex');
g.Position = [600 40 280 210];
xlabel('z','Interpreter', 'latex') ;
str='$-\Sigma_z$';
ylabel(str,'Interpreter','latex') ;
hold on

g=figure(21);
plot(zmesh_new,Psz,'LineWidth',1)
str='$P_{sz} $ vs z';
title(str,'Interpreter','latex');
xlabel('z','Interpreter', 'latex') ;
str='$P_{sz}$';
ylabel(str,'Interpreter','latex') ;
g.Position = [900 40 280 210];
hold on

g=figure(22);
plot(zmesh_new,Pfz,'LineWidth',1)
str='$-P_{fz} $ vs z';
title(str,'Interpreter','latex');
xlabel('z','Interpreter', 'latex') ;
str='$-P_{fz}$';
ylabel(str,'Interpreter','latex') ;
g.Position = [10 340 280 210];
hold on


g = figure(23);
plot(zmesh_new, Qs,'LineWidth',1);%zmesh_new, Qf,'-',
str='$Q_{s}$ vs z';
title(str,'Interpreter','latex');
xlabel('z','Interpreter', 'latex') ;
str='$Q_{s}$';
ylabel(str,'Interpreter','latex') ;
g.Position = [300 340 280 210];
hold on

g=figure(24);
plot(zmesh_new, Qf,'LineWidth',1);%zmesh_new, Qf,'-',
str='$Q_{f}$ vs z';
title(str,'Interpreter','latex');
xlabel('z','Interpreter', 'latex') ;
str='$Q_{f}$';
ylabel(str,'Interpreter','latex') ;
g.Position = [600 340 280 210];
hold on

%end

end

%% BVP function definition
function dydx=bvpfunc(z,y,K,k1,Q,Kalpha,etaf,b,d)


A = K/(2*etaf);

% c = 0.5;
% m=5;
% b = -log(1-c)/m;
% phim = c + exp(-m.*(z+b));
% phimz = -m*exp(-m.*(z+b));

a = 1-b;  
phim = a + b*abs(erf((z-0.5)^3*d));
phimz = 6*b*d*exp(-(d*(z-0.5)^3)^2)*(z-0.5)^2*erf((z-0.5)^3*d)/(sqrt(pi)*abs(erf((z-0.5)^3*d)));
% phim = a - b*(erf(d*(z-0.5)));
% phimz = -(2*b*d*exp(-d^2*(z - 1/2).^2))/pi^(1/2);

%phim = 1 + c*z;                                % phi_m variation in z 
f = ((y(1))^2/(phim-y(1))^2);

% Bingham Fluid
%Qs = 2*pi*( (y(2)/4)*(1-2/y(2))^2*(1/2 - ((1-2/y(2))^2)/4 -(2/y(2))*(1-2/y(2))/3) + y(2)/(4*alpha) );
%Qsg = pi*(0.125 + 0.5/alpha - 2*(y(2))^(-4));

% Newtonian Fluid
Qs = 2*pi*(y(2)/(f*16) + y(2)*Kalpha/(4));
Qsg = pi*(0.125/f + 0.5*Kalpha);
Qsf = -(pi/8)*(y(2)/((f^2)));


dfm = (-2*(y(1))^2/(phim-y(1))^3)*phimz;
df  = 2*phim*y(1)/(phim-y(1))^3;

phiz = -(A*2*pi*(y(2) + k1*dfm) + Qs - Q)/(k1*2*pi*A*df);


%dydx=zeros(2,1);                            % phi_z and G_z

dydx =[ 
        phiz
        %(Qf + Qs- Q)/(k1 * A * df/dphi)
        
        
        %(Qs/(y(1)*Qsg))*((A*2*pi*(y(2) + k1*dfm) + Qs  - Q)/(k1*2*pi*A*df)) %dphi/dz
        ((Qs + y(1)*Qsf*df )*(-phiz) - y(1)*Qsf*dfm)/(y(1)*Qsg) 
      
      ];

% figure(1)
% plot(z,dydx(1),'Marker','*');
% hold on


end

%% Boundary condition function

function res = bcfun(K,k1,Q,Kalpha,etaf,phin,b,d)
%syms g;

% A = K/(2*etaf);
% c = 0.5;
% m = 5;
% b = -log(1-c)/m;
% phim = c + exp(-m.*(1+b));
% phimz = -m*exp(-m.*(1+b));

%For dphi/dz = 0 at z=0  
%phim = 1;
%dfm = (-2*(phin)^2/(phim-phin)^3)*c;
%df = 2*phim*phin/(phim-phin)^3; 
%dfm = (-2*(phin)^2/(phim-phin)^3)*phimz;

a = 1-b;
phim = a + b*abs(erf((-0.5)^3*d));
phimz = 6*b*d*exp(-(d*(-0.5)^3)^2)*(-0.5)^3*erf((-0.5)^2*d)/(sqrt(pi)*abs(erf((-0.5)^3*d)));

% phim = a - b*(erf(d*(-0.5)));
% phimz = -(2*b*d*exp(-d^2*( - 1/2).^2))/pi^(1/2);

f = ((phin)^2/(phim-phin)^2);
df = 2*phim*phin/(phim-phin)^3; 
dfm = (-2*(phin)^2/(phim-phin)^3)*phimz;


% Bingham Fluid
%Qs = 2*pi*( (yb(2)/4)*(1-2/yb(2))^2*(1/2 - ((1-2/yb(2))^2)/4 -(2/yb(2))*(1-2/yb(2))/3) + yb(2)/(4*alpha) );
%Qsg = pi*(0.125 + 0.5/alpha - 2*(yb(2))^(-4));

% Newtonian Fluid
% Qs = 2*pi*(yb(2)/16 + yb(2)*Kalpha/(4));
% Qsg = pi*(0.125 + 0.5*Kalpha);

%Newtonian fluid for at z=0
%Qs = 2*pi*(g/16 + g*Kalpha/(4));



% dfm = (-2*(yb(1))^2/(phim-yb(1))^3)*c;    % f = phi^2/(phim-phi)^2 and phimz = c
% df = 2*phim*yb(1)/(phim-yb(1))^3;


% Qs = Qf (no relative velocity)
% F = fsolve(@(g) g + k1*dfm - k1*df*(A*2*pi*(g + k1*dfm) +  2*pi*(g/(f*16) + g*Kalpha/(4) - Q))/(k1*2*pi*A*df)-1e-2, 0 );

% Rough initial guess from bvp
% s = fsolve(@(g)   2*pi*(K/(2*etaf))*g + 2*pi*(g/(f*16) + Kalpha*g/(4)) - Q, 1);

% phi_z = 0
 C = fsolve(@(g)   2*pi*(K/(2*etaf))*(g + k1*dfm) + 2*pi*(g/(f*16) + Kalpha*g/(4)) - Q, 1);
%C = -0.00001;
P = phin;

res = [P ; C];



end