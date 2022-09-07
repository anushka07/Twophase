
%close all
clear all
clear('f_phi')

Zmesh = 5000;

K=0.001;
k1=1e-5;                    % Ps/L scales with k1
%k1=0;
Kalpha=0.00;                   % 1/alpha = 0, No slip condition
Q=1.0;
etaf= 1.0;
%phinr = linspace(0.10,0.49,20);
phinr = 0.1;
tend = zeros(1,length(phinr));
% Initial conditions

z = linspace(0,1,Zmesh)';
h = 1/(Zmesh-1);
phiz = zeros(1,length(z));

A = K/(2*etaf);

b = 0.5;        % Depth of constriction
b = b/2;
a =1-b;d = 100;
% phim = a + b*abs(erf((z-0.5).^3*d));
% phimz = 6*b*d*exp(-(d*(z-0.5).^3).^2).*(z-0.5).^2.*erf((z-0.5).^3*d)./(sqrt(pi).*abs(erf((z-0.5).^3*d)));

phim = a - b*(erf(d*(z-0.5)));
phimz = -(2*b*d*exp(-d^2*(z - 1/2).^2))/pi^(1/2);


I = ones(Zmesh,1);

for J = 1:length(phinr)

phin = phinr(J);

ff = (phin^2/(1-phin)^2);             % Fr t<0; phim is constant at 1
%dfm = (-2*phin^2/(1-phin)^3)*phimz;
G = fsolve(@(g)   2*pi*(K/(2*etaf))*(g) + 2*pi*(g/(ff*16) + Kalpha*g/(4)) - Q, 1);

G = I*G;
phi = I*phin;


%path1 = ['C:\Users\anush\OneDrive - University College London\UCL Phd\Mathematica files\Time dependent solution\'];
% mkdir(path1);
% fid=fopen('MyFile.txt','wt');

% h = figure(1);
% plot(z,phi,'k-',z,phim,'--','LineWidth',2);
% h.Position = [10 40 280 210];

% %Time step
% dt = 5*1e-6;
% h = dt;
% t=0;
% NSTEP = 1000000;

%% Initialize video
%% myVideo = VideoWriter('b_05_phin035_upwindscheme.avi'); %open video file
% myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
% open(myVideo)


%% ode15s Solver - RK4 with variable time step
%options = odeset(odeset('MaxStep',.001));

%Opt    = odeset('Events', @myEvent);%, 'MaxStep',1e-2);
%Opt    = odeset('RelTol', 1e-4);

sol = ode15s(@(t,y) f_phi(t,y,Zmesh,G,K,k1,Q,Kalpha,etaf,b,d),[1e-5 1e2],phi);%,Opt);
t = sol.x;
phi = sol.y;
phi(end,:) = (4*phi(end-1,:) - phi(end-2,:))/3;

%  t = linspace(1e-5,sol.x(end), 10);
tend(J) = sol.x(end);
%  phi = deval(t,sol) ;
    
%% Extract variables - G, phi, Qs

% f = phi.^2./(phim -phi).^2;
% r = pi*(0.125./f + 0.5*Kalpha);      %r = Qs/G
% 
% % f(phi,phim) derivatives wrt phi and phim
% dfm = ((-2*(phi).^2)./(phim - phi).^3).*phimz;
% df = 2*phim.*phi./(phim - phi).^3;
% 
% g = zeros(Zmesh,length(t));
% 
% % Pressure gradient from Qs+Qf = Q
% g(2:end-1,:) = (Q - 2*pi*A*(k1*dfm(2:end-1,:) + k1*df(2:end-1,:).*(phi(3:end,:) - phi(1:end-2,:))/(2*h)))./(2*pi*A + r(2:end-1,:));
% g(end,:) = (Q - 2*pi*A*(k1*dfm(end,:) +k1*df(end,:).*(0)/(2*h)))./(2*pi*A + r(end,:));
% g(1,:) = (Q - 2*pi*A*(k1*dfm(1,:) +k1*df(1,:).*(0)/(2*h)))./(2*pi*A + r(1,:));
% 
% Qs = 2*pi*(g./(f*16) + g*Kalpha/(4));
% Qsg = pi*(0.125./f + 0.5*Kalpha);
% Qsf = -(pi/8)*(g./(f.^2));
 
% S = Qs + phi.*Qsf.*df;
% figure(100);
% hold on
% plot(z,S(:,1:length(t)));

%% Store in txt file
% 
% T = table(z,'VariableNames',{'z'});
% 
% for v = 1:length(t)
% T1 = table(phi(:,v), 'VariableNames', {num2str(v)} );
% T = [T T1];
% end
% 
% % Write data to text file
%    
% whereToStore=fullfile(path1,['phi_vs_t' num2str(k1) '_' num2str(phin) '.txt']);
% writetable(T, whereToStore);
% fclose(fid);

%% Plot variables

    c = figure(J);
    plot(z,phi(:,1:length(t)),'LineWidth',2);
    %plot(z,phi(:,end),'LineWidth',2);
    str1='$\phi $ vs z, $t_{end} =$ ';
%   title([str1,num2str(t),str2, num2str(istep)],'Interpreter','latex');
    title([str1,num2str(t(end))],'Interpreter','latex');
    %hold on
    c.Position = [600 40 580 350];

% for i = 1:length(t)
% n = figure(2);
% n.Position = [600 40 580 350];
% plot(z,phi(:,i),'LineWidth',1.5);
% title(['$t=$',num2str(t(i))],'Interpreter','latex');
% pause(0.1);
% frame = getframe(gcf); %get frame
% L = frame2im(frame);        
% fpath = 'C:\Users\anush\OneDrive - University College London\UCL Phd\Mathematica files\Time dependent solution\ode45 method\';
% %whereToStore=fullfile(fpath,['phi_z' num2str(b) '.png']);
% %imwrite(L, whereToStore);
% writeVideo(myVideo,L);
% end
% close(myVideo);           

    
    %savefig(c,[path1 'b' num2str(b) '_d' num2str(d) '_' 'phi_vs_z.fig']);

%     o = figure(2);
%     plot(z,Qs(:,1:length(t)),'LineWidth',2);
%     %plot(z,Qs(:,end),'LineWidth',2);
%     str1='$Qs $ vs z at $t = $';
%     str2 = '$ NStep =$';
% %   title([str1,num2str(t),str2, num2str(istep)],'Interpreter','latex');
%     title([str1,num2str(t(end))],'Interpreter','latex');   
%     hold on
%     o.Position = [40 40 250 250];
%     
end

% figure(1000);
% plot(phinr, tend,'LineWidth',2);


function dydt= f_phi(t,y,Zmesh,G,K,k1,Q,Kalpha,etaf,b,d)

z = linspace(0,1,Zmesh)';
h = 1/(Zmesh-1);

phiz = zeros(Zmesh,1);
Gz = zeros(Zmesh,1);
dydt = zeros(Zmesh,1);
A = K/(2*etaf);

a =1-b;%d = 100;
% phim = a + b*abs(erf((z-0.5).^3*d));
% phimz = 6*b*d*exp(-(d*(z-0.5).^3).^2).*(z-0.5).^2.*erf((z-0.5).^3*d)./(sqrt(pi).*abs(erf((z-0.5).^3*d)));

phim = a - b*(erf(d*(z-0.5)));
phimz = -(2*b*d*exp(-d^2*(z - 1/2).^2))/pi^(1/2);


%Time derivatives of phi at discrete z points
y(end) = (4*y(end-1) - y(end-2))/3;

f = y.^2./(phim -y).^2;
r = pi*(0.125./f + 0.5*Kalpha);      % r = Qs/G

% f(phi,phim) derivatives wrt phi and phim
dfm = ((-2*(y).^2)./(phim - y).^3).*phimz;
df = 2*phim.*y./(phim - y).^3;

% Pressure gradient from Qs+Qf = Q
%phiz(2:end-1) = (y(2:end-1) - y(1:end-2))/h;
phiz(2:end-1) = (y(3:end) - y(2:end-1))/h;
phiz(1) = 0;

%G(2:end-1) = (Q - 2*pi*A*(k1*dfm(2:end-1) + k1*df(2:end-1).*(y(3:end) - y(1:end-2))/(2*h)))./(2*pi*A + r(2:end-1));
G(2:end-1) = (Q - 2*pi*A*(k1*dfm(2:end-1) + k1*df(2:end-1).*(phiz(2:end-1))))./(2*pi*A + r(2:end-1));

Qs = 2*pi*(G./(f*16) + G*Kalpha/(4));
Qsg = pi*(0.125./f + 0.5*Kalpha);
Qsf = -(pi/8)*(G./(f.^2));

%phiz(2:end-1) =  -(A*2*pi*(G(2:end-1) + k1*dfm(2:end-1)) + Qs(2:end-1) - Q)./(k1*2*pi*A*df(2:end-1)) ;   
% phiz(2:end-1) =  (y(3:end) - y(1:end-2))/(2*h);
% phiz(1) = (y(2) - y(1))/(2*h);



% Using phiz(1) = phiz(end) = 0
G(end) = (Q - 2*pi*A*(k1*dfm(end) +k1*df(end).*(0)/(h)))./(2*pi*A + r(end));
G(1) = (Q - 2*pi*A*(k1*dfm(1) + k1*df(1).*phiz(1)))./(2*pi*A + r(1));


%Gz(2:end-1) = (G(3:end) - G(1:end-2))/(2*h);

%Gz(2:end-1) = (G(2:end-1) - G(1:end-2))/(h);
Gz(2:end-1) = (G(3:end) - G(2:end-1))/(h);

%Central difference
dydt(2:end-1) = -(y(2:end-1).*Qsg(2:end-1).*Gz(2:end-1) + y(2:end-1).*Qsf(2:end-1).*dfm(2:end-1) + (Qs(2:end-1) + y(2:end-1).*Qsf(2:end-1).*df(2:end-1)).*(phiz(2:end-1)) ) ;


%To remove out of equation
dydt(end)=0;

dydt(1) = 0;
dydt(2) = 0;


% for i=1:Zmesh
%     f = ((y(i))^2/(phim(i)-y(i))^2);
%     r = pi*(0.125/f + 0.5*Kalpha);      % r = Qs/G
%     
%     % f(phi,phim) derivatives wrt phi and phim
%     dfm = ((-2*(y(i))^2)/(phim(i) - y(i))^3)*phimz(i);
%     df = 2*phim(i)*y(i)/(phim(i) - y(i))^3;
%     
%     for j=2:Zmesh-1
%         % Pressure gradient from Qs+Qf = Q
%     
%     G(j) = (Q - 2*pi*A*(k1*dfm +k1*df*(y(j+1) - y(j-1))/(2*h)))/(2*pi*A + r);
%     end
%     
%     if(i>1&&i<Zmesh)
%     Gz = (G(i+1) - G(i-1))/(2*h);
%     end
%     
%     % Newtonian Fluid
%     Qs = 2*pi*(G(i)/(f*16) + G(i)*Kalpha/(4));
%     Qsg = pi*(0.125/f + 0.5*Kalpha);
%     Qsf = -(pi/8)*(G(i)/((f^2)));
%     
%     
%     if(i>1&&i<Zmesh)
%     dydt(i) = -(y(i)*Qsg*Gz + (Qs + y(i)*Qsf*df )*(y(i+1) - y(i-1))/(2*h) + y(i)*Qsf*dfm) ;
%     else
%     dydt(i) = -(y(i)*Qsf*dfm) ;%phi_z, g_z = 0 at boundaries
%     end
% end

% b = figure(1);
%     plot(z,y,'LineWidth',2);
%     str1='$\phi $ vs z at $t = $';
%     title([str1,num2str(t)],'Interpreter','latex');
%     hold on
%     b.Position = [600 40 580 350];

end

function [value, isterminal, direction] = myEvent(T, Y)
if(sum(Y(2:end) - Y(1:end-1)>1e-3)>0)
    value=1;
else
    value=0;
end
%value      = (sum(Y(:) <= 0));
isterminal = 1;   % Stop the integration
direction  = 0;
end

