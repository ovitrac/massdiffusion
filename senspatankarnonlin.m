function res = senspatankarnonlin(F,ploton,dispon)
%SENSPATANKAR simulates transfer through n layers using a modified Patankar Method (see p 45)

% MS-MATLAB-WEB 1.0 - 07/04/08 - Olivier Vitrac - rev.

% Revision history
% 01/10/07 improve speed

% definitions
global timeout
timeout		= 800; % s
% options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Stats','yes','Initialstep',1e-5,'Maxstep',.05,'Maxorder',5);
options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.1,'Maxorder',2,'bdf','on');
D           = 1e-12*[1 10 6]; %[1e-11 1e-11 1e-11]; %;
P           = [1 0.3 0.25 10]; %[1 1 1 1]; %
xD          = [0.3690 0.6310];
xP          = [0.3690 0.6310];
Fdefault	= 	struct(...
				'Bi'		, 	1e3,...	Biot [hm.L1/D]
                'D'         ,   Dmodel([0 xD 1]*2.7,D),... diffusion coefficient
                'I'         ,   Imodel([0 xP 1],P),... sorption isotherm
                'ke'        ,   1,...10
                'pe'        ,   .6,... external concentration
                'l'         ,   100e-6,...
				'L'			,	100000,...	dilution factor (respectively to iref)
				'C0'		,	0,...	initial concentration
				'options'	,	options...
					);
Fdefault.t	= [0:.00001:.005 .01:.01:.1 .11:.1:50]; %0:.00001:.005; %[0:.00001:.005 .01:.01:.1 .11:.1:5]';
if Fdefault.Bi<10, Fdefault.t = Fdefault.t/(Fdefault.Bi/10); end
method		= 'cubic'; %'cubic';
ploton_default = false;
dispon_default = false;
nmesh_default  = 100; % number of nodes for a layer of normalized thickness 1
n_default	= 1e4;
MaxStepmax	= 10;
nstepchoice = 200;
zero = 1000*eps;

% arg check
initon = false;
if ~nargin, initon = true; end
if nargin<1, F = []; end
if nargin<2, ploton = []; end
if nargin<3, dispon = []; end
if isempty(F), F = Fdefault; end
if ~isfield(F,'autotime'), F.autotime = 'on'; end
if ~isfield(F,'n'), F.n = n_default; end
if ~isfield(F,'nmesh'), F.nmesh = nmesh_default; end
if ~isfield(F,'options'), F.options = options; end
if ~nargout, ploton=true; dispon=true; end
if isempty(ploton), ploton = ploton_default; end
if isempty(dispon), dispon = dispon_default; end

% physical check
if initon
    res = Fdefault;
    if dispon, disp('... init mode'), end
    return
end

if strcmpi(F.autotime,'on')
	ti		= linspace(min(F.t),max(F.t),F.n)';
else
	ti		= F.t;
end
MaxStepmax = max(MaxStepmax,ceil(F.t(end)/nstepchoice));

% normalization
Dref  = fnmin(F.D); % min D value
F.k_1   = fnder(F.I.C,1); % solubility or 1/k
F.k   = fnder(F.I.p,1); % solubility or 1/k
kref_1  = fnmin(F.k_1); % max slope
kref  = -fnmin(fncmb(F.k,-1)); % max
D     = fncmb(F.D,1/Dref);
%k_1   = fncmb(F.k_1,1/kref_1);
k     = fncmb(F.k,1/kref_1);
k0    = F.ke/kref; %k0    = F.ke*kref_1;
% a     = @(C) fnval(D,C).*fnval(k_1,C);
a     = @(C) fnval(D,C)./fnval(k,C);
a0    = F.Bi/k0;

% init
global A
C0       = zeros(F.nmesh+1,1);
C0(1)    = F.pe/F.ke;
A        = zeros(F.nmesh+1,F.nmesh+1);
de       = 1/(2*F.nmesh);
dw       = de;
vol_1    = 1/(dw+de);
xmesh    = linspace(dw,1-de,F.nmesh);
R        = @(C) [1/a0;de./a(C(2:end))]; % normalized resistances
I        = @(C) [F.ke*C(1);fnval(F.I.p,C(2:end))];

% plots (debug purposes)
clf
hs = subplots([.6 1],1,.1); p0 = get(hs(2),'position'); delete(hs)
hs = subplots([.6 1],[.8 .8 .4],.1,.1,'alive',1:3);
hs = [hs;subplots(1,[1 1],.1,.1,'position',p0)];
pi = linspace(0,1,1000);
Ci = fnval(F.I.C,pi);

subplot(hs(1)), plot(pi,Ci,'linewidth',2); xlabel('p'), ylabel('C')
subplot(hs(2)), plot(Ci,fnval(F.D,Ci),'linewidth',2); xlabel('C'), ylabel('D'), set(hs(2),'yscale','log')
subplot(hs(4)), hplot{1} = plot([C0 I(C0)],'-','linewidth',2);
legend(hplot{1},{'C' 'p'}); xlabel('position')
subplot(hs(5)), hplot{2} = plot(C0,'linewidth',2); xlabel('position'), ylabel('dC/dt')
subplot(hs(3)), hplot{3} = plot(R(C0)); xlabel('position'), ylabel('R')


% integration
% J = @(t,C)(A);
% F.options.Jacobian = J;
F.options.Vectorized = 'on';
[t,C] = ode15s(@dCdt,F.t,C0,F.options,R,I,F.L,vol_1,hplot,hs(3)); % integration
CF = C(:,1);
C  = C(:,2:end);

% % interpolation at each interface
% Ce = zeros(length(t),F.nmesh);
% Cw = zeros(length(t),F.nmesh);
% for i=1:length(t)
%     Ce(i,1:end-1) = C(i,1:end-1) - (de(1:end-1).*he(1:end-1).*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(1:end-1) )';
%     Ce(i,end)     = C(i,end-1);
%     Cw(i,2:end)   = C(i,2:end)   + (dw(2:end)  .*hw(2:end)  .*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(2:end)   )';
%     Cw(i,1)       = C(i,1)       +  dw(1)       *hw(1)       *( F.k0       /k(1)      * CF(i)         - C(i,1)      )  / D(1);
% end
% Cfull = reshape([Cw;C;Ce],length(t),3*F.nmesh);
% xw = xmesh-dw+zero;
% xe = xmesh+de-zero;
% xfull = [xw';xmesh';xe']; xfull = xfull(:);
% kfull = [k';k';k']; kfull = kfull(:);
% 
% % outputs
% res.C = interp1(t,trapz(xfull,Cfull,2)*C0eq,ti,method)/xfull(end); % av conc
% res.t = ti; % time base
% res.p = NaN; %interp1(t,repmat(kfull',length(t),1).*Cfull*C0eq,ti,method); % to fast calculations
% res.peq = F.peq;
% res.V  = F.lrefc(end); %*F.l(F.iref);
% res.C0 = C0*C0eq;
% res.F  = F;
% res.x  = xfull; %*F.l(F.iref);
% res.Cx = Cfull*C0eq; %interp1(t,Cfull*C0eq,ti,method);
% res.tC = t;
% res.CF = interp1(t,CF*C0eq,ti,method);
% res.fc = res.CF/F.L;
% res.f  = interp1(t,hw(1) * ( F.k0/k(1) * CF - C(:,1) ) * C0eq,ti,method);
% res.timebase = res.F.l(res.F.iref)^2/(res.F.D(res.F.iref));
% plot(res.x,res.Cx')

function r=dCdt(t,C,R,I,L,vol_1,hplot,hcontrol)
global A
persistent lasttime niter
if isempty(lasttime) || (t==0), lasttime=0; niter=0; end
niter = niter + 1;
[m,nC] = size(C);
r = zeros(m,nC);
for j=nC
    C(C<0) = 0;
    Ri = R(C(:,j));
    hw = 1./(Ri(1:end-1)+Ri(2:end));
    he = [hw(2:end);0];
    n = length(hw); %F.nmesh
    A(1,1:2) = 1/L*hw(1) * [-1 1];
    for i=1:n-1
        A(i+1,i:i+2) = vol_1 * [ hw(i), -hw(i)-he(i), he(i) ];
    end
    A(end,end-1:end) = vol_1 * [ hw(n) -hw(n)];
    r(:,j) = A*I(C(:,j));
end
if nC==1 && (t-lasttime)>1e-3
    dispf('%g',t)
    set(hplot{1}(1),'Ydata',C(:,1));
    set(hplot{1}(2),'Ydata',I(C(:,1)));
    set(hplot{2},'Ydata',r(:,1));
    set(hplot{3},'Ydata',Ri);
    drawnow
    lasttime =t;
    
end
