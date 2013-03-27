function res = sensn(F,ploton,dispon)
% SENSN simulates transfer through n layers (generalization of sens to n layers)
%
% PHYSICAL FORMULATION
% pi = ki.Ci (Ci in mass/volume)
% p1=p2= => K12 = C1/C2 = k2/k1 at equilibrium
%
% continuity at interface: D1.r1.dC1/dx = D2.r2.dC2/dx
% => D1*.dp1/dx = D2*.dp2/dx with Di* = (ri/ki).Di
%
% transport eq: dCi/dt = Di.d2Ci/dx2
% with: Di* = Di/ki and pi = ki.Ci:
% ri/ki.dpi/dt = Di*.d2pi/dx2
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The current formulation (june 2006) assumes:
%   that all conc are expressed in mass/volume (density already corrected)
%     e.g. p = keq.Ceq with keq = k/r and Ceq = C*r/C0max
%   the reference layer is chosen as the layer with the lowest Di*/li value (mass transport resistance)
%   the equivalent dilution factor respectively to the reference layer must be found with the following approach
%       crit   = inline('sum(r.*(l0+cumsum(l)).^(g+1)-(l0+[0 cumsum(l(1:end-1))]).^(g+1))/l0-L','l0','l','r','g','L');
%       l0     = lref*fzero(crit,100,[],l/lref,r,g,L);
%           where  L is the mass dilution factor between P and F
%                  l0 is the equivalent thickness for the food
%                  g is the geometry number (0=cartesian, 1=cylindrical, 2=spherical)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% MS-MATLAB-WEB 1.0 - 17/01/06 - Olivier Vitrac - rev. 28/06/06

% definitions
global timeout
timeout		= 800; % s
% options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Stats','yes','Initialstep',1e-5,'Maxstep',.05,'Maxorder',5);
options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Stats','yes','Initialstep',1e-8,'Maxstep',.01,'Maxorder',2);
Fdefault	= 	struct(...
				'Bi'		, 	1e3,...	Biot [hm.L1/D]
				'k'			,	[1 1 1 1],...[0.5 3 2],...	ki, i=1 (layer in contact with the liquid)
                'D'         ,   [1e-16 1e-14 1e-14 1e-14],... diffusion coefficient
                'k0'        ,   1,... 0 = liquid
                'l'         ,   [50 20 10 120]*1e-6,...[50 20 10 120]*1e-6,... m
				'L'			,	200/1800,...	dilution factor (respectively to iref)
				'C0'		,	[0 500 500 500],...	initial concentration in each layer
				'options'	,	options...
					); % if iref is mission, it is indentified
Fdefault.t	= 0:.00001:.005; %[0:.00001:.005 .01:.01:.1 .11:.1:5]';
method		= 'cubic';
BiKminode	= 0.1;
MaxStepmax	= 10;
nstepchoice = 200;
n_default	= 1e4;
if Fdefault.Bi<10, Fdefault.t = Fdefault.t/(Fdefault.Bi/10); end
ploton_default = false;
dispon_default = false;
nmesh_default  = 50; % number of nodes for a layer of normalized thickness 1
geometry_default = 0; % slab


% arg check
initon = false;
if ~nargin, initon = true; end
if nargin<1, F = []; end
if nargin<2, ploton = []; end
if nargin<3, dispon = []; end
if isempty(F), F = Fdefault; end
if ~isfield(F,'autotime'), F.autotime = 'on'; end
if ~isfield(F,'n'), F.n = n_default; end
if ~isfield(F,'geometry'), F.geometry = geometry_default; end
if ~isfield(F,'nmesh'), F.nmesh = nmesh_default; end
if ~nargout, ploton=true; dispon=true; end
if isempty(ploton), ploton = ploton_default; end
if isempty(dispon), dispon = dispon_default; end

% physical check
m = Inf;
for prop = {'D' 'k' 'C0'};
    if strcmp(prop{1},'C0')
        F.(prop{1}) = F.(prop{1})(F.(prop{1})>=0);
    else
        F.(prop{1}) = F.(prop{1})(F.(prop{1})>0);
    end
    m = min(m,length(F.(prop{1})));
end
F.m = m;
% init
if initon
    res = Fdefault;
    if dispon, disp('... init mode'), end
    return
end
if strcmp(lower(F.autotime),'on')
	ti		= linspace(min(F.t),max(F.t),F.n)';
else
	ti		= F.t;
end
MaxStepmax = max(MaxStepmax,ceil(F.t(end)/nstepchoice));

% renormalization (fields X1->Xn)
F.k = F.k/F.k0; F.k0 = 1;   % by convention
a = F.D(1:m)./F.k(1:m);
if isfield(F,'iref')
    iref = F.iref;
else
    [crit,iref] = min( a./F.l(1:m) );
    F.iref = iref;
end
F.a = a./a(iref);
F.lref  = F.l(1:m)./F.l(iref);
F.lrefc = cumsum(F.lref);
F.C0    = F.C0(1:m);
F.l     = F.l(1:m);
F.D     = F.D(1:m);

% check if restrictions applied %added 24/01/07
if isfield(F,'ivalid')
    F.ivalid = F.ivalid(F.ivalid<m);
    F.a = F.a(F.ivalid);
    F.lref = F.lref(F.ivalid);
    F.lrefc = F.lrefc(F.ivalid);
    F.C0 = F.C0(F.ivalid);
    F.l = F.l(F.ivalid);
    F.D = F.D(F.ivalid);
    F.k = F.k(F.ivalid);
    m = length(F.ivalid);
    F.m = m;
end
% mesh generation
xmesh = unique([linspace(0,sum(F.lref),F.lrefc(end)*F.nmesh*F.m) F.lrefc])'; %mesh1D(300,[0 1],[.4 .2]); %

% update mesh to geometry
if F.geometry>0 && xmesh(1)==0
    pos0 = F.L;
    F.lrefc = F.lrefc + pos0;
    xmesh = xmesh + pos0;
else
    pos0 = 0;
end
F.pos0 = pos0;
F.pos = [pos0-eps F.lrefc]; % only used for extrapolation

% equilibrium concentration
if m>1
    C0eq = trapz(xmesh,F.C0(ilayer(F,xmesh))'.*xmesh.^F.geometry)/( (1/F.L).^(F.geometry+1) + trapz(xmesh,1./F.k(ilayer(F,xmesh))'.*xmesh.^F.geometry) );
else % added 23/01/07
    C0eq = trapz(xmesh,F.C0(ilayer(F,xmesh)).*xmesh.^F.geometry)/( (1/F.L).^(F.geometry+1) + trapz(xmesh,1./F.k(ilayer(F,xmesh)).*xmesh.^F.geometry) );
end
F.peq = F.k0 * C0eq;


% solving
if dispon, disp('SENSN: solving...'), end
errpde = false;
tic, [sol,Fca,Fcb] = pdesolve1D(F.geometry,@pdediff_PDE,@pdediff_IC,@pdediff_BC,xmesh,F.t,F.options,F);
F.x = xmesh;

% extraction of the solution
if ~errpde
    if dispon, disp(['SENSN: end in ' num2str(toc) ' s']), end
    nt      = length(F.t);
    p		= sol(:,:,1)*F.peq;
    if F.m>1
        C   = p./repmat(F.k(ilayer(F,xmesh)),nt,1);
    else
        C   = p./F.k;
    end
    V       = trapz(xmesh,xmesh.^F.geometry);
    V0      = (1/F.L).^(F.geometry+1);
    Cm		= trapz(F.x,C.*repmat(xmesh',size(C,1),1).^F.geometry,2)/V;
    if m>1
        C0      = trapz(xmesh,F.C0(ilayer(F,xmesh))'.*xmesh.^F.geometry)/V; %     sum(F.C0.*F.lref)/sum(F.lref);
    else % added 23/01/07
        C0      = trapz(xmesh,F.C0(ilayer(F,xmesh)).*xmesh.^F.geometry)/V; %     sum(F.C0.*F.lref)/sum(F.lref);
    end
    C1		= C(:,1);
    fc		= pdeloss(F.x,C,0,F.geometry);
    Cmi		= interp1(F.t,Cm,ti,method);
    fci		= interp1(F.t,fc,ti,method);
    C1i		= interp1(F.t,C1,ti,method);
    p1i		= interp1(F.t,p(:,1),ti,method);
    fi		= F.Bi*p1i-F.Bi*(F.L/F.lrefc(end))*fci; % S.Bi*((Ua-0) - S.L * -Fc(1) )
    if fi(end)<0, fi = fi-min(fi); end
    res.p = p;
    res.peq = F.peq;
    res.V  = V;
    res.C0 = C0;
    res.F  = F;
    res.x  = F.x-pos0;
    res.Cx = C;
    res.tC = F.t;
    res.t  = ti;
    res.C  = Cmi;
    res.CF = (res.C0-res.C)*V/V0;
    res.fc = fci;
    res.f  = fi;
    res.fD = fi/F.Bi;
    
else
    if dispon, disp(['SENSN is TIMEOUT => check your parameters']), end
    res = [];
end

% plots
if ploton
    figure(1)
    subplot(231), hold on, plot(F.x,C), xlabel('x'), ylabel('C')
    subplot(232), hold on, plot(ti,Cmi,'b-'), xlabel('Fo')
    subplot(233), hold on, plot(sqrt(ti),Cmi,'b-'), xlabel('Fo^{1/2}')
    subplot(234), hold on, plot(ti,fi,'b-'), ylabel('j')
    subplot(236), hold on, plot(Cmi,fi,'b-'), xlabel('C'), ylabel('j')
    figure(2)
    subplot(121), hold on, plot(ti/F.Bi,Cmi,'k:'), xlabel('Fo/Bi'), ylabel('C')
    subplot(122), hold on, plot(Cmi,1/F.Bi*fi,'k:'), xlabel('C'), ylabel('j/Bi')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P H Y S I C A L   F O R M U L A T I O N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i = ilayer(S,x)
% ilayer identifies the position of any layer according to the abscissae x
i = min(floor(interp1q(cumsum([0 S.lref])',(0:S.m)',x(:)-S.pos0))+1,S.m); % NEW LINE (24/02/06 OV)

% PDE problem :  [C,F,S] = PDEFUN(X,T,U,DUDX,P1,P2...)
function [c,f,s] = pdediff_PDE(x,t,U,dUdx,S)
global timeout
if toc>timeout, eval('erron'), end
i = find(floor(S.pos-x)<0,1,'last'); %i = find(floor([-eps S.lrefc]-x)<0,1,'last');
c = 1/S.k(i);
f = S.a(i)*dUdx;
s = 0;

% IC : U = ICFUN(X,P1,P2...)
function U0 = pdediff_IC(x,S)
%i = min(floor(interp1q(cumsum([0 S.lref])',(0:S.m)',x-S.lrefc(1)))+1,S.m); % OLD LINE
% i = min(floor(interp1q(cumsum([0 S.lref])',(0:S.m)',x-S.pos0))+1,S.m); % NEW LINE (24/02/06 OV)
i = ilayer(S,x);
U0 = S.C0(i).*S.k(i)/S.peq;
if size(U0,1)==size(x,2), U0 = U0'; end


% BC : [PL,QL,PR,QR] = BCFUN(XL,UL,XR,UR,T,P1,P2...)
function [pa,qa,pb,qb,purge] = pdediff_BC(xa,Ua,xb,Ub,t,Fc,S)
global Hsmooth
% pa = 0;
% qa = 1;
% pb = S.Bi*((Ub-0) - S.L * - Fc(2) );
% qb = 1;
pa = S.Bi*((Ua-0) - S.L * -Fc(1) );   % K=S.k(1)/S.k0 with S.k0=1
qa = -1;
pb = 0;
qb = 1;

% General 1D PDE formulation (without advection)%
%  c(x,t,u,Du/Dx) * Du/Dt = x^(-m) * D(x^m * f(x,t,u,Du/Dx))/Dx + s(x,t,u,Du/Dx)
%  --------------                            --------------       --------------
%   capacitance                                   flux                source
% (diagonal positive nxn matrix)                nx1 vector           nx1 vector
%
%  u(t,x) can be defined in R^n on a<=x<=b during t0 <= t <= tf
%  m = 0, 1, or 2 => {slab, cylindrical, or spherical} symmetry
%  subjected to BC:  p(x,t,u) + q(x,t) * f(x,t,u,Du/Dx) = 0
%                               -------
%                           diagonal nxn matrix
%
