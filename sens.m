function [resode,respde] = sens(mode,F,ploton,dispon)
% SENS simule un transfert RD+RK+RH lors d'un contact solide liquide
%   Syntax: [resode,respde] = sens(mode,F[,ploton,dispon])

% SENS 1.0 - 12/07/02 - Olivier Vitrac - rev. 19/01/06

% revision history
% 14/01/06 add dispon
% 18/01/05 add example for left and right BC (see pdediff_BC)
% 19/01/06 add field geometry

% TO BE MODIFIED ACCORDING TO THE PROCESSOR SPEED
global timeout_ode timeout_pde
timeout_ode		= 1; % s
timeout_pde		= 600; % s

% Definitions
xdefault	= mesh1D(50,[0 1],[.4 .2]);
options		= odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','yes','Initialstep',1e-8,'Maxstep',.1,'Maxorder',5);
Fdefault	= 	struct(...
				'x'			, 	xdefault,...
				'Bi'		, 	1e3,...	Biot [hm.L1/D]
				'K'			,	1,...	Partition coefficient [CL/CP]
				'L'			,	1/9,...	dimensionless thickness [L1/L2]
				'C0'		,	500,...	dimensionless initial concentration []
				'xcrit'		,	0,...	dimensionless critical front position []
				'options'	,	options...
					);
Fdefault.t	= [0:.00001:.005 .01:.01:.1 .11:.1:5]';
method		= 'cubic';
Bi_critode	= 7e4;
Bi_maxode	= 1e18;
BiKminode	= 0.1;
MaxStepmax	= 10;
nstepchoice = 200;
n_default	= 1e4;
if Fdefault.Bi<10, Fdefault.t = Fdefault.t/(Fdefault.Bi/10); end
ploton_default = false;
dispon_default = false;
geometry_default = 0; % slab


% Input control
if nargin<1 & ~nargout, mode = 'auto'; elseif nargin<1,	mode = 'init'; end
if nargin<2, F = Fdefault; end
if nargin<3, ploton = []; end
if nargin<4, dispon = []; end
if ~nargout, mode = 'auto', end
if ~isfield(F,'autotime'), F.autotime = 'on'; end
if ~isfield(F,'n'), F.n = n_default; end
if ~isfield(F,'geometry'), F.geometry = geometry_default; end
if F.geometry>0 && (strcmp(mode,'ode') || strcmp(mode,'both'))
    if dispon, disp('... switch to PDE mode for non cartesian geometries'), end
    mode = 'pde';
end
if isempty(ploton), ploton_default; end
if isempty(dispon), dispon_default; end
switch lower(mode)
case 'auto',	ploton = true; mode = 'both'; if dispon, disp('... auto mode'), end
case 'init',	resode = Fdefault; if dispon, disp('... init mode'), end, return
case 'ode',		if dispon, disp('... ode mode'), end
case 'pde',		if dispon, disp('... pde mode'), end
case 'both',	if nargout>1, if dispon, disp('... pde+ode mode'), end, else, mode = 'ode'; if dispon, disp('... ode mode'), end, end
otherwise,		error('... unkown command')
end
if strcmp(lower(F.autotime),'on')
	ti		= linspace(min(F.t),max(F.t),F.n)';
else
	ti		= F.t;
end
MaxStepmax = max(MaxStepmax,ceil(F.t(end)/nstepchoice));

% update mesh to geometry
if F.geometry>0 && F.x(1)==0
    F.x = F.x + F.x(end)/F.L;
end

% Solution
switch lower(mode), case {'both' 'pde'}
	if dispon, disp('SENSpde: solving...'), end
    %  	if F.Bi*F.K>50, F.options.MaxOrder= 2; end
    errpde = 0;
    try
        tic, [sol,Fca,Fcb] = pdesolve1D(F.geometry,@pdediff_PDE,@pdediff_IC,@pdediff_BC,F.x,F.t,F.options,F);
    catch
        if findstr(lasterr,'erron')
            if dispon, disp('... TIME OUT'), end
            timeout_pdeold = timeout_pde; timeout_pde = 15;
            F.options = odeset(F.options,'bdf','on','maxorder',2,'maxstep',1e5);
            try
                tic, [sol,Fca,Fcb] = pdesolve1D(F.geometry,@pdediff_PDE,@pdediff_IC,@pdediff_BC,F.x,F.t,F.options,F);
            catch
                if findstr(lasterr,'erron'), if dispon, disp('... TIME OUT (2)'), end, else, if dispon, disp(sprintf('...unexpected error ''%s'' for F =',lasterr)), disp(F), end; end
                errpde = 1;
            end
            timeout_pde = timeout_pdeold;
        else
            if dispon, disp('... decrease maxstep'), end
            timeout_pdeold = timeout_pde; timeout_pde = 30;
            F.options.MaxStep = F.options.MaxStep/100;
            try
                tic, [sol,Fca,Fcb] = pdesolve1D(F.geometry,@pdediff_PDE,@pdediff_IC,@pdediff_BC,F.x,F.t,F.options,F);
            catch
                if findstr(lasterr,'erron'), if dispon, disp('... TIME OUT'), end, else, if dispon, disp(sprintf('...unexpected error ''%s'' for F =',lasterr)), disp(F); end, end
                errpde = 1;
            end
            timeout_pde = timeout_pdeold;        
        end
    end
    if ~errpde
        if dispon, disp(['SENSpde: end in ' num2str(toc) ' s']), end
        % Extraction
        C		= sol(:,:,1);
        C1		= C(:,end);
        F.x     = F.x(:)';
        if F.geometry>0
            Cm		= trapz(F.x,C.*repmat(F.x,size(C,1),1).^F.geometry,2) / trapz(F.x,F.x.^F.geometry);
            fc		= pdeloss(F.x,C,0,F.geometry);
        else
            Cm		= trapz(F.x,C,2);
            fc		= pdeloss(F.x,C);
        end
        Cmi		= interp1(F.t,Cm,ti,method);
        fci		= interp1(F.t,fc,ti,method);
        C1i		= interp1(F.t,C1,ti,method);
        fi		= F.Bi*F.K*C1i-F.Bi*F.L*fci;
        if fi(end)<0, fi = fi-min(fi); end
        % 	fi		= ndf(ti,fci);
        % 	fi(1)	= F.K*F.Bi*F.C0;
    else
        if dispon, disp(['SENSpde is TIMEOUT => switch to SENSode']), end
    end
end
switch lower(mode), case {'both' 'ode'}
    tic, if dispon, disp('SENSode: solving...'), end
    nsteptypical = ceil(F.t(end)/F.options.MaxStep);
    if (nsteptypical>nstepchoice) | (F.Bi*F.K<BiKminode)
        if dispon, disp(sprintf('... increase maxstep up to %0.2g',MaxStepmax)), end
        factor = ceil(MaxStepmax/F.options.MaxStep);
        F.options = odeset(F.options,'maxstep',MaxStepmax,'Initialstep',factor*F.options.InitialStep);
    elseif F.Bi*F.K<10*BiKminode
        step = MaxStepmax*BiKminode/(F.Bi*F.K);
        if dispon, disp(sprintf('... increase maxstep up to %0.2g',step)), end
        factor = ceil(step/F.options.MaxStep);
        F.options = odeset(F.options,'maxstep',step,'Initialstep',factor*F.options.InitialStep);
    end
    if F.Bi<Bi_maxode
        try
            if F.Bi>=Bi_critode
                F.options = odeset(F.options,'bdf','on','maxorder',2,'maxstep',.1);
                if dispon, disp('... bdf is ''on'''), end
            end
            tic
            [tC,Code]	= ode15s(@odediff,F.t,[F.C0],F.options,F);
        catch
            if findstr(lasterr,'erron')
                if dispon, disp('... TIME OUT'), end
                F.options = odeset(F.options,'maxstep',1e5,'Initialstep',.1);
                timeout_odeold = timeout_ode; timeout_ode = 10;
                tic
                try
                    [tC,Code]	= ode15s(@odediff,F.t,[F.C0],F.options,F);
                catch
                    if dispon, disp(sprintf('...unexpected error ''%s'' for F =',lasterr)), disp(F), end
                    if dispon, disp(sprintf('==>   STOP at %s', datestr(now))), end
                    tC		= F.t;
                    Code	= ones(size(tC));
                end
                timeout_ode = timeout_odeold;
            else
                if dispon, disp('... decrease maxstep'), end
                F.options.MaxStep = F.options.MaxStep/100;
                tic
                [tC,Code]	= ode15s(@odediff,F.t,[F.C0],F.options,F);
            end
        end
    else
        tC = F.t;
        Code = F.C0 * ones(size(F.t));
    end
    if dispon, disp(['SENSode: end in ' num2str(toc) ' s']), end
    Codei	= interp1(tC,Code,ti,method);
    %fodei	= - ndf(ti,Codei);
    timeout_odeold = timeout_ode; timeout_ode = 10; tic
    fodei	= - odediff(Codei,F);
    timeout_ode = timeout_odeold;
end

% Outputs
switch lower(mode), case {'both' 'ode'}
	if nargout>0
		resode.t = ti;
		resode.C = Codei;
		resode.f = fodei;
		resode.fD = fodei/F.Bi;
	end
end
switch lower(mode), case {'both' 'pde'}
	if nargout>0 & strcmp(lower(mode),'pde'), resode	= []; end
    if nargout>1
        if ~errpde
            respde.x = F.x;
            respde.Cx = C;
            respde.tC = F.t;
            respde.t = ti;
            respde.C = Cmi;
            respde.fc = fci;
            respde.f = fi;
            respde.fD = fi/F.Bi;
        else
            switch lower(mode)
            case 'pde', respde = sens('ode',F,ploton);
            case 'both', respde = resode;
            end
        end
	end
end
	
if ploton
	switch lower(mode), case {'both' 'pde'}	
		figure(1)
		subplot(231), hold on, plot(F.x,C), xlabel('x'), ylabel('C')
		subplot(232), hold on, plot(ti,Cmi,'b-'), xlabel('Fo')
		subplot(233), hold on, plot(sqrt(ti),Cmi,'b-'), xlabel('Fo^{1/2}')
		subplot(234), hold on, plot(ti,fi,'b-'), ylabel('j')
		subplot(236), hold on, plot(Cmi,fi,'b-',Codei,fodei,'m-'), xlabel('C'), ylabel('j')
		figure(2)
		subplot(121), hold on, plot(ti/F.Bi,Cmi,'k:'), xlabel('Fo/Bi'), ylabel('C')
		subplot(122), hold on, plot(Cmi,1/F.Bi*fi,'k:'), xlabel('C'), ylabel('j/Bi')
	end
	switch lower(mode), case {'both' 'ode'}
		figure(1)
		subplot(232), hold on, plot(ti,Codei,'r-')
		subplot(233), hold on, plot(sqrt(ti),Codei,'r-')
		subplot(234), hold on, plot(ti,fodei,'r-')
		subplot(236), hold on, plot(Codei,fodei,'r-')
		figure(2)
		subplot(121), hold on, plot(ti/F.Bi,Codei,'k-')
		subplot(122), hold on, plot(Codei,1/F.Bi*fodei,'k-')
	end
	switch lower(mode), case 'both'
		figure(1)
		subplot(235), hold on, plot(Cmi,Codei,'ro',[min(Cmi) 1],[min(Cmi) 1],'k-')
		xlabel('C_{pde}'), ylabel('C_{pde}')
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numerical differentiation
% ordre1 [h6y(8)/140]
% ordre2 [h6y(8)/560]
function dydt = ndf(t,y,ordre,dydt0)
if nargin<3, ordre = 1; end
if nargin<4, dydt0 = []; end
dt	= t(2)-t(1);
switch ordre
case 1
	D	= [-1 +9 -45 0 +45 -9 +1]'/(60*dt);
	D0	= [-137 300 -300 200 -75 12]'/(60*dt);
	D1	= -flipud(D0);
case 2
	D	= [2 -27 +270 -490 270 -27 +2]'/(180*dt.^2);
	D0	= [45 -154 +214 -156 61 -10]'/(12*dt.^2);
	D1 = flipud(D0);
end
% yfull	= [repmat(y(1),3,1) ; y ; repmat(y(end),3,1) ];
% dydt	= [yfull(1:end-6) yfull(2:end-5) yfull(3:end-4) yfull(4:end-3) yfull(5:end-2) yfull(6:end-1)  yfull(7:end)] * D;
dydt = 	[
		[y(1:3) y(2:4) y(3:5) y(4:6) y(5:7) y(6:8) ] * D0;
		[y(1:end-6) y(2:end-5) y(3:end-4) y(4:end-3) y(5:end-2) y(6:end-1)  y(7:end)] * D
		[y(end-7:end-5) y(end-6:end-4) y(end-5:end-3) y(end-4:end-2) y(end-3:end-1)  y(end-2:end)] * D1
	];
if any(dydt0), dydt(1) = dydt0; end
dy = diff(dydt);
i = intersect(find(dy*sign(mean(sign(dy)))<0),2:10);
if any(i), ni = setdiff(1:length(t),i); dydt(i) = interp1(t(ni),dydt(ni),t(i),'cubic','extrap'); end
%  dydt 	= interp1(t(6:end-5),dydt(6:end-5),t,'cubic','extrap');
% if mean(diff(dydt(end-3:end)))*mean(diff(dydt(5:end-4)))<0
% 	disp('**')
% end

% ODE
function dCdt = odediff(t,C,S)
global timeout_ode
isvectorized = 0;
if toc>timeout_ode, eval('erron'), end
if nargin<3
	isvectorized = 1;
	S = C;
	C = t;
end
if ~isvectorized
	j1		= S.Bi * ( (S.K+S.L)*C(1) - S.L*S.C0 ) / ( 1 + 2/6 * S.Bi * S.K);
	b		= S.Bi * S.K * sqrt( 3/2 * (S.C0 - C(1) ) );
	c		= S.Bi * ( (S.K-S.L) * S.C0 + S.L * C(1) );
	d		= b^2 + 4 * c;
    j2		= (( -b + sqrt(d) ) / 2)^2;
    if j2<=0
        dCdt = -j1;
    else
        if S.C0-sqrt(3/2*(S.C0-C)/j2)>S.xcrit,		dCdt = -j2;
        else,										dCdt = -j1; end
    end
else
	j1		= S.Bi * ( (S.K+S.L)*C - S.L*S.C0 ) / ( 1 + 2/6 * S.Bi * S.K);
	b		= S.Bi * S.K * sqrt( 3/2 * (S.C0 - C ) );
	c		= S.Bi * ( (S.K-S.L) * S.C0 + S.L * C );
	d		= b.^2 + 4 * c;
	j2		= (( -b + sqrt(d) ) / 2).^2;
    j2(j2<=0) = NaN;
	i2		= S.C0-sqrt(3/2*(S.C0-C)./j2)>S.xcrit;
	i1		= S.C0-sqrt(3/2*(S.C0-C)./j2)<=S.xcrit;
	dCdt	= zeros(size(C));
	dCdt(i2) = -j2(i2);
	dCdt(i1) = -j1(i1);
end


% PDE problem :  [C,F,S] = PDEFUN(X,T,U,DUDX,P1,P2...)
function [c,f,s] = pdediff_PDE(x,t,U,dUdx,S)
global timeout_pde
if toc>timeout_pde, eval('erron'), end
c = 1;
f = dUdx;
s = 0;

% IC : U = ICFUN(X,P1,P2...)
function U0 = pdediff_IC(x,S)
U0 = ones(size(x))*S.C0;

% BC : [PL,QL,PR,QR] = BCFUN(XL,UL,XR,UR,T,P1,P2...)
function [pa,qa,pb,qb,purge] = pdediff_BC(xa,Ua,xb,Ub,t,Fc,S)
global Hsmooth
% pervious contact on the left
% pb = 0;
% qb = 1;
% pa = S.Bi*( S.K*(Ua-0) - S.L * -Fc(1) );
% qa = -1;
% pervious contact on the right
pa = 0;
qa = 1;
pb = S.Bi*( S.K*(Ub-0) - S.L * -Fc(2) );   
qb = 1;

% General 1D PDE formulation (without advection)
%
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
