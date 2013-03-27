function res = senspatankarT(F,ploton,dispon)
%SENSPATANKART simulates non-isotherm transfer through n layers using a modified Patankar Method (see p 45)
%   Syntax: res = senspatankarT(F [,ploton,dispon]) 
%
% The dimensionless formulation is similar to SENSPATANKAR except:
%   k0_T0, k_T0 and D_T0 replace k0, k and D respectively
%   Eak and EaD are the activation energies of k (including k0) and D respectively
%   the normalization is based on on values at T0
%   T0: is the starting temperature (reference temperature for all subsequent normalization)
%   NOTE: the true equilibrium value cannot conventionaly be guessed a priori as T(t->+oo) remains unknown
%
%   To assign a temperature profile (temperature variation around T0) use either:
%       F.diffT = user function in t (assumed to be enough smooth in t)
%      or
%       F.diffT = @temp_profile; %which uses an internal function constructed with F.temp_profile
%       example: F.temprofile = struct('start_temp',{0 110 110},'final_temp',{110 110 0},'duration',{0.3 1 0.3})
%       NB: all temperatures are relative t0 and durations assumed a dimensionless formulation
%
%   To restart a simulation from a previously calculated solution, use: 
%       F.restart = struct('x',some values,'C',some values);
%       F.restart.method defines the method to interpolate from the previous solution (default = 'linear')
%
%   Transtion temperatures can be assigned by describing explicitely all transition temperatures (denoted T1)
%       F.T1 = mxn array, where m = number of transition temperatures and n is the number of layers 
%              T1 values must be increasing and T1 bounds must include T0
%       F.Eak and F.EaD are defined accordingly as mxn arrays
%           where Eak(i) and EaD(i) stand for the activation energies between T1(i) and T1(i+1)
%
% Example: using transition temperatures
%           S = senspatankarT;
%           S.Bi    = 0;
%           S.k_T0  = [1 1 1];
%           S.Eak   = [30 30 30 ; 30 90 30]*1e3;
%           S.D_T0  = [2e-15 1e-16 2e-15];
%           S.EaD   = [85 85 85 ; 85 200 85]*1e3;
%           S.k0_T0 = 1;
%           S.Eak0  = 30e3;
%           S.l     = [50 20 120]*1e-6;
%           S.C0    = [1000 0 0];
%           S.T0    = 273.15+20;
%           S.T1    = [293 293 293;383 383 383];
%           S.temp_profile = struct('start_temp',{0 110 110},'final_temp',{110 110 0},'duration',{0.3 1 0.3});
%           S.diffT = '@temp_profile';
%           res = senspatankarT(S);
%   
%   


% AMCOR 1.0 - 04/03/09 - INRA\Olivier Vitrac - rev. 31/03/09

% Revision history
% 16/03/09 add temp_profile, restart
% 30/03/09 add transition
% 31/03/09 add several transitions
% 21/04/09 add F.Bi_T0 and F.EaBi

% definitions
global timeout
ploton_default = false;
dispon_default = false;
timeout		= 800; % s
% options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Stats','yes','Initialstep',1e-5,'Maxstep',.05,'Maxorder',5);
options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.01,'Maxorder',2);
Fdefault	= 	struct(...
				'Bi_T0'		, 	1e3,...	Biot [hm.L1/D]
                'EaBi'      ,   60e3,...
				'k_T0'      ,	[1 1 1],...[0.5 3 2],...	ki, i=1 (layer in contact with the liquid)
                'Eak'       ,   [30 30 30]*1e3,...
                'D_T0'      ,   [1e-16 1e-16 1e-14],... diffusion coefficient
                'EaD'       ,   [85 85 85]*1e3,... Activation energy
                'k0_T0'     ,   1,... 0 = liquid
                'Eak0'      ,   30e3, ...
                'l'         ,   [50 20 120]*1e-6,...[50 20 10 120]*1e-6,... m
				'L'			,	200/1800,...	dilution factor (respectively to iref)
				'C0'		,	[0 500 0],...	initial concentration in each layer
                'T0'        ,   273.15+20,...
                'temp_profile', struct('start_temp',{0 110 110},'final_temp',{110 110 0},'duration',{0.3 1 0.3}), ... 
                'diffT'     ,   '@temp_profile',... @(t) 0*t
                'restart'   ,   [], ... restart solution (example: struct('x',some values,'C',some values))
				'options'	,	options...
					); % if iref is missing, it is indentified
Fdefault.t	= [0:.00001:.005 .01:.01:.1 .11:.1:5]; %0:.00001:.005; %[0:.00001:.005 .01:.01:.1 .11:.1:5]';
if Fdefault.Bi_T0<10, Fdefault.t = Fdefault.t/(Fdefault.Bi_T0/10); end
method		= 'cubic'; %'cubic';
nmesh_default  = 250; % number of nodes for a layer of normalized thickness 1
nmeshmin_default    = 10;
n_default	= 1e4;
zero = 1000*eps;
R = 8.31447215;

% Declaration for nested functions (to be sure that they are shared by all nested functions)
[D,k,k0,T,D_T0,k_T0,EaD,Eak] = deal([]);

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
if ~isfield(F,'nmeshmin'), F.nmeshmin = nmeshmin_default; end
if ~isfield(F,'options'), F.options = options; end
if isempty(ploton), ploton = ploton_default; end
if isempty(dispon), dispon = dispon_default; end

% temp_profile check
if ischar(F.diffT), F.diffT = eval(F.diffT); end
if strcmpi(func2str(F.diffT),[mfilename '/temp_profile'])
    temp_profile_on = true;
    F.tevents = cumsum([0 [F.temp_profile.duration]]);
else
    temp_profile_on = false;
end

% physical check
m = Inf;
for prop = {'D_T0' 'EaD' 'k_T0' 'Eak' 'C0'};
    if strcmp(prop{1},'C0')
        F.(prop{1}) = F.(prop{1})(F.(prop{1})>=0);
    elseif (strcmp(prop{1},'EaD') && size(F.(prop{1}),1)>1) || (strcmp(prop{1},'Eak') && size(F.(prop{1}),1)>1)
        if any(F.(prop{1})<0), error('inconsistent ''%s'' values',prop{1}), end
    else
        F.(prop{1}) = F.(prop{1})(F.(prop{1})>0);
    end
    m = min(m,length(F.(prop{1})));
end
F.m = m;
% check whether additional restrictions apply            %added 24/01/07
if isfield(F,'ivalid')
    F.ivalid = F.ivalid(F.ivalid<m);
    F.a = F.a(F.ivalid);
    F.lref = F.lref(F.ivalid);
    F.lrefc = F.lrefc(F.ivalid);
    F.C0 = F.C0(F.ivalid);
    F.l = F.l(F.ivalid);
    F.D = F.D_T0(F.ivalid);
    F.k = F.k(F.ivalid);
    m = length(F.ivalid);
    F.m = m;
end

% transition temperature T1
nlayer = size(F.D_T0,2);
F.nexttransitionatT0 = ones(1,nlayer);
if isfield(F,'T1')
    nT1 = size(F.T1,1);
    if size(F.EaD,1)~=nT1, error('the number of rows in EaD, %d, must be equal the number of transitions %d',size(F.EaD,1),nT1); end
    if size(F.Eak,1)~=nT1, error('the number of rows in Eak, %d, must be equal the number of transitions %d',size(F.Eak,1),nT1); end
    if (size(F.T1,2)~=nlayer), error('the number of transitions does not match the number of layers %d',nlayer), end
    if (size(F.EaD,2)~=nlayer) || (size(F.Eak,2)~=nlayer), error('the number of activation energies Eak or EaD does not match the number of layers %d',nlayer), end
    F.D_T1 = repmat(F.D_T0,nT1,1);
    F.k_T1 = repmat(F.k_T0,nT1,1);   
    for ilayer = 1:nlayer
        if any(diff(F.T1(:,ilayer))<=0), error('layer %d: all transition temperatures T1 must be increasing',ilayer), end
        istart = find(F.T0>=F.T1(:,ilayer),1,'first');
        if isempty(istart), error('layer %d: the diagram is not defined for temperatures lower than T0=%0.4g',ilayer,F.T0); end
        F.nexttransitionatT0(ilayer) = istart;
        % search transitions >= T0
        ind = find(F.T1(:,ilayer)>=F.T0);
        if isempty(ind), error('layer %d: the diagram is not defined for temperatures greater than T0=%0.4g',ilayer,F.T0); end
        Dtemp = F.D_T0(ilayer) * exp( -F.EaD(istart,ilayer)/R .* (1/F.T1(istart,ilayer) - 1/F.T0 ) );
        ktemp = F.k_T0(ilayer) * exp( -F.Eak(istart,ilayer)/R .* (1/F.T1(istart,ilayer) - 1/F.T0 ) );
        for iT1=1:length(ind)
            F.D_T1(ind(iT1),ilayer) = Dtemp * exp( -F.EaD(ind(iT1)-1,ilayer)/R .* (1/F.T1(ind(iT1),ilayer) - 1/F.T1(ind(iT1)-1,ilayer)) );
            F.k_T1(ind(iT1),ilayer) = ktemp * exp( -F.Eak(ind(iT1)-1,ilayer)/R .* (1/F.T1(ind(iT1),ilayer) - 1/F.T1(ind(iT1)-1,ilayer)) );
            Dtemp = F.D_T1(ind(iT1),ilayer);
            ktemp = F.k_T1(ind(iT1),ilayer);
        end
        % search transitions < T0
        ind = find(F.T1(:,ilayer)<F.T0);
        for iT1=length(ind):-1:1
            F.D_T1(ind(iT1),ilayer) = F.D_T1(ind(iT1)+1,ilayer) * exp( -F.EaD(ind(iT1),ilayer)/R .* (1/F.T1(ind(iT1),ilayer) - 1/F.T1(ind(iT1)+1,ilayer)) );
            F.k_T1(ind(iT1),ilayer) = F.k_T1(ind(iT1)+1,ilayer) * exp( -F.Eak(ind(iT1),ilayer)/R .* (1/F.T1(ind(iT1),ilayer) - 1/F.T1(ind(iT1)+1,ilayer)) );
        end
    end
    istransitiondefined = true;
else
    istransitiondefined = false;
end

% constructor
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

% Initial parameters
[k0,k,D,Bi,a,iref,lref,lrefc] = init_activation(0,F);

% false equilibrium value
C0eq = sum(lref*F.L.*F.C0)/(1+sum((k0./k .* lref*F.L)));
peq = k0 * C0eq;

% mesh generation, populate normalized transport properties at each node
X = ones(F.m,1);
for i=2:m
    X(i) = X(i-1)*(a(i-1)*lref(i))/(a(i)*lref(i-1));
end
X = max(F.nmeshmin,ceil(F.nmesh*X/sum(X)));
X = round((X/sum(X))*F.nmesh); % update
F.nmesh = sum(X); % roundoff
[xmesh,D_T0,k_T0,EaD,Eak,de,dw,he,hw,C0] = deal(zeros(F.nmesh,1));
allnodes = (1:F.nmesh)' ;
if istransitiondefined
    nT1 = size(F.T1,1);
    [Eak,EaD,T1,D_T1,k_T1] = deal(repmat(Eak,1,nT1));
    isnextT1 = false(F.nmesh,1);
    inextT1 = xmesh;
end
j = 1; x0 = 0;
for i=1:F.m
    ind = j+(0:X(i)-1);
    de(ind) = lref(i)/(2*X(i));
    dw(ind) = lref(i)/(2*X(i));
    D_T0(ind)= D(i);
    k_T0(ind)= k(i);
    if istransitiondefined
        for iT1 = 1:nT1
            EaD(ind,iT1) = F.EaD(iT1,i);
            Eak(ind,iT1) = F.Eak(iT1,i);
            D_T1(ind,iT1) = F.D_T1(iT1,i)/F.D_T0(iref);
            k_T1(ind,iT1) = F.k_T1(iT1,i)/F.k0_T0;
            T1(ind,iT1)   = F.T1(iT1,i);
        end
    else
        EaD(ind) = F.EaD(i);
        Eak(ind) = F.Eak(i);
    end
    C0(ind) = F.C0(i)/C0eq;
    xmesh(ind) = linspace(x0+dw(ind(1)),lrefc(i)-de(ind(end)),X(i));
    x0 = lrefc(i);
    j = ind(end)+1;
end

% use a previous solution (if any)
if isfield(F,'restart') && ~isempty(F.restart) && isstruct(F.restart) && isfield(F.restart,'x') && isfield(F.restart,'C')
    if ~isfield(F.restart,'method'), F.restart.method = 'linear'; end
    C0 = interp1(F.restart.x/F.restart.x(end),F.restart.C,xmesh/xmesh(end),F.restart.method)/C0eq; 
end

% Integration
dCdt = @(t,C)(rigidity(t)*C); % function to integrate
F.options.Vectorized = 'on';
if temp_profile_on, F.options.Events = @events; end
[t,C] = ode15s(dCdt,F.t,[0;C0],F.options); % integration
CF = C(:,1);
C  = C(:,2:end);

% interpolation at each interface
Ce = zeros(length(t),F.nmesh);
Cw = zeros(length(t),F.nmesh);
for i=1:length(t)
    rigidity(t(i)); %[k0,k,D]=update_transport(t(i));
    Ce(i,1:end-1) = C(i,1:end-1) - (de(1:end-1).*he(1:end-1).*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(1:end-1) )';
    Ce(i,end)     = C(i,end-1);
    Cw(i,2:end)   = C(i,2:end)   + (dw(2:end)  .*hw(2:end)  .*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(2:end)   )';
    Cw(i,1)       = C(i,1)       +  dw(1)       *hw(1)       *( k0       /k(1)      * CF(i)         - C(i,1)      )  / D(1);
end
Cfull = reshape([Cw;C;Ce],length(t),3*F.nmesh);
xw = xmesh-dw+zero;
xe = xmesh+de-zero;
xfull = [xw';xmesh';xe']; xfull = xfull(:);

% outputs
F.peq = peq;
F.lref = lref;
F.lrefc = lrefc;

res.C = interp1(t,trapz(xfull,Cfull,2)*C0eq,ti,method)/xfull(end); % av conc
res.t = ti; % time base
res.p = NaN; %interp1(t,repmat(kfull',length(t),1).*Cfull*C0eq,ti,method); % to fast calculations
res.peq = peq;
res.V  = F.lrefc(end); %*F.l(F.iref);
res.C0 = C0*C0eq;
res.F  = F;
res.x  = xfull; %*F.l(F.iref);
res.Cx = Cfull*C0eq; %interp1(t,Cfull*C0eq,ti,method);
res.tC = t;
res.CF = interp1(t,CF*C0eq,ti,method);
res.fc = res.CF/F.L;
res.f = zeros(size(ti));
for i=1:length(ti)
    [k0,k,D,Bi]=update_transport(ti(i));
    hw(1) = 1/( (k0/k(1))/(Bi) + dw(1)/(D(1))   );
    res.f(i)  = interp1q(t,hw(1) * ( k0/k(1) * CF - C(:,1) ) * C0eq,ti(i)); %interp1(t,hw(1) * ( k0/k(1) * CF - C(:,1) ) * C0eq,ti(i),method);
end
res.timebase = res.F.l(iref)^2/(res.F.D_T0(iref));

% ------------------------------------------------------
% Nested functions
% (see documentation how to use them, works only with Matlab versions later than 7.5 )
%------------------------------------------------------

function [value,isterminal,direction] = events(t,y)
    value = min(abs(F.tevents-t));
    isterminal = 0; % the integration is not stopped
    direction = 0;  % all zeros have to be located
end
    

function T = temp_profile(t)
    iseq = find(t<=F.tevents,1,'first')-1;
    if any(iseq) && iseq>0
        T = F.temp_profile(iseq).start_temp + (t-F.tevents(iseq)) * ...
            (F.temp_profile(iseq).final_temp - F.temp_profile(iseq).start_temp) / F.temp_profile(iseq).duration;
    else
        T = 0;
    end
    %dispf('t=%0.4g, T=%0.4g',t,T); pause %for debug
end


function [k0,k,D,Bi]=update_transport(t)
    % assume that all data were normalized
    T = F.T0 + F.diffT(t);
	k0 = exp(-F.Eak0/R*(1/T-1/F.T0)); % 1 at T=T0
    Bi = F.Bi_T0*exp(-F.EaBi/R*(1/T-1/F.T0)); % Bi at T=T0
    if istransitiondefined
          isnextT1(:) = false;
          for eachtransition=nT1:-1:1
              newnextT1 = (T>=T1(:,eachtransition)) & ~isnextT1;
              isnextT1(newnextT1) = true;
              inextT1(newnextT1) = eachtransition;
          end
          indT1ref = sub2ind([F.nmesh,nT1],allnodes,inextT1);
          k =  k_T1(indT1ref)  .* exp(-Eak(indT1ref)/R.*(1/T-1./T1(indT1ref)));
          D =  D_T1(indT1ref)  .* exp(-EaD(indT1ref)/R.*(1/T-1./T1(indT1ref)));
%         state = (sign(T1-F.T0)*T)>(sign(T1-F.T0)*T1);
%         k = state .* (k_T0  .* exp(-Eak(:,1)/R*(1./T1-1/F.T0)) .* exp(Eak(:,2)/R*(1./T-1/F.T1))) + ...
%         (1-state) .* (k_T0  * exp(-Eak(:,1)/R*(1./T-1/F.T0))) ;
%         D = state .* (D_T0  .* exp(-EaD(:,1)/R*(1./T1-1/F.T0)) .* exp(EaD(:,2)/R*(1./T-1/F.T1))) + ...
%         (1-state) .* (D_T0  * exp(-EaD(:,1)/R*(1./T-1/F.T0))) ;
    else
        k =  k_T0  .* exp(-Eak/R*(1/T-1/F.T0));
        D =  D_T0  .* exp(-EaD/R*(1/T-1/F.T0));
    end
end


function A=rigidity(t)
    [k0,k,D,Bi]=update_transport(t);
    % equivalent conductances
    he = zeros(F.nmesh,1); hw = he;
    hw(1) = 1/( (k0/k(1))/(Bi) + dw(1)/(D(1))   ); % pervious contact
    for inode=2:F.nmesh
        hw(inode) = 1/( (de(inode-1)/D(inode-1))*(k(inode-1)/k(inode)) + dw(inode)/D(inode)  );
    end
    he(1:F.nmesh-1) = hw(2:F.nmesh); %flux continuity
    he(end) = 0; % impervious BC
    if ~nargout, return, end
    % Assembling
    A = zeros(F.nmesh+1,F.nmesh+1);
    % fluid phase = node 0 (position 1)
    A(1,1:2) = F.L*hw(1) * [-k0/k(1)
                         1           ]';
    % node 1 (position 2)
    A(2,1:3) = 1/(dw(1)+de(1)) * [ hw(1)*k0/k(1)
                                  -hw(1)-he(1)*k(1)/k(2)
                                   he(1) ]';                     
    % node inode<n (position inode+1)
    for inode=2:F.nmesh-1
        A(inode+1,inode:inode+2) = 1/(dw(inode)+de(inode)) * [ hw(inode)*k(inode-1)/k(inode)
                                         -hw(inode)-he(inode)*k(inode)/k(inode+1)
                                          he(inode) ]';
    end
    % node n (position n+1)
    inode = F.nmesh;
    A(end,end-1:end) = 1/(dw(inode)+de(inode)) * [ hw(inode)*k(inode-1)/k(inode)
                                          -hw(inode)]';
end


function [k0,k,D,Bi,a,iref,lref,lrefc] = init_activation(t,F)
% initialize the normalization of data
    T = F.T0 + F.diffT(t);
    % renormalization (fields X1->Xn)
	k0 = F.k0_T0 .* exp(-F.Eak0/R*(1/T-1/F.T0))  /F.k0_T0; % exactly 1 when T=T0
    Bi = F.Bi_T0 .* exp(-F.EaBi/R*(1/T-1/F.T0))  /F.k0_T0; % exactly Bi when T=T0
    k =  F.k_T0  .* exp(-F.Eak(sub2ind(size(F.Eak),F.nexttransitionatT0,1:nlayer))/R*(1/T-1/F.T0))   /F.k0_T0;
    D =  F.D_T0  .* exp(-F.EaD(sub2ind(size(F.EaD),F.nexttransitionatT0,1:nlayer))/R*(1/T-1/F.T0));
    a = D./k;
    % reference layer
    if isfield(F,'iref'), iref = F.iref; else [crit,iref] = min( a./F.l(1:F.m) ); end
    a = a./a(iref);
    lref  = F.l./F.l(iref);
    lrefc = cumsum(lref);
    D = D/D(iref);
end

end