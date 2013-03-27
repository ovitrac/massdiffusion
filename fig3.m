% FIGURE 3: compares Energy, Erfc and Sagiv's solutions
% INRA\Daniel Goujot, Olivier Vitrac (c)2013 - Last revision 27/03/2013

% This script illustrates how to use sensc for calculating solutions

%% Calculations
F=sensc('init');
F.Bi=100; K0=1; F.L=1; filename='fig3';
%F.Bi=1000; K=10; F.L=.1;% sensca ne marche pas dans ce contexte :-( TODO bug a corriger.
F.x=1;
F.nmax=10;
nmaxerfc=10;
F.iso.S.Cs=[.0 5];% 5 pour etre sur qu'il n'y a pas de coin.
F.iso.S.K=[1/K0 1/K0]; % sensc se trompe dans le sens de K, pas sensca_shorttimes.
F.t=10.^[[-9:.1:-3.1] [-3:.01:2]]';
[resode,respde,cosinus]=sensc('anAlytic',F,false,false);% le A sert à ce que long.Cx contienne aussi les temps relatifs à f.t ne soit pas vidé.
[a,b,saguiv]=sensc('anALyTic',F,false,false);
Fref=F;Fref.nmax=5000;[a,b,saguivref]=sensc('anALyTic',Fref,false,false);
erfcs=sensca_shorttime(F.x,0,cosinus.tC',F.Bi,K0,cosinus.Cx(end)*K0,0,cosinus.Cx(end),1,F.L,nmaxerfc);
erfcs(~isreal(erfcs))=NaN;
if F.x==1 && nmaxerfc>=33
    varC=(1-cosinus.Cx(end));% amplitude de variation de C_s^t, cas p=0.
    for terf=[cosinus.tC';erfcs'];
        controle_erfcs=(sensca_shorttime(cosinus.x(end)',0,terf(1),F.Bi,K0,cosinus.Cx(end)/K0/varC,0,cosinus.Cx(end)/varC)-cosinus.Cx(end)/varC)*varC+cosinus.Cx(end);
        if abs(controle_erfcs-terf(2))>eps^.2, error('no convergence, ask Daniel or Olivier'), end
    end
end

%% Do plot
close all
% common definitions
formatfig(gcf,'figname','fig3','paperposition',[ 2.5225    7.4226   15.9391   14.8323])
leg = {{'Energy solution' 'N=10'},{'Erfc solution' 'with N=10'},{'Sagiv solution' 'with N=10'},{'Sagiv solution' 'with N=5000' '(reference)'}};
prop = plotpub(     'color', rgb({'Red','Blue','Green','Purple'}),...
                'linewidth', [2 1 3 4],...
                'linestyle', {'--','-','-.',':'},...
                   'xscale', 'log',...
                     'xlim', F.t([1 end]),...
                   'nmarker', 10 );
xplot = {cosinus.tC,cosinus.tC,saguiv.tC,saguivref.tC};
yplot = {cosinus.Cx(:,end),erfcs,saguiv.Cx(:,end),saguivref.Cx(:,end)};
% main plot
main = gca;
hp = plotpub(xplot,yplot,prop);
ylabel(sprintf('c(%0.3g,t)',F.x),'fontsize',18)
xlabel('t-t_0','fontsize',18)
formatax(gca,'fontsize',12)
hl = legendpub(hp,leg,axes('position',[0.15    0.33    0.2000    0.35]),'','fontsize',12);
% inset
inset = axes('position',[0.4911    0.5810    0.3893    0.3095]);
plotpub(xplot,yplot,prop,'nmarker',50,'xlim',[.001 .1],'ylim',[.1 .3]);
formatax(inset,'fontsize',8)
% titles
titles([main inset],[],'suffix',')','x',[0.95 .9],'y',[-0.02 .1],'fontsize',12)

%% print
pdfoutput = [get(gcf,'filename') '.pdf'];
print_pdf(600,pdfoutput,[],'nocheck')
if isunix, system(['evince ' pdfoutput ' &']); else winopen(pdfoutput), end