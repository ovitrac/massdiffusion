function [c_x_t_EQ46,errEQ46,c_x_t_EQ46_x,c_x_t_EQ46_mat,c_x_t_EQ46_x_mat,c_x_t_EQ46_t,c_x_t_EQ46_x_t,c_x_t_EQ46_x_t_t,c_x_t_EQ46_x_x,c_t_EQ46,C_f_t]=sensca_shorttime(x,t0,t,B,K0,Cftinfty,p,ctinfty,varargin)% varargin: ct0, L, nmax.
% Compute the concentration at point x at time t in condensed.
% Based on formulas (60), (61) and (62) of senscnonlin_v77.docx
% The condensed region is between 0<x<1.
% The fluid region is between 1<x<2. if x<50, then 100 is added to x and formula (27) of senscnonlin_v87.pdf is used instead.
% ctinfty and Cftinfty are the final concentration in condensed and fluid parts.
% At time t0, the initial concentration profile is ctinfty+exp(-p*x).
% t should be close to initial time t0, the error is in exp(-(1./(t-t0))).
% B is the Biot number, multiplied by L.
% K0 is the inverse of (the coefficient de partage, divided by L.)
% p=0 is a constant initial profile. p may be complex (useful for cosine, ...)
% Dilution factor is 1. See senscnonlin_v68.docx for details to use this for any other dilution factor.
% if x==1, then bulk concentrtions C_f_t/L and C_s_t are also computed.
% use matlab2007b to make analytical tests.
% c_x_t_EQ46=sum(c_x_t_EQ46_mat(:,1)) to make partial sums, terms contain exp(c_x_t_EQ46_mat(:,2))*erfc(c_x_t_EQ46_mat(:,3)) which may be used to classify them.
% if varargin is present, it is interpreted as ct0, L, nmax; the initial concentration profile is ctinfty+A*exp(-p*x) with A such that the mean initial concentration is ct0, and x and t may be vectors. L, if given, is automatically taken in account, and B, K0 and C_f_t are no more dependent of L. nmax is, if positive, the number of erfc computed, and if negative, the number of terms of EQ46 that are computed.

if ~isempty(varargin)
  ct0=varargin{1};
  L=1; if length(varargin)>1; L=varargin{2}; if L~=1; warning('marche pas encore pour L different de 1, il y a un bug'); keyboard; end; end;
  nmax=1e9; if length(varargin)>2; nmax=varargin{3}; end;% if >0, nb of erfc, if <0, nb of monomials in EQ46.
  C_f_t=NaN(length(t),1);
  c_t_EQ46=NaN(length(t),1);
  if isempty(t) || isempty(x); return; end;
  for not=1:length(t);
    for nox=1:length(x);
      varCx0=ct0-ctinfty;
      if abs(p)>=eps^.8
        varCx0=varCx0/((1-exp(-p))/p);
      end
      clear a;% gestion d'un nb variable d'arguments
      if x(nox)<1;
        a{min(8,nargout)}=1;
      else
        a{nargout}=1;
      end
      if length(nmax)==1 && nmax==1e9 ; else a{5}=1; end;
      [a{:}]=sensca_shorttime(x(nox),t0,t(not),B*L,K0*L,Cftinfty/L/varCx0,p,ctinfty/varCx0);
      if length(nmax)==1 && nmax==1e9 ; a{5}=1; end;
      [c_x_t,err,c_x_t_x,c_x_t_EQ46_mat,c_x_t_EQ46_x_mat]=deal(a{1:5});
      if nargout>5 && length(nmax)==1 && nmax==1e9; c_x_t_EQ46_t(not,nox)=a{6}*varCx0; else c_x_t_EQ46_t=[]; end;
      if nargout>6 && length(nmax)==1 && nmax==1e9; c_x_t_EQ46_x_t(not,nox)=a{7}*varCx0; else c_x_t_EQ46_x_t=[]; end;
      if nargout>7 && length(nmax)==1 && nmax==1e9; c_x_t_EQ46_x_t_t(not,nox)=a{7}*varCx0; else c_x_t_EQ46_x_t_t=[]; end;
      if nargout>8 && length(nmax)==1 && nmax==1e9; c_x_t_EQ46_x_x(not,nox)=a{8}*varCx0; else c_x_t_EQ46_x_x=[]; end;
      c_x_t_EQ46_x(not,nox,1)=c_x_t_x*varCx0;
      % [c_x_t_EQ46,errEQ46,c_x_t_EQ46_x,c_x_t_EQ46_mat,c_x_t_EQ46_x_mat,c_x_t_EQ46_t,c_x_t_EQ46_x_t,c_x_t_EQ46_x_x,c_t_EQ46,C_f_t]
      errEQ46(not,nox)=err*varCx0;
      c_x_t_EQ46(not,nox,1)=c_x_t*varCx0;
      c_x_t_EQ46_x(not,nox,1)=c_x_t_x*varCx0;
      if length(nmax)>1 || nmax<1e9
	[sortedlines,unu_sed,nomat2noargerfc]=unique([-majorant_logerfc(c_x_t_EQ46_mat(:,3)),c_x_t_EQ46_mat(:,3)],'rows');
	[unused,noargerfc2nomat]=sortrows([-majorant_logerfc(c_x_t_EQ46_mat(:,3)),c_x_t_EQ46_mat(:,3),-(c_x_t_EQ46_mat(:,2))]);
	% x=-20:20;y=-20:20;[X,Y]=meshgrid(x,y);Z=X;Z(:)=erfc_complex(X(:)+1i*Y(:));[h,c]=contour(X,Y,log10(eps+abs(Z)));

        for nbnmax=1:length(nmax)
          if nmax(nbnmax)>0;
            % nb of erfc, sorted by abs(their argument)
            c_x_t_EQ46(not,nox,nbnmax)=sum(c_x_t_EQ46_mat(nomat2noargerfc<=nmax(nbnmax),1))*varCx0;
          else
            % nb of monomials. sorted by abs(argument of erfc), with minor sort on (argument of exp)
            c_x_t_EQ46(not,nox,nbnmax)=sum(c_x_t_EQ46_mat(noargerfc2nomat<=-nmax(nbnmax),1))*varCx0;
          end
        end
	[sortedlines,unu_sed,nomat2noargerfc]=unique([-majorant_logerfc(c_x_t_EQ46_x_mat(:,3)),c_x_t_EQ46_x_mat(:,3)],'rows');
	[unused,noargerfc2nomat]=sortrows([-majorant_logerfc(c_x_t_EQ46_x_mat(:,3)),c_x_t_EQ46_x_mat(:,3),-(c_x_t_EQ46_x_mat(:,2))]);
        for nbnmax=1:length(nmax)
          if nmax(nbnmax)>0;
            % nb of erfc, sorted by abs(their argument)
            c_x_t_EQ46_x(not,nox,nbnmax)=sum(c_x_t_EQ46_x_mat(nomat2noargerfc<=nmax(nbnmax),1))*varCx0;
          else
            % nb of monomials. sorted by abs(argument of erfc), with minor sort on (argument of exp)
            c_x_t_EQ46_x(not,nox,nbnmax)=sum(c_x_t_EQ46_x_mat(noargerfc2nomat<=-nmax(nbnmax),1))*varCx0;
          end
        end
      end
      if x(nox)==1;
        C_f_t(not,1,1:length(nmax))=(K0*c_x_t_EQ46(not,nox,:)+ c_x_t_EQ46_x(not,nox,:)/B)*L;
        c_t_EQ46(not,1,1:length(nmax))=(Cftinfty-C_f_t(not,1,1:length(nmax)))/L+ctinfty;
      end;
      % if nox==54; keyboard; end;% pour debug: mat2str2D=@(c)max(-999.9999,min(999.9999,c))+0*(1./(c.^2+1)), plus joli que mat2str pour une matrice 2D ...
    end
  end
  if length(t)*length(x)>1; c_x_t_EQ46_mat=[]; c_x_t_EQ46_x_mat=[]; end;
  return;
end

if ischar(x)
  % type sensca_shorttime('check_Int',0,0,0,0,0,0,0) to check various parts of this file.
  % type sensca_shorttime('check_Int',0,0,0,0,0,0,0) to check various parts of this file.
  c_x_t_EQ46=feval(x,t0,t,B,K0,Cftinfty,p,ctinfty);
  return
end
% calcul de Cft0 par conservation de la matiere
if abs(p)<eps^.8
  Cft0=Cftinfty-1;
else
  Cft0=Cftinfty-(1-exp(-p))/p;
end
% l'erreur est proportionnelle à erreEQ46.
u=t-t0;

racine=(B.^2.*K0.^2-4.*B).^(1./2);
sqrtr=[(-1./2.*K0.*B+1./2.*racine);(-1./2.*K0.*B-1./2.*racine)];
r=sqrtr.^2;

c1=B.^2.*K0./(r([2 1])-r);
c2=-B.*(r+B)./(r([2 1])-r);
c3=2.*K0.*B.*(r+B)./(r([2 1])-r);
c4=-2.*r.*B.^2.*K0.^2./(r([2 1])-r);
% note: erfc(x)~exp(-x^2)/x/sqrt(Pi)
if abs(p)>0 && x>-50
  x=x-100;% on utilise l'ancienne méthode, avec Int_exp_a_fois_erfc_b ... qui sait gérer p non nul.
end
if x>-50
  r=sqrtr;
%  mathtypetolatex2.09:
%   \frac{{c(x,{t_0} + {\rm{Fo}}) - {c^{t \to \infty }}}}{{{c^{{t_0}}} - {c^{t \to \infty }}}} =  \\
%   \frac{1}{2}\frac{{{C_f}({t_0}) - C_f^{t \to \infty }}}{{{c^{{t_0}}} - {c^{t \to \infty }}}}\sum\limits_{{n_5} = 1}^2 {\sum\limits_{{n_6} = 1}^2 {\sum\limits_{{n_4} = 1}^2 {\left( {{c_{1,{n_5}}} - {{( - 1)}^{{n_4}}}\frac{{{c_{2,{n_5}}}}}{{{r_{{n_5}}}}}} \right){e^{r_{{n_5}}^2{\rm{Fo}} + {{( - 1)}^{{n_4}}}\left( {1 + {{( - 1)}^{{n_6}}}x} \right){r_{{n_5}}}}}{\rm{erfc}}\left( {\frac{{1 + {{( - 1)}^{{n_6}}}x}}{{2\sqrt {{\rm{Fo}}} }} + {{( - 1)}^{{n_4}}}{r_{{n_5}}}\sqrt {{\rm{Fo}}} } \right)} } }  +  \\ 
%   \sum\limits_{{n_6} = 1}^3 {\sum\limits_{{n_7} = 1}^2 {\frac{{{{( - 1)}^{{n_7}}}}}{2}\frac{{2{n_6} - 5}}{{|2{n_6} - 5|}}{\rm{erfc}}\left( {\frac{{2{n_6} - 5}}{{|2{n_6} - 5|}}\frac{{1 - {{( - 1)}^{{n_7}}}}}{{4\sqrt {{\rm{Fo}}} }} + \frac{{2 + {{( - 1)}^{{n_6}}}x}}{{2\sqrt {{\rm{Fo}}} }}} \right)} }  +  \\
%   \sum\limits_{{n_5} = 1}^2 {\sum\limits_{{n_6} = 1}^3 {\sum\limits_{{n_4} = 1}^2 {\sum\limits_{{n_7} = 1}^4 {\left( {{{( - 1)}^{{n_4}}}\frac{{{c_{4,{n_5}}}}}{{{r_{{n_5}}}}} - {c_{3,{n_5}}}} \right)\frac{{{{( - 1)}^{{n_7}}}{{( - 1)}^{{n_4}}}}}{{4{r_{{n_5}}}}}\frac{{2{n_7} - 5}}{{|2{n_7} - 5|}}\frac{{2{n_6} - 5}}{{|2{n_6} - 5|}} \times } } } }  \\
%   {e^{r_{{n_5}}^2{\rm{Fo}} + {{( - 1)}^{{n_4}}}\left( {2 + {{( - 1)}^{{n_6}}}x} \right){r_{{n_5}}} + \frac{{1 + {{( - 1)}^{{n_7}}}}}{2}\left( {( - {{( - 1)}^{{n_4}}}{r_{{n_5}}})\left( {(2 + {{( - 1)}^{{n_6}}}x\sqrt {{\rm{Fo}}}  + 2{{( - 1)}^{{n_4}}}{\rm{Fo}}{r_{{n_5}}})} \right) + r_{{n_5}}^2} \right) + \frac{{1 - {{( - 1)}^{{n_7}}}}}{4}(1 - \frac{{3 - 2{n_7}}}{{|3 - 2{n_7}|}}){{( - 1)}^{{n_4}}}\frac{{2{n_6} - 5}}{{|2{n_6} - 5|}}{r_{{n_5}}}}} \times  \\
%   {\rm{erfc}}\left( {\frac{{2 + {{( - 1)}^{{n_6}}}x}}{{2\sqrt {{\rm{Fo}}} }} + (1 + \frac{{2{n_7} - 5}}{{|2{n_7} - 5|}})\frac{{2{n_6} - 5}}{{|2{n_6} - 5|}}\frac{1}{{4\sqrt {{\rm{Fo}}} }} + {{( - 1)}^{{n_4}}}\frac{{1 - {{( - 1)}^{{n_7}}}}}{2}{r_{{n_5}}}\sqrt {{\rm{Fo}}} } \right) +  \\
%   \sum\limits_{{n_7} = 1}^2 {\frac{{{{( - 1)}^{{n_7}}}}}{2}{\rm{erfc}}\left( {\frac{{1 - {{( - 1)}^{{n_7}}}}}{{4\sqrt {{\rm{Fo}}} }} - \frac{x}{{2\sqrt {{\rm{Fo}}} }}} \right)}  + \sum\limits_{{n_7} = 1}^2 {\frac{{{{( - 1)}^{{n_7}}}}}{2}{\rm{erfc}}\left( {\frac{{1 - {{( - 1)}^{{n_7}}}}}{{4\sqrt {{\rm{Fo}}} }} + \frac{x}{{2\sqrt {{\rm{Fo}}} }}} \right)}  + {e^{ - (\frac{{1 - \varepsilon }}{{{\rm{Fo}}}})}}O(1) \\
%  tex2maple:
%  '; EQ1 := [([[(c((x, t[0] + Fo)) - c[t,"=",infty])/((c@@())*(t[0]) - c[t,"=",infty]) =(1)/(2) * (C[f]((t[0])) - C[f][t,"=",infty])/((c@@())*(t[0]) - c[t,"=",infty]) * Sum(Sum(Sum((c[1, n[5]] - (-1)^(n[4]) * (c[2, n[5]])/(r[n[5]])) * exp(r[n[5]]^(2) * Fo + (-1)^(n[4]) * (1 + (-1)^(n[6]) * x) * r[n[5]]) * erfc(((1 + (-1)^(n[6]) * x)/(2 * sqrt(Fo)) + (-1)^(n[4]) * r[n[5]] * sqrt(Fo))), n[4]= 1..2), n[6]= 1..2), n[5]= 1..2)+Sum(Sum(((-1)^(n[7]))/(2) * (2 * n[6] - 5)/(abs(2 * n[6] - 5)) * erfc(((2 * n[6] - 5)/(abs(2 * n[6] - 5)) * (1 - (-1)^(n[7]))/(4 * sqrt(Fo)) + (2 + (-1)^(n[6]) * x)/(2 * sqrt(Fo)))), n[7]= 1..2), n[6]= 1..3)+Sum(Sum(Sum(Sum(((-1)^(n[4]) * (c[4, n[5]])/(r[n[5]]) - c[3, n[5]]) * ((-1)^(n[7]) * (-1)^(n[4]))/(4 * r[n[5]]) * (2 * n[7] - 5)/(abs(2 * n[7] - 5)) * (2 * n[6] - 5)/(abs(2 * n[6] - 5)) * ×, n[7]= 1..4), n[4]= 1..2), n[6]= 1..3), n[5]= 1..2)], [exp(r[n[5]]^(2) * Fo + (-1)^(n[4]) * (2 + (-1)^(n[6]) * x) * r[n[5]] + (1 + (-1)^(n[7]))/(2) * ((-(-1)^(n[4]) * r[n[5]]) * ((2 + (-1)^(n[6]) * x * sqrt(Fo) + 2 * (-1)^(n[4]) * Fo * r[n[5]])) + r[n[5]]^(2)) + (1 - (-1)^(n[7]))/(4) * (1 - (3 - 2 * n[7])/(abs(3 - 2 * n[7]))) * (-1)^(n[4]) * (2 * n[6] - 5)/(abs(2 * n[6] - 5)) * r[n[5]]) * erfc(((2 + (-1)^(n[6]) * x)/(2 * sqrt(Fo)) + (1 + (2 * n[7] - 5)/(abs(2 * n[7] - 5))) * (2 * n[6] - 5)/(abs(2 * n[6] - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n[4]) * (1 - (-1)^(n[7]))/(2) * r[n[5]] * sqrt(Fo)))+Sum(((-1)^(n[7]))/(2) * erfc(((1 - (-1)^(n[7]))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), n[7]= 1..2) + Sum(((-1)^(n[7]))/(2) * erfc(((1 - (-1)^(n[7]))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), n[7]= 1..2) + exp(-((1 - epsilon)/(F o))) * O((1))]])   ]');
%  '<,'>s:n\[\([1-9]\)]:n\1:gc
%  '<,'>s:c\[\([1-9]\), \([n1-9]*\)]:c\1(\2):gc
%  '<,'>s:\[:(:g
%  '<,'>s:]:):g
%  '<,'>s:\.\.:\::g
%  pour faire sum(@(n5), mettre dans le registre a: n%BhvBd%phr)bi@( 
  Fo=u;
  %c_x_t_EQ46 = ctinfty+ (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo))), 1:1),  1:2), 1:2)+Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * erfc_complex(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)))),  1:2),  1:3)+Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x * sqrt(Fo) + 2 * (-1)^(n4) * Fo * r(n5))) + r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))), 1:4), 1:1), 1:3), 1:2)+Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), 1:2) + Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), 1:2) ;
  % TODO: remplacer ceci par sum(c_x_t_EQ46_mat(:,1)), qui devrait aller plus vite en donnant la meme chose !!
  c_x_t_EQ46 = ctinfty+ (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo))), 1),  1:2), 1:2)+Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * erfc_complex(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)))),  1:2),  1:3)+Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))), 1:4), 1:1), 1:3), 1:2)+Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), 1:2) + Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), 1:2) ;
  epsilon=.01;% any positive number will do.
  errEQ46 = exp(-((1 - epsilon)./(u)));
  % save internals_sensca_l132
  if nargout<3
    return;
  end
  c_x_t_EQ46_x =  (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))]), 1:2),  1:2), 1:2)+Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * diff_erfc(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo))))*(( ( (-1)^(n6) )/(2 * sqrt(Fo)))),  1:2),  1:3)+Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b_ab(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)),[(-1)^(n4) * ( (-1)^(n6)) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * (( (-1)^(n6)))),(( (-1)^(n6))/(2 * sqrt(Fo)))]), 1:4), 1:2), 1:3), 1:2)+Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))*(( - (1)/(2 * sqrt(Fo)))), 1:2) + Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))))*(( (1)/(2 * sqrt(Fo)))), 1:2) ;
  if nargout==3
    return;
  end
  c_x_t_EQ46_mat = [ctinfty,0,0,ctinfty,0];% col 1: contribution additive, col 2: argument de exp, col 3: argument de erfc, col 4: coeff de proportionnalite, col 5: chaque chiffre est la valeur de l'indice de sommation.
  n5=kron([1],kron([1;1],[1;2]));
  n6=kron([1],kron([1;2],[1;1]));
  n4=kron([1],kron([1;1],[1;1]));
  c_x_t_EQ46_mat=[c_x_t_EQ46_mat;
    (1)./(2) .* (Cft0 - Cftinfty) .* (c1(n5) - (-1).^(n4) .* (c2(n5))./(r(n5))) .* exp_a_fois_erfc_b(r(n5).^(2) .* Fo + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* r(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (-1).^(n4) .* r(n5) .* sqrt(Fo))),r(n5).^(2) .* Fo + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* r(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (-1).^(n4) .* r(n5) .* sqrt(Fo)),(1)./(2) .* (Cft0 - Cftinfty) .* (c1(n5) - (-1).^(n4) .* (c2(n5))./(r(n5))),n5*100+n6*10+n4
  ];
  n6=kron([1;1],[1;2;3]);
  n7=kron([1;2],[1;1;1]);
  c_x_t_EQ46_mat=[c_x_t_EQ46_mat;
    ((-1).^(n7))./(2) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* erfc_complex(((2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)))),0*n6,((2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo))),((-1).^(n7))./(2) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)),n6*10+n7
  ];
  n5=kron([1;1;1;1],kron([1],kron([1;1;1],[1;2])));% 
  n6=kron([1;1;1;1],kron([1],kron([1;2;3],[1;1])));% 
  n4=kron([1;1;1;1],kron([1],kron([1;1;1],[1;1])));% 
  n7=kron([1;2;3;4],kron([1],kron([1;1;1],[1;1])));% 
  c_x_t_EQ46_mat=[c_x_t_EQ46_mat;
    ((-1).^(n4) .* (c4(n5))./(r(n5)) - c3(n5)) .* ((-1).^(n7) .* (-1).^(n4))./(4 .* r(n5)) .* (2 .* n7 - 5)./(abs(2 .* n7 - 5)) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* exp_a_fois_erfc_b(r(n5).^(2) .* Fo + (-1).^(n4) .* (2 + (-1).^(n6) .* x) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* ((2 + (-1).^(n6) .* x)) -Fo.* r(n5).^(2)) + (1 - (-1).^(n7))./(4) .* (1 - (3 - 2 .* n7)./(abs(3 - 2 .* n7))) .* (-1).^(n4) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* r(n5),((2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (1 + (2 .* n7 - 5)./(abs(2 .* n7 - 5))) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1)./(4 .* sqrt(Fo)) + (-1).^(n4) .* (1 - (-1).^(n7))./(2) .* r(n5) .* sqrt(Fo))),r(n5).^(2) .* Fo + (-1).^(n4) .* (2 + (-1).^(n6) .* x) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* ((2 + (-1).^(n6) .* x)) -Fo.* r(n5).^(2)) + (1 - (-1).^(n7))./(4) .* (1 - (3 - 2 .* n7)./(abs(3 - 2 .* n7))) .* (-1).^(n4) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* r(n5),((2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (1 + (2 .* n7 - 5)./(abs(2 .* n7 - 5))) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1)./(4 .* sqrt(Fo)) + (-1).^(n4) .* (1 - (-1).^(n7))./(2) .* r(n5) .* sqrt(Fo)),((-1).^(n4) .* (c4(n5))./(r(n5)) - c3(n5)) .* ((-1).^(n7) .* (-1).^(n4))./(4 .* r(n5)) .* (2 .* n7 - 5)./(abs(2 .* n7 - 5)) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)),n5*100+n6*10+n4+n7/10
  ];
  n7=[1;2];
  c_x_t_EQ46_mat=[c_x_t_EQ46_mat;
    ((-1).^(n7))./(2) .* erfc_complex(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) - (x)./(2 .* sqrt(Fo)))),0*n7,(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) - (x)./(2 .* sqrt(Fo)))),((-1).^(n7))./(2),n7+.1;
    ((-1).^(n7))./(2) .* erfc_complex(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (x)./(2 .* sqrt(Fo)))),0*n7,(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (x)./(2 .* sqrt(Fo)))),((-1).^(n7))./(2),n7+.2;
  ];
  if nargout==4
    return;
  end
  c_x_t_EQ46_x_mat = [0,0,0,0,0];% col 1: contribution additive, col 2: argument de exp, col 3: argument de erfc, col 4: coeff de proportionnalite, col 5: chaque chiffre est la valeur de l'indice de sommation.
  n5=kron([1],kron([1;1],[1;2]));
  n6=kron([1],kron([1;2],[1;1]));
  n4=kron([1],kron([1;1],[1;1]));
  c_x_t_EQ46_x_mat=[c_x_t_EQ46_x_mat;
    (1)./(2) .* (Cft0 - Cftinfty) .* (c1(n5) - (-1).^(n4) .* (c2(n5))./(r(n5))) .* exp_a_fois_erfc_b(r(n5).^(2) .* Fo + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* r(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (-1).^(n4) .* r(n5) .* sqrt(Fo))).*[ (-1).^(n4) .* ( (-1).^(n6)) .* r(n5)],r(n5).^(2) .* Fo + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* r(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (-1).^(n4) .* r(n5) .* sqrt(Fo)),(1)./(2) .* (Cft0 - Cftinfty) .* (c1(n5) - (-1).^(n4) .* (c2(n5))./(r(n5))).*[ (-1).^(n4) .* ( (-1).^(n6)) .* r(n5)],n5*100+n6*10+n4+.1;
    (1)./(2) .* (Cft0 - Cftinfty) .* (c1(n5) - (-1).^(n4) .* (c2(n5))./(r(n5))) .* (-2/sqrt(pi)).*moins_2_h2_exp_a_moins_carre_b_ab_sur_sqrt_pi(r(n5).^(2) .* Fo + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* r(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (-1).^(n4) .* r(n5) .* sqrt(Fo)),[ (-1).^(n4) .* ( (-1).^(n6)) .* r(n5),(( (-1).^(n6) )./(2 .* sqrt(Fo)))]),r(n5).^(2) .* Fo + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* r(n5)-((1 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (-1).^(n4) .* r(n5) .* sqrt(Fo)).^2,0*n5,(1)./(2) .* (Cft0 - Cftinfty) .* (c1(n5) - (-1).^(n4) .* (c2(n5))./(r(n5))) .* (-2/sqrt(pi)).*constante_moins_2_h2_exp_a_moins_carre_b_ab_sur_sqrt_pi(r(n5).^(2) .* Fo + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* r(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (-1).^(n4) .* r(n5) .* sqrt(Fo)),[ (-1).^(n4) .* ( (-1).^(n6)) .* r(n5),(( (-1).^(n6) )./(2 .* sqrt(Fo)))]),n5*100+n6*10+n4+.2
  ];
  n6=kron([1;1],[1;2;3]);
  n7=kron([1;2],[1;1;1]);
  c_x_t_EQ46_x_mat=[c_x_t_EQ46_x_mat;
    ((-1).^(n7))./(2) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* diff_erfc(((2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)))).*(( ( (-1).^(n6) )./(2 .* sqrt(Fo)))),-(((2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)))).^2,0*n7,((-1).^(n7))./(2) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* constante_diff_erfc(0).*(( ( (-1).^(n6) )./(2 .* sqrt(Fo)))),n6*10+n7
  ];
  n5=kron([1;1;1;1],kron([1],kron([1;1;1],[1;2])));% 
  n6=kron([1;1;1;1],kron([1],kron([1;2;3],[1;1])));% 
  n4=kron([1;1;1;1],kron([1],kron([1;1;1],[1;1])));% 
  n7=kron([1;2;3;4],kron([1],kron([1;1;1],[1;1])));% 
  c_x_t_EQ46_x_mat=[c_x_t_EQ46_x_mat;
    ((-1).^(n4) .* (c4(n5))./(r(n5)) - c3(n5)) .* ((-1).^(n7) .* (-1).^(n4))./(4 .* r(n5)) .* (2 .* n7 - 5)./(abs(2 .* n7 - 5)) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* exp_a_fois_erfc_b(r(n5).^(2) .* Fo + (-1).^(n4) .* (2 + (-1).^(n6) .* x) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* ((2 + (-1).^(n6) .* x)) -Fo.* r(n5).^(2)) + (1 - (-1).^(n7))./(4) .* (1 - (3 - 2 .* n7)./(abs(3 - 2 .* n7))) .* (-1).^(n4) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* r(n5),((2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (1 + (2 .* n7 - 5)./(abs(2 .* n7 - 5))) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1)./(4 .* sqrt(Fo)) + (-1).^(n4) .* (1 - (-1).^(n7))./(2) .* r(n5) .* sqrt(Fo))).*[(-1).^(n4) .* ( (-1).^(n6)) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* (( (-1).^(n6))))],r(n5).^(2) .* Fo + (-1).^(n4) .* (2 + (-1).^(n6) .* x) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* ((2 + (-1).^(n6) .* x)) -Fo.* r(n5).^(2)) + (1 - (-1).^(n7))./(4) .* (1 - (3 - 2 .* n7)./(abs(3 - 2 .* n7))) .* (-1).^(n4) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* r(n5),((2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (1 + (2 .* n7 - 5)./(abs(2 .* n7 - 5))) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1)./(4 .* sqrt(Fo)) + (-1).^(n4) .* (1 - (-1).^(n7))./(2) .* r(n5) .* sqrt(Fo)),((-1).^(n4) .* (c4(n5))./(r(n5)) - c3(n5)) .* ((-1).^(n7) .* (-1).^(n4))./(4 .* r(n5)) .* (2 .* n7 - 5)./(abs(2 .* n7 - 5)) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)).*[(-1).^(n4) .* ( (-1).^(n6)) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* (( (-1).^(n6))))],n5*100+n6*10+n4+n7/10+0.01;
    ((-1).^(n4) .* (c4(n5))./(r(n5)) - c3(n5)) .* ((-1).^(n7) .* (-1).^(n4))./(4 .* r(n5)) .* (2 .* n7 - 5)./(abs(2 .* n7 - 5)) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* moins_2_h2_exp_a_moins_carre_b_ab_sur_sqrt_pi(r(n5).^(2) .* Fo + (-1).^(n4) .* (2 + (-1).^(n6) .* x) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* ((2 + (-1).^(n6) .* x)) -Fo.* r(n5).^(2)) + (1 - (-1).^(n7))./(4) .* (1 - (3 - 2 .* n7)./(abs(3 - 2 .* n7))) .* (-1).^(n4) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* r(n5),((2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (1 + (2 .* n7 - 5)./(abs(2 .* n7 - 5))) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1)./(4 .* sqrt(Fo)) + (-1).^(n4) .* (1 - (-1).^(n7))./(2) .* r(n5) .* sqrt(Fo)),[(-1).^(n4) .* ( (-1).^(n6)) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* (( (-1).^(n6)))),(( (-1).^(n6))./(2 .* sqrt(Fo)))]),r(n5).^(2) .* Fo + (-1).^(n4) .* (2 + (-1).^(n6) .* x) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* ((2 + (-1).^(n6) .* x)) -Fo.* r(n5).^(2)) + (1 - (-1).^(n7))./(4) .* (1 - (3 - 2 .* n7)./(abs(3 - 2 .* n7))) .* (-1).^(n4) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* r(n5)-((2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (1 + (2 .* n7 - 5)./(abs(2 .* n7 - 5))) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1)./(4 .* sqrt(Fo)) + (-1).^(n4) .* (1 - (-1).^(n7))./(2) .* r(n5) .* sqrt(Fo)).^2,0*n7,((-1).^(n4) .* (c4(n5))./(r(n5)) - c3(n5)) .* ((-1).^(n7) .* (-1).^(n4))./(4 .* r(n5)) .* (2 .* n7 - 5)./(abs(2 .* n7 - 5)) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* constante_moins_2_h2_exp_a_moins_carre_b_ab_sur_sqrt_pi(r(n5).^(2) .* Fo + (-1).^(n4) .* (2 + (-1).^(n6) .* x) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* ((2 + (-1).^(n6) .* x)) -Fo.* r(n5).^(2)) + (1 - (-1).^(n7))./(4) .* (1 - (3 - 2 .* n7)./(abs(3 - 2 .* n7))) .* (-1).^(n4) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* r(n5),((2 + (-1).^(n6) .* x)./(2 .* sqrt(Fo)) + (1 + (2 .* n7 - 5)./(abs(2 .* n7 - 5))) .* (2 .* n6 - 5)./(abs(2 .* n6 - 5)) .* (1)./(4 .* sqrt(Fo)) + (-1).^(n4) .* (1 - (-1).^(n7))./(2) .* r(n5) .* sqrt(Fo)),[(-1).^(n4) .* ( (-1).^(n6)) .* r(n5) + (1 + (-1).^(n7))./(2) .* ((-(-1).^(n4) .* r(n5)) .* (( (-1).^(n6)))),(( (-1).^(n6))./(2 .* sqrt(Fo)))]),n5*100+n6*10+n4+n7/10+0.02
  ];
  n7=[1;2];
  c_x_t_EQ46_x_mat=[c_x_t_EQ46_x_mat;
    ((-1).^(n7))./(2) .* diff_erfc(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) - (x)./(2 .* sqrt(Fo)))).*(( - (1)./(2 .* sqrt(Fo)))),-(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) - (x)./(2 .* sqrt(Fo)))).^2,0*n7,((-1).^(n7))./(2) .* constante_diff_erfc(0).*(( - (1)./(2 .* sqrt(Fo)))),n7*10+1;
    ((-1).^(n7))./(2) .* diff_erfc(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (x)./(2 .* sqrt(Fo)))).*(( (1)./(2 .* sqrt(Fo)))) ,-(((1 - (-1).^(n7))./(4 .* sqrt(Fo)) + (x)./(2 .* sqrt(Fo)))).^2,0*n7,((-1).^(n7))./(2) .* constante_diff_erfc(0).*(( (1)./(2 .* sqrt(Fo)))),n7*10+2
  ];
  if nargout==5
    return;
  end
  c_x_t_EQ46_t = (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[r(n5)^(2),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * r(n5) / sqrt(Fo)*(1/2))]), 1:1),  1:2), 1:2)+Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * diff_erfc(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo))))*(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)*Fo)*(-1/2) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2))),  1:2),  1:3)+Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b_ab(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)),[r(n5)^(2) + (1 + (-1)^(n7))/(2) * ( - r(n5)^(2)),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) / sqrt(Fo)*(1/2))]), 1:4), 1:1), 1:3), 1:2)+Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))*(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))/Fo*(-1/2), 1:2) + Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))))*(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))))/Fo*(-1/2), 1:2) ;
  c_x_t_EQ46_x_t =  (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[0,(((-1)^(n6) )/(2 * sqrt(Fo)*Fo)*(-1/2))],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))],[r(n5)^(2),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * r(n5) / sqrt(Fo)*(1/2))]), 1:1),  1:2), 1:2)+Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * diff_erfc_plus_diff_diff_erfc(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo))),((((-1)^(n6) )/(2 * sqrt(Fo)*Fo)*(-1/2))),( ( (-1)^(n6) )/(2 * sqrt(Fo))),(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)*Fo)*(-1/2) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2)))),  1:2),  1:3)+Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)),[0,(( (-1)^(n6) )/(2 * sqrt(Fo)*Fo)*(-1/2))],[(-1)^(n4) * ( (-1)^(n6)) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * (( (-1)^(n6)))),(( (-1)^(n6))/(2 * sqrt(Fo)))],[r(n5)^(2) + (1 + (-1)^(n7))/(2) * ( - r(n5)^(2)),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) / sqrt(Fo)*(1/2))]), 1:4), 1:1), 1:3), 1:2)+Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))),(( (1)/(2 * sqrt(Fo))))/Fo*(-1/2),( - (1)/(2 * sqrt(Fo)))*(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))/Fo*(-1/2),(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))/Fo*(-1/2)), 1:2) + Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))),(((1)/(2 * sqrt(Fo)))),( (1)/(2 * sqrt(Fo))),(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))))/Fo*(-1/2)), 1:2) ; % bien égal à c_x_t_EQ46_t, OUF !
  c_x_t_EQ46_x_x =  (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[0,0],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))]), 1:1),  1:2), 1:2)+Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * diff_erfc_plus_diff_diff_erfc(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo))),0,( ( (-1)^(n6) )/(2 * sqrt(Fo))),( ( (-1)^(n6) )/(2 * sqrt(Fo)))),  1:2),  1:3)+Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)),[0,0],[(-1)^(n4) * ( (-1)^(n6)) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * (( (-1)^(n6)))),(( (-1)^(n6))/(2 * sqrt(Fo)))],[(-1)^(n4) * ( (-1)^(n6)) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * (( (-1)^(n6)))),(( (-1)^(n6))/(2 * sqrt(Fo)))]), 1:4), 1:1), 1:3), 1:2)+Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))),0,( - (1)/(2 * sqrt(Fo))),( - (1)/(2 * sqrt(Fo)))), 1:2) + Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))),0,( (1)/(2 * sqrt(Fo))),( (1)/(2 * sqrt(Fo)))), 1:2) ; % bien égal à c_x_t_EQ46_t, OUF !
  c_x_t_EQ46_x_t_t = NaN+ (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[0,(((-1)^(n6) )/(2 * sqrt(Fo)*Fo)*(-1/2))],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))],[r(n5)^(2),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * r(n5) / sqrt(Fo)*(1/2))]), 1:1),  1:2), 1:2)+Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * diff_erfc_plus_diff_diff_erfc(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo))),((((-1)^(n6) )/(2 * sqrt(Fo)*Fo)*(-1/2))),( ( (-1)^(n6) )/(2 * sqrt(Fo))),(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)*Fo)*(-1/2) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2)))),  1:2),  1:3)+Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)),[0,(( (-1)^(n6) )/(2 * sqrt(Fo)*Fo)*(-1/2))],[(-1)^(n4) * ( (-1)^(n6)) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * (( (-1)^(n6)))),(( (-1)^(n6))/(2 * sqrt(Fo)))],[r(n5)^(2) + (1 + (-1)^(n7))/(2) * ( - r(n5)^(2)),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) / sqrt(Fo)*(1/2))]), 1:4), 1:1), 1:3), 1:2)+Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))),(( (1)/(2 * sqrt(Fo))))/Fo*(-1/2),( - (1)/(2 * sqrt(Fo)))*(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))/Fo*(-1/2),(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))/Fo*(-1/2)), 1:2) + Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))),(((1)/(2 * sqrt(Fo)))),( (1)/(2 * sqrt(Fo))),(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))))/Fo*(-1/2)), 1:2) ; % % TODO a changer ... tout deriver en Fo.
  if x==1
    % - c_x_t_EQ46_x/B=K0*c_x_t_EQ46 - C_f_t
    C_f_t=K0*c_x_t_EQ46+ c_x_t_EQ46_x/B;
    c_t_EQ46=-C_f_t+Cftinfty+ctinfty;
    % maple value(subs(EQ15[1]=EQ15[3],c(y,t[0])=1+c[t,"=",infty],op([1,1,3],EQ19)))
    % -c_x_t_EQ46_x_t/B=K0*c_x_t_EQ46_t+c_x_t_EQ46_x selon eq 7 de adsorptionv57.tex
    % -s*c_s_t_x/B+c_s_x_t0/B=K0*s*c_x_t-K0*c_x_t0+K0*c_x_tinfty+c_x_t_EQ46_x
    [-c_x_t_EQ46_x_t/B, K0*c_x_t_EQ46_t+c_x_t_EQ46_x];
    % mat2str([1.054,ans,diff(ans),log10(max(abs(ans))),log10(abs(diff(ans)))])% vérification dérivée en temps de la condition de bord 
    % pour importer en latex le mathtype stocké dans les .doc et exporté en .tex par ooffice:
    %{
    '<,'>s:{\\textbackslash}:\\:g
    '<,'>s:\\left::g
    '<,'>s:\\right::g
    '<,'>s:\\\([{}]\):\1:g
    
    \begin{array}{*{20}c}\par \ \ \ { {}-
    \frac{1}{Bi}\frac{{\partial
    \^{}2 c}}{{\partial x\partial
    t}}(1,t) =
    \frac{{\rm{d}}}{{{\rm{d}}t}}(
    {\beta \_f [ {\beta
    \_s\^{}{ {}- 1} ( {c(
    {1,t} )} )}
    ]} ) +
    L\frac{{\partial
    c}}{{\partial x}}(1,t)} \&
    {{\rm{for}}} \& {t {\textgreater} t\_0 }
    \ \\\par \end{array}
    \par \]}
    %}
  end
  % disp(mat2str([1.053,c_x_t_EQ46_x_x,c_x_t_EQ46_t,c_x_t_EQ46_x_x-c_x_t_EQ46_t,log10(abs(c_x_t_EQ46_x_x)),log10(abs(c_x_t_EQ46_x_x-c_x_t_EQ46_t))]))% verification eq diffusion OK.
  % bon, second terme des deux vecteurs ci-dessous ne sont pas les meme, alors que le test ver103..ver106 affirme qu'ils devraient l'etre. On compare leurs calculs.
  % res=erfc_complex(1./2.*b)./p-erfc_complex(1./2.*a.*yspan+1./2.*b)./p./exp(p.*yspan)+exp(p.*b./a+p.^2./a.^2).*erfc_complex(1./2.*a.*yspan+1./2.*b+p./a)-exp(p.*b./a+p.^2./a.^2).*erfc_complex(1./2.*b+p./a))./p; fait des NaN de temps en temps.
  % ((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (3 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)))
  %mat2str([1.001,(1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo))), 1:2),  1:2), 1:2),Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * erfc_complex(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)))),  1:2),  1:3),Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))), 1:4), 1:2), 1:3), 1:2),Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), 1:2) , Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), 1:2),(1)/(2) * (Cft0 - Cftinfty)] )
  if verLessThan('symbolic','3.2.3')% spécial maple dans matlab2007b.
    save totonew
    maple('restart');
    maple interface(prettyprint=1)                                                                                                                          
    maple('exp_a_fois_erfc_b:=(a,b)->exp(a)*erfc(b);');%parse("proc" "(a,b)local res;res:=exp(a)*erfc(b);print([1.005,res]);res;en" "d proc");')
    maple('erfc_complex:=erfc;')
    maple('Somme:=(funct,l)->print([op(op(3,funct))]);');%sum(op(2,funct),op(1,funct)=l);');
    maple('c_x_t_EQ46_t_form := diff(ctinfty+ (1)/(2) * (Cft0 - Cftinfty) * Somme( n5->Somme(n6->Somme(n4->(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo))), 1..1),  1..2), 1..2)+Somme(n6->Somme(n7->((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * erfc_complex(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)))),  1..2),  1..3)+Somme(n5-> Somme(n6-> Somme(n4-> Somme(n7-> ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))), 1..4), 1..1), 1..3), 1..2)+Somme(n7-> ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), 1..2) + Somme(n7-> ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), 1..2),u) ;')
    maple('c_x_t_EQ46_x_form := diff(ctinfty+ (1)/(2) * (Cft0 - Cftinfty) * Somme( n5->Somme(n6->Somme(n4->(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo))), 1..1),  1..2), 1..2)+Somme(n6->Somme(n7->((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * erfc_complex(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)))),  1..2),  1..3)+Somme(n5-> Somme(n6-> Somme(n4-> Somme(n7-> ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))), 1..4), 1..1), 1..3), 1..2)+Somme(n7-> ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), 1..2) + Somme(n7-> ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), 1..2),x,x) ;')
% a commenter une fois le debug fini.
    a=load('totonew')
    for f=fieldnames(a)'% load pour maple...
    v=a.(f{1});if isfloat(v);if length(v)==1;s=[f{1},':=',mat2str(v)];maple(s);else;for i=1:length(v);s=[f{1},'(',num2str(i),'):=',mat2str(v(i))];;maple(s);end;end;end;
    end
    maple('exp_a_fois_erfc_b:=parse("proc" "(a,b)local res;res:=exp(a)*erfc(b);print([1.005,res]);res;en" "d proc");')
    maple('Somme:=parse("proc" "(funct,l)local c,res,tmp,bres,k;res:=0;bres:=[seq(k,k=l)];print(bres);fo" "r c from op(1,l) to op(2,l) do;tmp:=evalf(funct(c));if nops([tmp])<>1 then;print([1.00004,funct,c,l]);fi;bres[c]:=tmp;res:=res+bres[c];od;print([1.006,op(2,l),res,bres,funct]);print(bres);res;en" "d proc");Somme(a->cos(a+1),[1,2]);')
  % :.s:@(\(n[0-9]\)):\1->:g
  % :.s/:/../g
  % :.s/Sum/Somme/g
    maple('c_x_t_EQ46 = ctinfty+ (1)/(2) * (Cft0 - Cftinfty) * Somme( n5->Somme(n6->Somme(n4->(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo))), 1..1),  1..2), 1..2)+Somme(n6->Somme(n7->((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * erfc_complex(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)))),  1..2),  1..3)+Somme(n5-> Somme(n6-> Somme(n4-> Somme(n7-> ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))), 1..4), 1..1), 1..3), 1..2)+Somme(n7-> ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), 1..2) + Somme(n7-> ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), 1..2) ;')
    mat2str([1.0531, (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[r(n5)^(2),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * r(n5) / sqrt(Fo)*(1/2))]), 1:1),  1:2), 1:2),Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * diff_erfc(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo))))*(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)*Fo)*(-1/2) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2))),  1:2),  1:3),Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b_ab(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)),[r(n5)^(2) + (1 + (-1)^(n7))/(2) * ( - r(n5)^(2)),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) / sqrt(Fo)*(1/2))]), 1:4), 1:1), 1:3), 1:2),Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))*(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))))/Fo*(-1/2), 1:2) , Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))))*(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))))/Fo*(-1/2), 1:2)])
    mat2str([1.0532, (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[0,0],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))]), 1:1),  1:2), 1:2),Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * diff_erfc_plus_diff_diff_erfc(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo))),0,( ( (-1)^(n6) )/(2 * sqrt(Fo))),( ( (-1)^(n6) )/(2 * sqrt(Fo)))),  1:2),  1:3),Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo)),[0,0],[(-1)^(n4) * ( (-1)^(n6)) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * (( (-1)^(n6)))),(( (-1)^(n6))/(2 * sqrt(Fo)))],[(-1)^(n4) * ( (-1)^(n6)) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * (( (-1)^(n6)))),(( (-1)^(n6))/(2 * sqrt(Fo)))]), 1:4), 1:1), 1:3), 1:2),Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo))),0,( - (1)/(2 * sqrt(Fo))),( - (1)/(2 * sqrt(Fo)))), 1:2) , Sum(@(n7) ((-1)^(n7))/(2) * diff_erfc_plus_diff_diff_erfc(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo))),0,( (1)/(2 * sqrt(Fo))),( (1)/(2 * sqrt(Fo)))), 1:2)])
  
    mat2str([1.001,ctinfty, (1)/(2) * (Cft0 - Cftinfty) * Sum( @(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo))), 1:1),  1:2), 1:2),Sum(@(n6)Sum(@(n7)((-1)^(n7))/(2) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * erfc_complex(((2 * n6 - 5)/(abs(2 * n6 - 5)) * (1 - (-1)^(n7))/(4 * sqrt(Fo)) + (2 + (-1)^(n6) * x)/(2 * sqrt(Fo)))),  1:2),  1:3),Sum(@(n5) Sum(@(n6) Sum(@(n4) Sum(@(n7) ((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * exp_a_fois_erfc_b(r(n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5),((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))), 1:4), 1:1), 1:3), 1:2),Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) - (x)/(2 * sqrt(Fo)))), 1:2) , Sum(@(n7) ((-1)^(n7))/(2) * erfc_complex(((1 - (-1)^(n7))/(4 * sqrt(Fo)) + (x)/(2 * sqrt(Fo)))), 1:2),1.0011, (1)/(2) * (Cft0 - Cftinfty) ])
    maple('a:=(2*n6-5)/abs(2*n6-5)/sqrt(Fo);')
    maple('b:=(2+(-1)^n6*x)/sqrt(Fo)+2*(-1)^n4*r(n5)*sqrt(Fo);')
    maple('p:=-(-1)^n4*(2*n6-5)/abs(2*n6-5)*r(n5);')
    maple('tmp1:=factor(subs(r(n5)=r(n5)^2,sqrtr(n5)=r(n5),u=Fo,r(n5) * u + (-1)^(n4) * (2 + (-1)^(n6) * x + (2*n6 - 5)/(abs(2*n6 - 5)) * y) * sqrtr(n5)) + p*y);')
    maple('ver2:=expand(subs(r(n5)=r(n5)^2,sqrtr(n5)=r(n5),u=Fo, ((((2 + (-1)^(n6) * x + (2*n6 - 5)/(abs(2*n6 - 5)) * y))/(2 * sqrt(u)) + (-1)^(n4) * sqrtr(n5)*sqrt(u))))-(a*y+b)/2);')
  % '; EQ48 := [([[Int((erfc(((a * y + b)/(2))))/(exp(p * y)), y= 0..y[prime]),"=",((([[(erfc(((b)/(2))))/(p) - (erfc(((a * y[prime] + b)/(2))))/(p * exp(p * y[prime])) + exp((p * b)/(a) + (p^(2))/(a^(2))) * (erfc(((a * y[prime] + b)/(2) + (p)/(a))) - erfc(((b)/(2) + (p)/(a))))/(p), p ,"<>", 0], [(a * y[prime] + b)/(a) * erfc(((a * y[prime] + b)/(2))) - (b)/(a) * erfc(((b)/(2))) + 2 * (exp(-b^(2)/4) - exp(-(a * y[prime] + b)^(2)/4))/(a * sqrt(Pi)), p,"=",0]]))),"=",Sum(((-1)^(n[7]))/(p) * (2*n[7] - 5)/(abs(2*n[7] - 5)) * exp((1 + (-1)^(n[7]))/(2) * (p * b * a + p^(2))/(a^(2)) - (1 - (-1)^(n[7]))/(2) * (1 - (3 - 2*n[7])/(abs(3 - 2*n[7]))) * (p * y[prime])/(2)) * erfc(((b)/(2) + (1 + (-1)^(n[7]))/(2) * (p)/(a) + (1 + (2*n[7] - 5)/(abs(2*n[7] - 5))) * (a * y[prime])/(4))), n[7]= 1..4)]]) ]'); % de senscnonlin_v68.m
  % '; EQ48 := [([[Int((erfc(((a * y + b)/(2))))/(exp(p * y)), y= 0..y[prime]),"=",((([[(erfc(((b)/(2))))/(p) - (erfc(((a * y[prime] + b)/(2))))/(p * exp(p * y[prime])) + exp((p * b)/(a) + (p^(2))/(a^(2))) * (erfc(((a * y[prime] + b)/(2) + (p)/(a))) - erfc(((b)/(2) + (p)/(a))))/(p), p ,"<>", 0], [(a * y[prime] + b)/(a) * erfc(((a * y[prime] + b)/(2))) - (b)/(a) * erfc(((b)/(2))) + 2 * (exp(-b^(2)/4) - exp(-(a * y[prime] + b)^(2)/4))/(a * sqrt(Pi)), p,"=",0]]))),"=",Sum(((-1)^(n7))/(p) * (2*n7 - 5)/(abs(2*n7 - 5)) * exp((1 + (-1)^(n7))/(2) * (p * b * a + p^(2))/(a^(2)) - (1 - (-1)^(n7))/(2) * (1 - (3 - 2*n7)/(abs(3 - 2*n7))) * (p * y[prime])/(2)) * erfc(((b)/(2) + (1 + (-1)^(n7))/(2) * (p)/(a) + (1 + (2*n7 - 5)/(abs(2*n7 - 5))) * (a * y[prime])/(4))), n7= 1..4)]]) ]'); % puis y[prime] devient 1
    maple('tmp3:=simplify((n5)^(2) * Fo + (-1)^(n4) * (2 + (-1)^(n6) * x) * r(n5) + (1 + (-1)^(n7))/(2) * ((-(-1)^(n4) * r(n5)) * ((2 + (-1)^(n6) * x)) -Fo* r(n5)^(2)) + (1 - (-1)^(n7))/(4) * (1 - (3 - 2 * n7)/(abs(3 - 2 * n7))) * (-1)^(n4) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * r(n5) - ((1 + (-1)^(n7))/(2) * (p * b * a + p^(2))/(a^(2)) - (1 - (-1)^(n7))/(2) * (1 - (3 - 2*n7)/(abs(3 - 2*n7))) * (p * 1)/(2))-tmp1);')
    maple combine(subs(n4=1,ver3))
    maple subs(n4=2,ver3) % OK!
    maple('ver4:=factor(((2 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (1 + (2 * n7 - 5)/(abs(2 * n7 - 5))) * (2 * n6 - 5)/(abs(2 * n6 - 5)) * (1)/(4 * sqrt(Fo)) + (-1)^(n4) * (1 - (-1)^(n7))/(2) * r(n5) * sqrt(Fo))-(((b)/(2) + (1 + (-1)^(n7))/(2) * (p)/(a) + (1 + (2*n7 - 5)/(abs(2*n7 - 5))) * (a * 1)/(4))));')
    maple('ver5:=factor(((-1)^(n7))/(p) * (2*n7 - 5)/(abs(2*n7 - 5))*subs(r(n5)=r(n5)^2,sqrtr(n5)=r(n5), 1/(4) * (c3(n5) - (-1)^(n4) * (c4(n5))/(sqrtr(n5))) ) -(((-1)^(n4) * (c4(n5))/(r(n5)) - c3(n5)) * ((-1)^(n7) * (-1)^(n4))/(4 * r(n5)) * (2 * n7 - 5)/(abs(2 * n7 - 5)) * (2 * n6 - 5)/(abs(2 * n6 - 5))));')
    maple eval(subs(n4=2,n6=1,ver5))
    maple eval(subs(n4=1,n6=1,ver5))
    maple eval(subs(n4=2,n6=3,ver5))
    maple eval(subs(n4=1,n6=3,ver5))
    % on regarde si EQ36[1][1][3] verifie l'équation de diffusion
    maple restart
    senscnonlin_v68
    maple tmp101:=value(subs(sqrt(z/u)=x/2/sqrt(u),sqrt(r*z)=sqrt(r)*x/2,sqrt(r*u)=sqrt(r)*sqrt(u),z=x^2/4,EQ36[1][1][3]))
    maple ver101:=simplify(diff(tmp101,x,x))-simplify(diff(tmp101,u))
    maple tmp102:=(value(subs(r[n[5]]=r,c[1,n[5]]=1,c[2,n[5]]=sqrt(r),n[6]=n,sqrt(r*u)=sqrt(r)*sqrt(u),op([1,1,3,2,-1,1,1],EQ46))))
    maple tmp105:=(value(subs(r[n[5]]=r,c[1,n[5]]=1,c[2,n[5]]=-sqrt(r),n[6]=n,sqrt(r*u)=sqrt(r)*sqrt(u),op([1,1,3,2,-1,1,1],EQ46))))
    maple ver103:=factor(subs(n=1,expand(diff(tmp102,u)-diff(tmp102,x,x))))
    maple ver104:=factor(subs(n=2,expand(diff(tmp102,u)-diff(tmp102,x,x))))
    maple ver106:=factor(subs(n=1,expand(diff(tmp105,u)-diff(tmp105,x,x))))
    maple ver107:=factor(subs(n=2,expand(diff(tmp105,u)-diff(tmp105,x,x))))
    % BUG le premier terme de la dérivée en t de la somme est différent, -5592.00867264452+i*7270.89141715664 chez matlab et 5592.0648421664349935710695260054+7270.6738735336845567871600567306*i chez maple ... trop de différence. tout ca à cause d'une confusion entre ' et .' ... :-(
    maple tmp109:=[(value(subs(Sum=seq,sqrt(r[n[5]]*u)=r(n[5])*sqrt(u),sqrt(r[n[5]])=r(n[5]),(r[n[5]])^(-1/2)=1/r(n[5]),c[1,n[5]]=c1(n[5]),c[2,n[5]]=c2(n[5]),r[n[5]]=r(n[5])^2,op([1,1,3,2,-1,1,1,1],EQ46))))]
    maple tmp110:=diff(tmp109,u);
    maple tmp111:=diff(tmp109,x,x);
    % vérification que EQ19 vérifie bien les conditions de bord:
    maple c_x_s:=subs(K[0]=K0,piecewise(x<1,0,1)=0,piecewise(-x<0,0,1)=0,value(subs(C[f](t[0])=Cft0,C[f][t,"=",infty]=Cftinfty,oop(0,EQ16[1])=unapply(EQ16[3],q),EQ15[1]=EQ15[5],c(y,t[0])=1+c[t,"=",infty],op([1,1,3],EQ19))))
    maple c_x_s_tout:=subs(K[0]=K0,piecewise(x<1,0,1)=0,piecewise(-x<0,0,1)=0,value(subs(C[f](t[0])=Cft0,C[f][t,"=",infty]=Cftinfty,oop(0,EQ16[1])=unapply(EQ16[3],q),EQ15[1]=EQ15[5],c(y,t[0])=1+c[t,"=",infty],op([1,1,3],EQ19))))
    maple c_x_s_1:=subs(K[0]=K0,piecewise(x<1,0,1)=0,piecewise(-x<0,0,1)=0,eval(value(subs(C[f](t[0])=Cft0,C[f][t,"=",infty]=Cftinfty,oop(0,EQ16[1])=unapply(EQ16[3],q),EQ15[1]=EQ15[5],c(y,t[0])=1+c[t,"=",infty],op([1,1,3,1],EQ19)))))
    maple c_x_s_2:=subs(K[0]=K0,piecewise(x<1,0,1)=0,piecewise(-x<0,0,1)=0,value(subs(op([1,1,3,2,1,2],EQ19)=op([1,1,3,2,1,2,1],EQ19),C[f](t[0])=Cft0,C[f][t,"=",infty]=Cftinfty,oop(0,EQ16[1])=unapply(EQ16[3],q),EQ15[1]=EQ15[5],c(y,t[0])=1+c[t,"=",infty],op([1,1,3,2],EQ19))))
    maple c_x_s_3:=subs(K[0]=K0,piecewise(x<1,0,1)=0,piecewise(-x<0,0,1)=0,value(subs(op([1,1,3,2,1,2],EQ19)=op([1,1,3,2,1,2,2],EQ19),C[f](t[0])=Cft0,C[f][t,"=",infty]=Cftinfty,oop(0,EQ16[1])=unapply(EQ16[3],q),EQ15[1]=EQ15[5],c(y,t[0])=1+c[t,"=",infty],op([1,1,3,2],EQ19))))
    maple c_x_s_4:=subs(K[0]=K0,piecewise(x<1,0,1)=0,piecewise(-x<0,0,1)=0,value(subs(op([1,1,3,2,1,2],EQ19)=op([1,1,3,2,1,2,3],EQ19),C[f](t[0])=Cft0,C[f][t,"=",infty]=Cftinfty,oop(0,EQ16[1])=unapply(EQ16[3],q),EQ15[1]=EQ15[5],c(y,t[0])=1+c[t,"=",infty],op([1,1,3,2],EQ19))))
    maple tmp222:=expand(factor(c_x_s-c_x_s_1))
    maple tmp224:=1/(e[3](q)*exp(q)-e[3](-q)/exp(q))=1/(e[3](q)*exp(q)^2-e[3](-q))*exp(q)       
    maple  ver223:=subs(tmp224,expand(factor(tmp222)))-expand(factor(c_x_s_2+c_x_s_3+c_x_s_4))
    maple c_x_s_x:=diff(c_x_s,x)
    maple c_s:=int(c_x_s,x=0..1)
    maple ver141:=subs(x=0,c_x_s_x)% cond de gauche juste.
    maple Cft0:=Cftinfty-1;
    maple K[0]:=K0;
    maple tmp142:=simplify(eval(subs(x=1,op(0,EQ16[1])=unapply(EQ16[3],q),-c_x_s_x/B)))
    maple tmp143:=simplify(eval(subs(x=1,op(0,EQ16[1])=unapply(EQ16[3],q),-(K0*c_x_s))))
    maple tmp144:=simplify(eval(subs(x=1,op(0,EQ16[1])=unapply(EQ16[3],q),-(c_s))))
    maple ver145:=combine(1/op(-1,tmp144)-expand(exp(q)/op(-1,tmp143)))% nul car meme denominateurs.
    % TODO: calculer la transfor inverse de laplace à la main pour notre valeur de t.
    %maple c_x_t_marchepas:=1/(2*Pi)*int(c_x_s(s=c+1i*u)*exp((c+1i*u)*t),u=-infinity..infinity)
    maple int(1/(c+i*z)*exp((c+i*z)*t),z=-infinity..infinity)/2/Pi
    maple('Undef:=parse("proc" "(expr,s);infinity/infinity;en" "d proc");')
    maple('invl:=parse("proc" "(expr,s,t,c,D,T)local res,a,nT,nC,Undef,debug;debug:=0;Undef:=(a,b)->infinity/infinity;res:=[];fo" "r Digits from D+debug to D+1 do;fo" "r nC from 1+debug to 2 do;fo" "r nT from 1+debug to 2 do;res:=[op(res),evalf(Int(subs(s=c*nC+i*z,expr)*exp((c*nC+i*z)*t),z=-T*nT..T*nT)/2/Pi)];print(eval(subs(Int=Undef,op(-1,res))));od;od;od;eval(subs(Int=Undef,[op(res),`D`=D,seq(op(a+4,res)-op(a,res),a=1..4-10*debug),`T`=T,seq(op(a*2,res)-op(a*2-1,res),a=1..4-10*debug),`c`=c,seq(op(a+2,res)-op(a,res),a=1..2-10*debug),seq(op(a+2,res)-op(a,res),a=5..6)]));en" "d proc");')% selon formule (2.18) p19 de Levine??theControl
    % maple invl(100/(s^2+1e4),s,1e-3,2,18,1e5)% marche ! donne bien sin(.1)
    % maple invl(2/(s^2+4),s,1e-3,2,16,1e4)% marche ! donne bien cos(2e-3)
    maple('Somme:=parse("proc" "(funct,l)local c,res,tmp,bres,k;res:=0;bres:=[seq(k,k=l)];print(bres);fo" "r c from op(1,l) to op(2,l) do;tmp:=evalf(funct(c));if nops([tmp])<>1 then;print([1.00004,funct,c,l]);fi;bres[c]:=tmp;res:=res+bres[c];od;print([1.006,op(2,l),res,bres,funct]);print(bres);res;en" "d proc");Somme(a->cos(a+1),[1,2]);')
    %maple('_EnvAllSolutions:=evalb(1=1);int(1/(12+i*c)*exp((c+i*c)*3),c=-infinity..infinity)/2/Pi;')% marche pas :-(
    maple eval(inttrans[invlaplace](1/S,S,t));
    maple c_x_t_inv:=value(inttrans[invlaplace](subs(q=sqrt(s),op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s),s,t))
    maple c_t_inv:=value(inttrans[invlaplace](subs(q=sqrt(s),op(0,EQ16[1])=unapply(EQ16[3],q),c_s),s,t))
    maple c_x_t_inv_x:=value(inttrans[invlaplace](subs(q=sqrt(s),op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s_x),s,t))
    maple c_x_s_t:=s*c_x_s-subs(x=0,c_x_s);
    maple c_x_s_x_t:=s*c_x_s-subs(x=0,c_x_s_x);
    maple c_x_t_inv_t:=value(inttrans[invlaplace](subs(q=sqrt(s),c_x_s_t),s,t))
    maple c_x_t_inv_x_t:=value(inttrans[invlaplace](subs(q=sqrt(s),c_x_s_x_t),s,t))
    maple tmp146:=subs(op(-1,tmp144)=op(-1,tmp143)/exp(q),tmp144)
    maple ver147:=simplify(tmp142*exp(q)/op(-1,tmp143)+tmp143*exp(q)/op(-1,tmp143)+tmp144/op(-1,tmp144))% cond de droite juste.
    maple subs(x=1,c_x_s_x)          
    a=load('totonew')
    for f=fieldnames(a)'% load pour maple...
    v=a.(f{1});if isfloat(v);if length(v)==1;s=[f{1},':=',mat2str(v)];maple(s);else;for i=1:length(v);s=[f{1},'(',num2str(i),'):=',mat2str(v(i))];;maple(s);end;end;end;
    end
    % warning('ATTENTION ceci fausse tout');maple x:=0;x=0;% seulement pour tester condition en 0 à gauche ... bon c'est fait ci-dessus plus proprement.
    if K0==.5; maple K0:=1/2; end
    maple tmp109 % affiche bien 0.022467808767594467202753105212248 - 0.1493175372570296 10^(-31)   i, a comparer a 0.0224678087675945 de chez matlab : juste !
    n5=1;n6=1;n4=1;
    maple n5:=1:n6:=1:n4:=1;
    maple n[5]:=1:n[6]:=1:n[4]:=1;
    maple [1.0110,evalf(tmp110)]
    maple [1.0111,evalf(tmp111)]
    mat2str((c1(n5) - (-1)^(n4) * (c2(n5))/(r(n5))) ) % 0-i*5.16397779494322
    maple evalf(c1(n5)-(-1)^(n4)*(c2(n5))/(r(n5))) % 0.13651109142201291765890657380973e-14    - 5.1639777949432209425896326998908 i, OK.
    maple exp_a_fois_erfc_b_ab:=diff(exp(a)*erfc(b),a)*h[1]+diff(exp(a)*erfc(b),b)*h[2]
    maple exp_a_fois_erfc_b_ab_plus_ab_ab:=diff(exp(a)*erfc(b),a,a)*h1[1]*h1[1]+2*diff(exp(a)*erfc(b),a,b)*h1[2]*h1[1]+diff(exp(a)*erfc(b),b,b)*h1[2]*h1[2]+diff(exp(a)*erfc(b),a)*h[1]+diff(exp(a)*erfc(b),b)*h[2]
    % debug derivee en t
    mat2str([r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[r(n5)^(2),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * r(n5) / sqrt(Fo)*(1/2))]])% [2.5e-06-i*9.68245836551854e-06 0.0025-i*0.00193649167310371 2.5-i*9.68245836551854 1250-i*968.245836551854]
    maple('a:=r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5);b:=((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo));h:=[r(n5)^(2),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * r(n5) / sqrt(Fo)*(1/2))]') %      a := 0.24999999999999939681755842359e-5  - 0.9682458365518550e-5   i; b := 0.0025 - 0.00193649167310371 i;h := [2.4999999999999939681755842359 - 9.682458365518550 i, 1250.0000000000000000000000000000 - 968.24583655185500000000000000000 i], OK.
    mat2str(exp_a_fois_erfc_b_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[r(n5)^(2),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)*Fo)*(-1/2) + (-1)^(n4) * r(n5) / sqrt(Fo)*(1/2))])) % -1408.00206853651-i*1082.88782304998
    maple evalf(exp_a_fois_erfc_b_ab) % -1407.9599413950666913846195667893 + 1082.8987002311311997155068782810 i, BUG ICI. bug dans exp_a_fois_erfc_b_ab ...
    % debug derivee seconde en x,x
    mat2str([r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[0,0],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))]])% [2.5e-06-i*9.68245836551854e-06 0.0025-i*0.00193649167310371 2.5-i*9.68245836551854 1250-i*968.245836551854]
    maple('a:=r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5);b:=((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo));h:=[0,0];h1:=[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))]')% a := 0.24999999999999939681755842359e-5   - 0.9682458365518550e-5i; b := 0.0025 - 0.00193649167310371i; h := [0, 0]; h1 := [-2.5 + 1.93649167310371i, -500.00000000000000000000000000000]
    mat2str(exp_a_fois_erfc_b_ab_plus_ab_ab(r(n5)^(2) * Fo + (-1)^(n4) * (1 + (-1)^(n6) * x) * r(n5),((1 + (-1)^(n6) * x)/(2 * sqrt(Fo)) + (-1)^(n4) * r(n5) * sqrt(Fo)),[0,0],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))],[ (-1)^(n4) * ( (-1)^(n6)) * r(n5),(( (-1)^(n6) )/(2 * sqrt(Fo)))]))
    maple evalf(exp_a_fois_erfc_b_ab_plus_ab_ab)%   -1407.9599413950666913846195667893 + 1082.8987002311311997155068782810i
    % verification que EQ19 satisfait la condition aux bord.
    % int(c_x_t-ctinfty,x)+C_f_t-Cftinfty=C_t0-ctinfty+Cft0-Cftinfty=0
    % int(c_x_s,x)+C_f_s=0
    % -c_x_s_x/B=K*c_x_s-C_f_s
    % -c_x_s_x/B=K*c_x_s+c_s
    % -c_x_s_x/B=K*c_x_s+c_s
    maple c_x_t_inv% marche pas.
    maple evalf(c_x_t_inv)% marche pas.
    maple c_x_s_1_sansq:=subs(q=sqrt(s),eval(subs(op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s_1)));
    maple c_x_s_2_sansq:=subs(q=sqrt(s),eval(subs(op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s_2)));
    maple c_x_s_3_sansq:=subs(q=sqrt(s),eval(subs(op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s_3)));
    maple c_x_s_4_sansq:=subs(q=sqrt(s),eval(subs(op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s_4)));
    maple c_x_s_sansq:=subs(q=sqrt(s),eval(subs(op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s)));
    % maple invl(c_x_s_1_sansq,s,t,20,18,1e5)% -0.3701 pour t=0.00151878149142892, x=1.
    % maple invl(c_x_s_2_sansq,s,t,20,18,1e5)% 0.31
    % maple invl(c_x_s_3_sansq,s,t,20,18,1e5)% <1e-21
    % maple invl(c_x_s_4_sansq,s,t,20,18,1e5)% 0.500 
    % total des 4: 0.4399 ... compatible avec 0.4447 solution analytique temps courts.
    % maple invl(c_x_s_sansq,s,1e-9,20,18,1e5)% ??
    % maple invl(c_x_s_sansq,s,2,20,18,1e5)% ??
    % maple invl(c_x_s_sansq,s,20,20,18,1e5)% ??
    % decomposition à comparer à ceci:
    [1.001 0.5 -0.37019344723613 0.5 -3.01735705410263 0.5 7.13743353094711e-74 1.0011 -0.5]
    % c'est quoi ce facteur -3.01735705410263 ????
    % FAIREla meme chose à EQ40 et EQ44 que à 


    %maple invl(subs(q=sqrt(s),eval(subs(op(0,EQ16[1])=unapply(EQ16[3],q),c_x_s))),s,t,20,18,1e5)% 0.44
    maple tmp140:=-c_x_s_x/B-(K0*c_x_s+c_s)
    maple('evalf([seq(c_x_s_x/B,q=1..20)])')
    maple('evalf([seq(K0*c_x_s+c_s,q=1..20)])')

    % C_f_tinfty=
    % C_f_s=Cftinfty+
    [-c_x_t_EQ46_x_t/B, K0*c_x_t_EQ46_t+c_x_t_EQ46_x];
    mat2str([1.054,ans,diff(ans),log10(max(abs(ans))),log10(abs(diff(ans)))])% vérification dérivée en temps de la condition de bord 
  end
  return
end
% ancienne solution ... avec Int_exp_a_fois_erfc_b ... je n'en ai plus besoin maintenant, le truc ci-dessus est le meme à .1^14 près.
x=x+100;

%{
% EQ46 := [([[c((x, u + t[0])),"=",c[t,"=",infty]+(C[f]((t[0])) - C[f][t,"=",infty])/(2) * Sum(Sum(Sum((c[1, n[5]] - (-1)^(n[4]) * (c[2, n[5]])/(sqrt(r[n[5]]))) * exp_a_fois_erfc_b(r[n[5]] * u + (-1)^(n[4]) * (1 + (-1)^(n[6]) * x) * sqrt(r[n[5]]), ((1 + (-1)^(n[6]) * x)/(2 * sqrt(u)) + (-1)^(n[4]) * sqrt(r[n[5]] * u))), n[4]= 1..2), n[6]= 1..2), n[5]= 1..2)+(1)/(2) * Sum(Int(exp(-p * y) * (exp(-(1)/(4*u) * (2 + (-1)^(n[6]) * x + (2*n[6] - 5)/(abs(2*n[6] - 5)) * y)^(2)))/(sqrt(pi * u)), y= 0..1), n[6]= 1..3) + Int((exp(-p * y))/(4) * Sum(Sum(Sum((c[3, n[5]] - (-1)^(n[4]) * (c[4, n[5]])/(sqrt(r[n[5]]))) * exp_a_fois_erfc_b(r[n[5]] * u + (-1)^(n[4]) * (2 + (-1)^(n[6]) * x + (2*n[6] - 5)/(abs(2*n[6] - 5)) * y) * sqrt(r[n[5]]), (((2 + (-1)^(n[6]) * x + (2*n[6] - 5)/(abs(2*n[6] - 5)) * y))/(2 * sqrt(u)) + (-1)^(n[4]) * sqrt(r[n[5]] * u))), n[4]= 1..2), n[6]= 1..3), n[5]= 1..2), y= 0..1)+(1)/(2) * Int(exp(-p * y) * (exp(-(1)/(4*u) * (x - y)^(2)))/(sqrt(pi * u)), y= 0..1) + (1)/(2) * Int(exp(-p * y) * (exp(-(1)/(4*u) * (x + y)^(2)))/(sqrt(pi * u)), y= 0..1) + exp(-((1 - epsilon)/(u))) * O((1))]]) ]');
%}
c_x_t_EQ46 = ctinfty+(Cft0 - Cftinfty)./(2) .* Sum(@(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1).^(n4) .* (c2(n5))./(sqrtr(n5))) .* exp_a_fois_erfc_b(r(n5) .* u + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* sqrtr(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(u)) + (-1).^(n4) .* sqrtr(n5).*sqrt( u))), 1:1), 1:2), 1:2)+(1)./(2) .* Sum(@(n6)Int_exp(@(y)-p .* y-(1)./(4.*u) .* (2 + (-1).^(n6) .* x + (2.*n6 - 5)./(abs(2.*n6 - 5)) .* y).^(2), 0:1)./(sqrt(pi .* u)), 1:3) + Sum(@(n5)Sum(@(n6)Sum(@(n4)1./(4) .* (c3(n5) - (-1).^(n4) .* (c4(n5))./(sqrtr(n5))) .* Int_exp_a_fois_erfc_b(@(y)(-p .* y + r(n5) .* u + (-1).^(n4) .* (2 + (-1).^(n6) .* x + (2.*n6 - 5)./(abs(2.*n6 - 5)) .* y) .* sqrtr(n5)) , @(y)((((2 + (-1).^(n6) .* x + (2.*n6 - 5)./(abs(2.*n6 - 5)) .* y))./(2 .* sqrt(u)) + (-1).^(n4) .* sqrtr(n5).*sqrt(u))), 0:1),  1:1), 1:3), 1:2)+(1)./(2) .* Int_exp(@(y)(-p .* y-(1)./(4.*u) .* (x - y).^(2)), 0:1)./(sqrt(pi .* u)) + (1)./(2) .* Int_exp(@(y)(-p .* y-(1)./(4.*u) .* (x + y).^(2)),  0:1)./(sqrt(pi .* u));

mat2str([1.002,ctinfty,(Cft0 - Cftinfty)./(2) .* Sum(@(n5)Sum(@(n6)Sum(@(n4)(c1(n5) - (-1).^(n4) .* (c2(n5))./(sqrtr(n5))) .* exp_a_fois_erfc_b(r(n5) .* u + (-1).^(n4) .* (1 + (-1).^(n6) .* x) .* sqrtr(n5),((1 + (-1).^(n6) .* x)./(2 .* sqrt(u)) + (-1).^(n4) .* sqrtr(n5).*sqrt( u))), 1:1), 1:2), 1:2),(1)./(2) .* Sum(@(n6)Int_exp(@(y)-p .* y-(1)./(4.*u) .* (2 + (-1).^(n6) .* x + (2.*n6 - 5)./(abs(2.*n6 - 5)) .* y).^(2), 0:1)./(sqrt(pi .* u)), 1:3) , Sum(@(n5)Sum(@(n6)Sum(@(n4)1./(4) .* (c3(n5) - (-1).^(n4) .* (c4(n5))./(sqrtr(n5))) .* Int_exp_a_fois_erfc_b(@(y)(-p .* y + r(n5) .* u + (-1).^(n4) .* (2 + (-1).^(n6) .* x + (2.*n6 - 5)./(abs(2.*n6 - 5)) .* y) .* sqrtr(n5)) , @(y)((((2 + (-1).^(n6) .* x + (2.*n6 - 5)./(abs(2.*n6 - 5)) .* y))./(2 .* sqrt(u)) + (-1).^(n4) .* sqrtr(n5).*sqrt(u))), 0:1),  1:1), 1:3), 1:2),(1)./(2) .* Int_exp(@(y)(-p .* y-(1)./(4.*u) .* (x - y).^(2)), 0:1)./(sqrt(pi .* u)) , (1)./(2) .* Int_exp(@(y)(-p .* y-(1)./(4.*u) .* (x + y).^(2)),  0:1)./(sqrt(pi .* u)) ,1.0021,(Cft0 - Cftinfty)./(2)  ,(1)./(2)   , 0,(1)./(2)  ./(sqrt(pi .* u)) , (1)./(2)./(sqrt(pi .* u)) ])
epsilon=.01;% any positive number will do.
errEQ46 = exp(-((1 - epsilon)./(u)));
save totoold

function bidon=maple_check_article_shorttime(varargin)
maple restart
run senscnonlin_v68
maple translatedlaplace:=[C[prime][f](t)=bar[C][f](s)+C[f][t,"=",infty],c(x,t)=bar[c](x,s)+c[t,"=",infty]];
maple verEQ8no16[6]:=expand(subs(op(subs(x=1,translatedlaplace)),EQ6[1]=EQ6[3],K[prime][0]=K[0],c[prime][f0]=c[f0],subs(EQ5[1]=EQ5[3],EQ4[1,1,3])-EQ8[1,1,6]));
maple tmp9:=subs(0..1=0..x,abs(x-y)=x-y,value(op(3,EQ9)))+subs(0..1=x..1,abs(x-y)=y-x,value(op([3,3],EQ9))):;
maple s:=q^2;
maple verEQ9no25[vervisuelletmp3]:=factor(expand(combine(EQ9[3]-tmp9,ranges)));
maple verEQ9no25[1]:=expand(combine(expand(diff(tmp9,x,x))-expand(subs(EQ9[1]=tmp9,EQ7[1])),ranges));
maple tmp10:=subs(0..1=0..x,abs(x-y)=x-y,value(op(3,EQ10)))+subs(0..1=x..1,abs(x-y)=y-x,value(op([3,3],EQ10))):;
maple verEQ10no26:=expand(combine(expand(diff(tmp9,x))-expand(tmp10)),ranges);
maple verEQ11no27:=expand(value(subs(x=0,tmp10)+EQ11[1]));
maple verEQ12no28:=factor(expand(value(subs(abs(1-y)=1-y,eval(subs(subs(x=1,EQ10[1]=EQ10[3]),subs(x=1,EQ9[1]=EQ9[3]),EQ8[1,1,4]-EQ8[1,1,6]-EQ12[1]+EQ12[3]))))));
maple verEQ13no29:=subs(abs(1-y)=1-y,eval(subs(subs(x=1,EQ10[1]=EQ10[3]),subs(x=1,EQ9[1]=EQ9[3]),EQ8[1,1,7]-EQ8[1,1,9]-EQ13[1]+EQ13[3])));
maple verEQ14no30:=map(expand,evalm(op([1,1],EQ14)&*op([1,2],EQ14)-EQ14[3]+[[-EQ11[1]],[EQ12[3]-EQ12[1]],[EQ13[3]-EQ13[1]]]))
maple e[3]:=unapply(EQ16[3],q):;
maple verEQ15no31:=expand([EQ15[5]-EQ15[3],EQ15[3]-det(op([1,1],EQ14))]);
maple tmp17:=solve([-EQ11[1],EQ12[3]-EQ12[1],EQ13[3]-EQ13[1]],[e[1](s),e[2](s),bar[C][f](s)]):;
maple verEQ17no32[5]:=combine(expand(EQ17[1,1,3]-EQ17[1,1,5]));
maple verEQ17no32[4]:=expand(diff(subs(0..1=0..u,expand(EQ17[1,1,5]-EQ17[1,1,7])),u));
maple verEQ17no32[3]:=expand(eval(value(subs(0..1=0..0,expand(EQ17[1,1,5]-EQ17[1,1,7])))));
maple verEQ17no32[2]:=factor(expand(diff(eval(subs(0..1=0..u,op([1,1,2],tmp17)*exp(q)*EQ15[5]*B-EQ17[1,1,7]*B)),u)));
maple verEQ17no32[1]:=factor(value(subs(0..1=0..0,op([1,1,2],tmp17)*exp(q)*EQ15[5]*B-EQ17[1,1,7]*B)));
maple verEQ18no33[5]:=combine(expand(EQ18[1,1,3]-EQ18[1,1,5]));
maple verEQ18no33[4]:=expand(diff(subs(0..1=0..u,expand(EQ18[1,1,5]-EQ18[1,1,7])),u));
maple verEQ18no33[3]:=expand(eval(value(subs(0..1=0..0,expand(EQ18[1,1,5]-EQ18[1,1,7])))));
maple verEQ18no33[2]:=factor(expand(diff(eval(subs(0..1=0..u,op([1,2,2],tmp17)*exp(-q)*EQ15[5]*B-EQ18[1,1,7]*B)),u)));
maple verEQ18no33[1]:=factor(value(subs(0..1=0..0,op([1,2,2],tmp17)*exp(-q)*EQ15[5]*B-EQ18[1,1,7]*B)));
maple tmp19:=subs(e[1](s)=EQ17[1,1,7]/EQ17[1,1,1]*e[1](s),e[2](s)=EQ18[1,1,7]/EQ18[1,1,1]*e[2](s),EQ9[3])-EQ19[1,1,3]:;
maple verEQ19no33[1]:=factor(value(subs(0..1=0..0,tmp19)));
maple verEQ19no33[2]:=simplify(factor(value(diff(subs(0..1=0..u,tmp19),u))));
maple verEQ20no35[1]:=EQ15[1]*EQ20[1]-1;
maple verEQ20no35[2]:=EQ15[5]*EQ20[3]-1;
maple verEQ20no35[3]:=simplify(expand(value(EQ20[5])-EQ20[3]));
maple verEQ20no35[4]:=simplify(expand(value(EQ20[7])-EQ20[3]));
% maple verEQ20no35[4]:=factor(EQ20[7]-EQ20[5])); marchera pas :-( pas grave, maple peut me calculer la somme exacte.
maple tmp21:=subs(0..1=0..u,simplify(value(EQ21[1,1,3])*EQ15[5]^2-subs(EQ15[1]=EQ15[5],EQ19[1,1,3])*EQ15[5]^2)):;
maple verEQ21no36[1]:=combine(factor(value(subs(u=0,tmp21))));
maple verEQ21no36[1]:=factor(expand(value(diff(tmp21,u))));
maple assume(z>0);
maple verEQ22no37:=inttrans[laplace](subs(t[0]=t-u,EQ22[1,1,3]-c[t,"=",infty]),u,s)-EQ22[1,1,7];
maple verEQ23no38:=inttrans[laplace](subs(t[0]=t-u,EQ23[1,1,3]-c[t,"=",infty]),u,s)-EQ23[1,1,7];
maple verEQ24no39:=[expand(EQ24[1]-EQ24[3]),expand(EQ24[5]-EQ24[3])];
maple verEQ26no41:=eval(inttrans[laplace](subs(t[0]=0,c[t,"=",infty]=0,t[prime]=tprime,eval(EQ26[1,1,3])),t,s)-subs(bar[c][6]=unapply(laplace(c[6](x, t), t, S),x,S),EQ26[1,1,7]));
'EQ27 non vérifiée ... mais vrai selon th VI de pdf171 p300 de Carslaw59conduction'
%maple eval(inttrans[laplace](subs(t[0]=0,c[t,"=",infty]=0,t[prime]=tprime,eval(EQ27[1,1,3])),t,s)-subs(bar[c][8]=unapply(laplace(c[8](x, t), t, S),x,S),EQ27[1,1,7]));
%maple eval(inttrans[invlaplace](subs(bar[c][8]=unapply(laplace(c[8](x, t), t, S),x,S),EQ27[1,1,7]),S,t));
% maple verEQ28no43:=eval(inttrans[laplace](eval(subs(u=t,c[t,"=",infty]=0,eval(EQ28[1,1,3]))),t,s)-EQ28[1,1,7]) marche pas.
maple verEQ28no43[1]:=expand(subs(c[8](x,t)=EQ22[1,1,3],EQ27[1,1,3])-subs(u=t-t[0],EQ28[1,1,3])-EQ27[1,1,1]+c[7](x,t));
maple verEQ28no43[2]:=eval(subs(bar[c][8]=unapply(subs(s=S,EQ22[1,1,7]),x,S),EQ27[1,1,7])-EQ28[1,1,7]);
maple verEQ29no44[1]:=expand(subs(c[8](x,t)=EQ23[1,1,3],EQ27[1,1,3])-subs(u=t-t[0],EQ29[1,1,3])-EQ27[1,1,1]+c[7](x,t));
maple verEQ29no44[2]:=eval(subs(bar[c][8]=unapply(subs(s=S,EQ23[1,1,7]),x,S),EQ27[1,1,7])-EQ29[1,1,7]);

maple assume(zu>0);
maple tmp30:=combine(subs(c[6](x,t[prime])=EQ28[1,1,3],u=tprime-t[0],t[prime]=tprime,c[tprime+u,"=",infty]=c[t,"=",infty],EQ26[1,1,3])-subs(u=t-t[0],student[changevar](u[prime]=tprime-t[0],EQ30[1,1,3],tprime))-EQ26[1,1,1]+c[5](x,t)):;
maple verEQ30no45[1]:=factor(subs(t=t0+zu,eval(factor(value(diff(subs(t[0]=t0,tmp30),t))))));
maple verEQ30no45[2]:=value(subs(t[0]=t0,t=t0,tmp30));
maple verEQ30no45[3]:=eval(subs(bar[c][6]=unapply(subs(s=S,EQ28[1,1,7]),x,S),EQ26[1,1,7])-EQ30[1,1,7]);

maple tmp31:=combine(subs(c[6](x,t[prime])=EQ29[1,1,3],u=tprime-t[0],t[prime]=tprime,c[tprime+u,"=",infty]=c[t,"=",infty],EQ26[1,1,3])-subs(u=t-t[0],student[changevar](u[prime]=tprime-t[0],EQ31[1,1,3],tprime))-EQ26[1,1,1]+c[5](x,t)):;
maple verEQ31no46[1]:=factor(value(diff(subs(t[0]=t0,tmp31),t)));
maple verEQ31no46[2]:=value(subs(t[0]=t0,t=t0,tmp31));
maple verEQ31no46[3]:=eval(subs(bar[c][6]=unapply(subs(s=S,EQ29[1,1,7]),x,S),EQ26[1,1,7])-EQ31[1,1,7]);

maple assume(ni::integer);
maple tmp32:=combine(subs(c[6](x,t[prime])=student[changevar](u[prime]=vprime-t[0],EQ32[1,1,3],vprime),u[prime]=uprime,u=tprime-t[0],t[prime]=tprime,c[tprime+u,"=",infty]=c[t,"=",infty],EQ26[1,1,3])-subs(u[prime]=uprime,u=t-t[0],student[changevar](u[prime]=tprime-t[0],subs(n=n+1,EQ32[1,1,3]),tprime))-EQ26[1,1,1]+c[5](x,t)):;
maple verEQ32no45[1]:=subs(n*(n-1)!=n!,expand(factor(combine(factor(subs(t=t0+zu,vprime=tprime,eval(factor(value(diff(subs(t[0]=t0,tmp32),t))))))))));
maple verEQ32no45[2]:=value(subs(t[0]=t0,t=t0,tmp32));
maple verEQ32no45[3]:=expand(factor(eval(subs(bar[c][6]=unapply(subs(s=S,EQ32[1,1,7]),x,S),n=ni,EQ26[1,1,7])-subs(n=ni+1,EQ32[1,1,7]))));
maple verEQ32no45[4]:=[value(eval(factor(combine(subs(n=1,EQ32[1,1,3])-EQ30[1,1,3])))),subs(n=1,EQ32[1,1,7])-EQ30[1,1,7]];

maple tmp33:=combine(subs(c[6](x,t[prime])=student[changevar](u[prime]=vprime-t[0],EQ33[1,1,3],vprime),u[prime]=uprime,u=tprime-t[0],t[prime]=tprime,c[tprime+u,"=",infty]=c[t,"=",infty],EQ26[1,1,3])-subs(u[prime]=uprime,u=t-t[0],student[changevar](u[prime]=tprime-t[0],subs(n=n+1,EQ33[1,1,3]),tprime))-EQ26[1,1,1]+c[5](x,t)):;
maple verEQ33no46[1]:=subs(n*(n-1)!=n!,expand(factor(combine(factor(subs(t=t0+zu,vprime=tprime,eval(factor(value(diff(subs(t[0]=t0,tmp33),t))))))))));
maple verEQ33no46[2]:=value(subs(t[0]=t0,t=t0,tmp33));
maple verEQ33no46[3]:=expand(factor(eval(subs(bar[c][6]=unapply(subs(s=S,EQ33[1,1,7]),x,S),n=ni,EQ26[1,1,7])-subs(n=ni+1,EQ33[1,1,7]))));
maple verEQ33no46[4]:=[value(eval(factor(combine(subs(n=1,EQ33[1,1,3])-EQ31[1,1,3])))),subs(n=1,EQ33[1,1,7])-EQ31[1,1,7]];

maple verEQ34no49[1]:=expand(subs(c[7](x,t)=EQ32[1,1,3],u[prime]=uprime,t[0]=t0,t=t0+u,n=1+ni,eval(solve(EQ27[1,1,1]=EQ27[1,1,3],c[8](x,t))))-subs(u[prime]=uprime,t=t0+u,n=1+ni,EQ34[1,1,3]));
maple verEQ34no49[2]:=expand(eval(subs(bar[c][7]=unapply(subs(s=S,EQ32[1,1,7]),x,S),eval(solve(subs(q^2=q^2-r,EQ27[1,1,5]=EQ27[1,1,7]),bar[c][8](x,s)))))-EQ34[1,1,7]);

maple verEQ35no50[1]:=value(combine(expand(subs(c[7](x,t)=EQ33[1,1,3],u[prime]=uprime,t[0]=t0,t=t0+u,n=1+ni,eval(solve(EQ27[1,1,1]=EQ27[1,1,3],c[8](x,t))))-subs(u[prime]=uprime,t=t0+u,n=1+ni,EQ35[1,1,3]))));
maple verEQ35no50[2]:=expand(eval(subs(bar[c][7]=unapply(subs(s=S,EQ33[1,1,7]),x,S),eval(solve(subs(q^2=q^2-r,EQ27[1,1,5]=EQ27[1,1,7]),bar[c][8](x,s)))))-EQ35[1,1,7]);

'verEQ36no51 non verifie, mais il n''est utilisé que pour EQ37no53 qui lui est vérifié indépendamment'
% maple tmp36:=inttrans[invlaplace](eval(subs(q^2=S,EQ36[1,1,7])),S,t)-subs(t[0]=0,u=t,c[t,"=",infty]=0,EQ36[1,1,3])
maple verEQ37no52[1]:=subs(sqrt(z*r)=sqrt(z)*sqrt(r),sqrt(z/u)=sqrt(z)/sqrt(u),sqrt(r*u)=sqrt(r)*sqrt(u),(z/u)^(-1/2)=sqrt(u)/sqrt(z),sqrt(Pi/z)=sqrt(Pi)/sqrt(z),eval(expand(value(diff(subs(u[prime]=uprime,EQ37[1]-EQ37[3]),u)))));
maple verEQ37no52[2]:=value(evalf(subs(u[prime]=uprime,u=1e-100,r=2+sqrt(-2),z=4,EQ37[1]-EQ37[3])));
maple verEQ37no52[3]:=value(evalf(subs(u[prime]=uprime,u=1e-100,r=-4-sqrt(-5),z=27,EQ37[1]-EQ37[3])));
maple verEQ38no53[1]:=subs(sqrt(z*r)=sqrt(z)*sqrt(r),sqrt(z/u)=sqrt(z)/sqrt(u),sqrt(r*u)=sqrt(r)*sqrt(u),(r*u)^(-1/2)=1/sqrt(r)/sqrt(u),sqrt(Pi/r)=sqrt(Pi)/sqrt(r),eval(expand(value(diff(subs(u[prime]=uprime,EQ38[1]-EQ38[3]),u)))));
maple verEQ38no53[2]:=value(evalf(subs(u[prime]=uprime,u=1e-100,r=2+sqrt(-2),z=4,EQ38[1]-EQ38[3])));
maple verEQ38no53[3]:=value(evalf(subs(u[prime]=uprime,u=1e-100,r=-4-sqrt(-5),z=27,EQ38[1]-EQ38[3])));

maple verEQ39no55[1]:=map(expand,factor(value(diff(subs(u[prime]=uprime,EQ39[1]-EQ39[3]),u))))
maple verEQ39no54[2]:=value(evalf(subs(u[prime]=uprime,u=1e-100,r=2+sqrt(-2),z=4,n=2.35,EQ39[1]-EQ39[3])));
maple verEQ39no54[3]:=value(evalf(subs(u[prime]=uprime,u=1e-100,r=-4-sqrt(-5),z=27,n=39.58,EQ39[1]-EQ39[3])));
%maple verEQ40no55:=series(value(subs(q=-log(lq),epsilon=0,x=x01,EQ40[1,1,3]-subs(infinity=2,c(y,t[0])=c[t,"to",infty]+exp(-p*y),EQ21[1,1,3]))),lq=0,2)
maple assume(x01::RealRange(0,1));
%maple verEQ40no55:=combine(simplify(subs(epsilon=0,x=x01,C[f](t[0])=C[f][t, "=", infty]+C,e[3](q)=e3,e[3](-q)=em,EQ40[1,1,3]-subs(infinity=0,c(y,t[0])=c[t,"to",infty]+exp(-p*y),EQ21[1,1,3]))));
maple verEQ40no55:=expand(combine(subs(keepInt=Int,value(subs(epsilon=0,x=x01,C[f](t[0])=C[f][t, "=", infty]+C,e[3](q)=e3,e[3](-q)=em,Int=keepInt,EQ40[1,1,3]-subs(infinity=0,c(y,t[0])=c[t,"=",infty]+exp(-p*y),EQ21[1,1,3]))))));
' on ne verifie pas que les termes correspondant a n1>0 sont en O comme il faut dans EQ21'
% maple verEQ40no55:=series(subs(q=-log(lq),expand(combine(subs(keepInt=Int,value(subs(epsilon=0,x=x01,C[f](t[0])=C[f][t, "=", infty]+C,e[3](q)=e3,e[3](-q)=em,Int=keepInt,EQ40[1,1,3]-subs(infinity=2,c(y,t[0])=c[t,"=",infty]+exp(-p*y),EQ21[1,1,3]))))))),lq=0,2);
maple verEQ41no56[1]:=factor(EQ41[1,1,3]-EQ41[1,1,1]);
maple r1:=subs(Roots(plusminus^(2) - 1)=1,EQ42[1]):;
maple r2:=subs(Roots(plusminus^(2) - 1)=-1,EQ42[1]):;
maple verEQ42:=expand([e[3](r1),e[3](r2)]);
maple r1:=r1^2:;
maple r2:=r2^2:;
maple verEQ41no56[2]:=factor(op([1,1,3],EQ41)-subs(r[1]=r1,r[2]=r2,op([1,1,5],EQ41)));
maple verEQ41no56[3]:=factor(value(op([1,1,7],EQ41)-op([1,1,5],EQ41)));
maple c1n5:=solve(op([1,1,9,1,1],EQ41)=op([1,1,7,1,1],EQ41),c[1,n[5]]):;
maple c2n5:=solve(op([1,1,9,2,2,1],EQ41)=op([1,1,7,2,2,1],EQ41),c[2,n[5]]):;
maple verEQ41no56[4]:=factor(value(op([1,1,7],EQ41)-subs(c[1,n[5]]=c1n5,c[2,n[5]]=c2n5,op([1,1,9],EQ41))));
%maple verEQ41no56[2.1]:=factor(op([1,1,3,1],EQ41)-subs(r[1]=r1,r[2]=r2,op([1,1,5,1],EQ41)+op([1,1,5,2],EQ41)));
%maple verEQ41no56[2.2]:=factor(op([1,1,3,2],EQ41)-subs(r[1]=r1,r[2]=r2,op([1,1,5,3],EQ41)));
maple verEQ43no57[1]:=factor(EQ43[1,1,3]-EQ43[1,1,1]);
maple verEQ43no57[2]:=factor(op([1,1,3],EQ43)-subs(r[1]=r1,r[2]=r2,op([1,1,5],EQ43)));
maple verEQ43no57[3]:=factor(value(op([1,1,7],EQ43)-op([1,1,5],EQ43)));
maple c3n5:=solve(op([1,1,9,2,1],EQ43)=op([1,1,7,2,1],EQ43),c[3,n[5]]):;
maple c4n5:=solve(op([1,1,9,3,2,1],EQ43)=op([1,1,7,3,2,1],EQ43),c[4,n[5]]):;
maple verEQ43no56[4]:=factor(value(op([1,1,7],EQ43)-subs(c[3,n[5]]=c3n5,c[4,n[5]]=c4n5,op([1,1,9],EQ43))));
maple verEQ44no58:=factor(subs(r[1]=r1,r[2]=r2,value(op([1,1,3],EQ40)-subs(c[1,n[5]]=c1n5,c[2,n[5]]=c2n5,c[3,n[5]]=c3n5,c[4,n[5]]=c4n5,op([1,1,3],EQ44)))));
maple assume(xp>0);
maple assume(yp>0);
maple assume(qp>0);
maple verEQ45no59[manuelle1]:=(expand(subs(abs(x-y)^2=(x-y)^2,expand(combine(subs(iint=Int,t[0]=t0,t=u+t0,value(subs(c[1,n[5]]=0,c[2,n[5]]=0,c[3,n[5]]=0,c[4,n[5]]=0,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,EQ23[1,1,7]*sqrt(s)=(EQ23[1,1,3]-c[t,"=",infty])*q),expo=[op([1,1,3,2,2,1,2,1,1,1,1],EQ44),op([1,1,3,3,2,1,2,1,1,1],EQ44),op([1,1,3,3,2,1,2,2,1,1],EQ44)]),Int=iint,EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3]))))))));
maple verEQ45no59[2]:=subs(Int(0,u[prime]=0..u)=0,eval(factor(subs(0!=1,sqrt((1+x)^2/pi)=(1+x)/sqrt(pi),sqrt((x-1)^2/pi)=(1-x)/sqrt(pi),sqrt((1-x)^2/pi)=(1-x)/sqrt(pi),1/sqrt(pi*u[prime])=1/sqrt(pi)/sqrt(u[prime]),u[prime]^(-3/2)=1/u[prime]/sqrt(u[prime]),eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ34[1,1,7]*(s-r)=(EQ34[1,1,3]-c[t,"=",infty])*(s-r)),expo=[op([1,1,3,1,2,1,1,1,2,1],EQ44)]),Int=iint,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[1,n[5]]))))))))));
maple verEQ45no59[3]:=subs(Int(0,u[prime]=0..u)=0,eval(factor(subs(0!=1,sqrt((1+x)^2/pi)=(1+x)/sqrt(pi),sqrt((1-x)^2/pi)=(1-x)/sqrt(pi),sqrt((x-1)^2/pi)=(1-x)/sqrt(pi),1/sqrt(pi*u[prime])=1/sqrt(pi)/sqrt(u[prime]),u[prime]^(-3/2)=1/u[prime]/sqrt(u[prime]),eval(eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ34[1,1,7]*(s-r)=(EQ34[1,1,3]-c[t,"=",infty])*(s-r)),expo=[op([1,1,3,1,2,1,1,1,2,1],EQ44)]),Int=iint,1..2=1..1,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[1,n[5]])))))))))));% la dérivée ci-dessus additionne en fait les dérivée en c[1,1] et c[1,2], donc on vérifie indpéndamment deux morceaux de dérivées ici et ci-dessous.
maple verEQ45no59[4]:=subs(Int(0,u[prime]=0..u)=0,eval(factor(subs(0!=1,sqrt((1+x)^2/pi)=(1+x)/sqrt(pi),sqrt((x-1)^2/pi)=(1-x)/sqrt(pi),1/sqrt(pi*u[prime])=1/sqrt(pi)/sqrt(u[prime]),u[prime]^(-3/2)=1/u[prime]/sqrt(u[prime]),eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ34[1,1,7]*(s-r)=(EQ34[1,1,3]-c[t,"=",infty])*(s-r)),expo=[op([1,1,3,1,2,1,1,1,2,1],EQ44)]),Int=iint,1..2=2..2,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[1,n[5]]))))))))));
maple verEQ45no59[5]:=subs(Int(Int(0,u[prime]=0..u),y=0..1)=0,eval(factor(subs(0!=1,sqrt((2-x-y)^2/pi)=(2-x-y)/sqrt(pi),sqrt((-2+x+y)^2/pi)=(2-x-y)/sqrt(pi),sqrt((2-x+y)^2/pi)=(2-x+y)/sqrt(pi),sqrt((-2+x-y)^2/pi)=(2-x+y)/sqrt(pi),sqrt((2+x-y)^2/pi)=(2+x-y)/sqrt(pi),1/sqrt(pi*u[prime])=1/sqrt(pi)/sqrt(u[prime]),u[prime]^(-3/2)=1/u[prime]/sqrt(u[prime]),eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ34[1,1,7]*(s-r)=(EQ34[1,1,3]-c[t,"=",infty])*(s-r)),expo=[op([1,1,3,2,2,1,2,2,1,1,2,1],EQ44)]),Int=iint,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[3,n[5]]))))))))));
maple verEQ45no59[6]:=subs(Int(Int(0,u[prime]=0..u),y=0..1)=0,eval(factor(subs(0!=1,sqrt((2-x-y)^2/pi)=(2-x-y)/sqrt(pi),sqrt((-2+x-y)^2/pi)=(2-x+y)/sqrt(pi),sqrt((-2+x+y)^2/pi)=(2-x-y)/sqrt(pi),sqrt((2-x+y)^2/pi)=(2-x+y)/sqrt(pi),sqrt((2+x-y)^2/pi)=(2+x-y)/sqrt(pi),1/sqrt(pi*u[prime])=1/sqrt(pi)/sqrt(u[prime]),u[prime]^(-3/2)=1/u[prime]/sqrt(u[prime]),eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ34[1,1,7]*(s-r)=(EQ34[1,1,3]-c[t,"=",infty])*(s-r)),expo=[op([1,1,3,2,2,1,2,2,1,1,2,1],EQ44)]),Int=iint,1..2=1..1,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[3,n[5]]))))))))));
maple verEQ45no59[7]:=subs(Int(Int(0,u[prime]=0..u),y=0..1)=0,eval(factor(subs(0!=1,sqrt((2-x-y)^2/pi)=(2-x-y)/sqrt(pi),sqrt((-2+x-y)^2/pi)=(2-x+y)/sqrt(pi),sqrt((-2+x+y)^2/pi)=(2-x-y)/sqrt(pi),sqrt((2-x+y)^2/pi)=(2-x+y)/sqrt(pi),sqrt((2+x-y)^2/pi)=(2+x-y)/sqrt(pi),1/sqrt(pi*u[prime])=1/sqrt(pi)/sqrt(u[prime]),u[prime]^(-3/2)=1/u[prime]/sqrt(u[prime]),eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ34[1,1,7]*(s-r)=(EQ34[1,1,3]-c[t,"=",infty])*(s-r)),expo=[op([1,1,3,2,2,1,2,2,1,1,2,1],EQ44)]),Int=iint,1..2=2..2,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[3,n[5]]))))))))));
maple verEQ45no59[8]:=subs(Int(0,u[prime]=0..u)=0,eval(factor(eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ35[1,1,7]*(s-r)*sqrt(s)=(EQ35[1,1,3]-c[t,"=",infty])*(s-r)*q),expo=[op([1,1,3,1,2,1,1,1,2,1],EQ44)]),Int=iint,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[2,n[5]])))))))));
maple verEQ45no59[9]:=subs(Int(0,u[prime]=0..u)=0,eval(factor(eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ35[1,1,7]*(s-r)*sqrt(s)=(EQ35[1,1,3]-c[t,"=",infty])*(s-r)*q),expo=[op([1,1,3,1,2,1,1,1,2,1],EQ44)]),Int=iint,1..2=2..2,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[2,n[5]])))))))));
maple verEQ45no59[10]:=subs(Int(0,u[prime]=0..u)=0,eval(factor(eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ35[1,1,7]*(s-r)*sqrt(s)=(EQ35[1,1,3]-c[t,"=",infty])*(s-r)*q),expo=[op([1,1,3,1,2,1,1,1,2,1],EQ44)]),Int=iint,1..2=1..1,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[2,n[5]])))))))));
maple verEQ45no59[11]:=factor(eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ35[1,1,7]*(s-r)*sqrt(s)=(EQ35[1,1,3]-c[t,"=",infty])*(s-r)*q),expo=[op([1,1,3,2,2,1,2,2,1,1,2,1],EQ44)]),Int=iint,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[4,n[5]])))))));
maple verEQ45no59[12]:=factor(eval(combine(subs(iint=Int,value(subs(n[6]=n6,seq(subs(sqrt(s*z)=-expo/2,z=(-expo/q/2)^2,n[6]=n6,n=1,r=r[n[5]],EQ35[1,1,7]*(s-r)*sqrt(s)=(EQ35[1,1,3]-c[t,"=",infty])*(s-r)*q),expo=[op([1,1,3,2,2,1,2,2,1,1,2,1],EQ44)]),Int=iint,1..2=2..2,diff(EQ45[1,1,3]-c[t,"=",infty]-EQ44[1,1,3],c[4,n[5]])))))));
% TODO dans verEQ45? on peut faire la différentiation en der plutot que les triplements de litnes ci-dessus.

maple verEQ46no60[1]:=(expand(subs(abs(x-y)^2=(x-y)^2,expand(combine(subs(iint=Int,t[0]=t0,t=u+t0,value(subs(c[1,n[5]]=0,c[2,n[5]]=0,c[3,n[5]]=0,c[4,n[5]]=0,Int=iint,EQ45[1,1,3]-EQ46[1,1,3]))))))));
maple verEQ45no59[2]:=combine(frontend(expand,[subs(iint=Int,diff(value(subs(seq(subs(sqrt(z/u)=expo,sqrt(Pi/z)=sqrt(Pi)/expo/sqrt(u),(Pi/z)^(-1/2)=expo*sqrt(u)/sqrt(Pi),sqrt(z*r)=expo*sqrt(u)*sqrt(r),z=u*expo^2,sqrt(r*u)=-(-1)^n[4]*sqrt(r*u),r=r[n[5]],Int=iint,op([3,1,2],value(EQ37[3]))=piecewise(n[4]=2,expand(op([3,1,2],value(EQ37[3]))),expand(solve(EQ37[1]=value(EQ37[3]),op([3,1,2],value(EQ37[3])))))),expo=[op([1,1,3,2,3,1,1,1,3,1,1],EQ46)]),c[1,n[5]]=c[1,n[5]]*der,c[2,n[5]]=0,c[3,n[5]]=0,c[4,n[5]]=0,Int=iint,EQ45[1,1,3]-c[t,"=",infty]-EQ46[1,1,3])),der))]));
% c[1,2]=0,c[1,1]=1,C[f][t, "=", infty]=0,C[f](t[0])=1,
maple verEQ46no60[3]:=combine(frontend(expand,[subs(iint=Int,diff(value(subs(seq(subs(sqrt(z/u)=expo,sqrt(Pi/z)=sqrt(Pi)/expo/sqrt(u),(Pi/z)^(-1/2)=expo*sqrt(u)/sqrt(Pi),sqrt(z*r)=expo*sqrt(u)*sqrt(r),(Pi/r)^(-1/2)=sqrt(r)/sqrt(Pi),(Pi*u[prime])^(-1/2)=1/sqrt(Pi)/sqrt(u[prime]),z=u*expo^2,sqrt(r*u)=-(-1)^n[4]*sqrt(r*u),r=r[n[5]],Int=iint,op([3,1,2],value(EQ37[3]))=piecewise(n[4]=2,expand(op([3,1,2],value(EQ37[3]))),expand(solve(EQ38[1]=value(EQ38[3]),op([3,1,2],value(EQ37[3])))))),expo=[op([1,1,3,2,3,1,1,1,3,1,1],EQ46)]),c[1,n[5]]=0,c[2,n[5]]=c[2,n[5]]*der,c[3,n[5]]=0,c[4,n[5]]=0,Int=iint,EQ45[1,1,3]-c[t,"=",infty]-EQ46[1,1,3])),der))]));
% maple('`diff/iint`:=parse("pr" "oc(f,xlim,x)iint(diff(f,x),xlim)+piecewise(op(0,xlim)=op(0,op=nops),apply(unapply(f,x),op([2,2],xlim))*diff(op([2,2],xlim),x)-apply(unapply(f,x),op([2,1],xlim))*diff(op([2,1],xlim),x),0) end pr" "oc:");') % non parait fauser les calculs :-(
maple tmp46:=subs(iint=Int,value(diff(value(subs(seq(subs(sqrt(z/u)=expo,sqrt(Pi/z)=sqrt(Pi)/expo/sqrt(u),(Pi/z)^(-1/2)=expo*sqrt(u)/sqrt(Pi),sqrt(z*r)=expo*sqrt(u)*sqrt(r),(Pi/r)^(-1/2)=sqrt(r)/sqrt(Pi),(Pi*u[prime])^(-1/2)=1/sqrt(Pi)/sqrt(u[prime]),z=u*expo^2,sqrt(r*u)=-(-1)^n[4]*sqrt(r*u),r=r[n[5]],Int=iint,op([3,1,2],value(EQ37[3]))=piecewise(n[4]=2,expand(op([3,1,2],value(EQ37[3]))),expand(solve(EQ37[1]=value(EQ37[3]),op([3,1,2],value(EQ37[3])))))),expo=[op([1,1,3,4,1,3,1,1,1,3,1,1],EQ46)]),c[1,n[5]]=0,c[2,n[5]]=0,c[3,n[5]]=c[3,n[5]]*der,c[4,n[5]]=0,Int=iint,u[prime]=v,EQ45[1,1,3]-c[t,"=",infty]-EQ46[1,1,3])),der))):;
maple tmp461:=subs(0..1=0..m,tmp46):;
maple verEQ46no60[4]:=simplify(diff(tmp461,m));
maple verEQ46no60[5]:=value(subs(0..1=0..0,tmp46));
maple tmp463:=subs(iint=Int,value(diff(value(subs(seq(subs(sqrt(z/u)=expo,sqrt(Pi/z)=sqrt(Pi)/expo/sqrt(u),(Pi/z)^(-1/2)=expo*sqrt(u)/sqrt(Pi),sqrt(z*r)=expo*sqrt(u)*sqrt(r),(Pi/r)^(-1/2)=sqrt(r)/sqrt(Pi),(Pi*u[prime])^(-1/2)=1/sqrt(Pi)/sqrt(u[prime]),z=u*expo^2,sqrt(r*u)=-(-1)^n[4]*sqrt(r*u),r=r[n[5]],Int=iint,op([3,1,2],value(EQ37[3]))=piecewise(n[4]=2,expand(op([3,1,2],value(EQ37[3]))),expand(solve(EQ38[1]=value(EQ38[3]),op([3,1,2],value(EQ37[3])))))),expo=[op([1,1,3,4,1,3,1,1,1,3,1,1],EQ46)]),c[1,n[5]]=0,c[2,n[5]]=0,c[3,n[5]]=0,c[4,n[5]]=c[4,n[5]]*der,Int=iint,u[prime]=v,EQ45[1,1,3]-c[t,"=",infty]-EQ46[1,1,3])),der))):;
maple tmp464:=subs(0..1=0..m,tmp463):;
maple verEQ46no60[6]:=simplify(diff(tmp464,m));
maple verEQ46no60[7]:=value(subs(0..1=0..0,tmp463));
maple verEQ47no61[1]:=expand(diff(EQ47[3]-EQ47[1],y[prime]))
maple verEQ47no61[2]:=value(expand(subs(y[prime]=0,EQ47[3]-EQ47[1])))
maple verEQ47no61[3]:=value(EQ47[3]-EQ47[5])
maple verEQ48no62[1]:=expand(diff(EQ48[1,1,3,1,1]-EQ48[1,1,1],y[prime]))
maple verEQ48no62[2]:=factor(expand(subs(p=0,diff(EQ48[1,1,3,2,1]-EQ48[1,1,1],y[prime]))))
maple verEQ48no62[3]:=value(expand(subs(y[prime]=0,EQ48[1,1,3,1,1]-EQ48[1,1,1])))
maple verEQ48no62[4]:=expand(value(EQ48[1,1,3,1,1]-EQ48[1,1,5]))

%%% horreur !!!
% maple 1/2*diff(iint(exp(-p*y)*(1/2*exp(r[1]*u)*(2-x-y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*v-1/4/v*(2-x-y)^2)/v^(3/2),v=0..u)+1/2*exp(r[1]*u)*(2+x-y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*v-1/4/v*(2+x-y)^2)/v^(3/2),v=0..u)+1/2*exp(r[1]*u)*(2-x+y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*v-1/4/v*(2-x+y)^2)/v^(3/2),v=0..u)+1/2*exp(r[2]*u)*(2-x-y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*v-1/4/v*(2-x-y)^2)/v^(3/2),v=0..u)+1/2*exp(r[2]*u)*(2+x-y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*v-1/4/v*(2+x-y)^2)/v^(3/2),v=0..u)+1/2*exp(r[2]*u)*(2-x+y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*v-1/4/v*(2-x+y)^2)/v^(3/2),v=0..u)),y=0..1),der) % est bien différentié.
% maple 1/2*diff(iint(exp(-p*y)*(1/2*exp(r[1]*u)*(2+x-y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*u[prime]-1/4/u[prime]*(2+x-y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[1]*u)*(2-x+y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*u[prime]-1/4/u[prime]*(2-x+y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[2]*u)*(2-x-y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*u[prime]-1/4/u[prime]*(2-x-y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[2]*u)*(2+x-y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*u[prime]-1/4/u[prime]*(2+x-y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[2]*u)*(2-x+y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*u[prime]-1/4/u[prime]*(2-x+y)^2)/u[prime]^(3/2),u[prime]=0..u)),y=0..1),der) % est bien différentié
% maple 1/2*diff(iint(exp(-p*y)*(1/2*exp(r[1]*u)*(2-x-y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*u[prime]-1/4/u[prime]*(2-x-y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[1]*u)*(2+x-y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*u[prime]-1/4/u[prime]*(2+x-y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[1]*u)*(2-x+y)/pi^(1/2)*c[3,1]*der*iint(exp(-r[1]*u[prime]-1/4/u[prime]*(2-x+y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[2]*u)*(2-x-y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*u[prime]-1/4/u[prime]*(2-x-y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[2]*u)*(2+x-y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*u[prime]-1/4/u[prime]*(2+x-y)^2)/u[prime]^(3/2),u[prime]=0..u)+1/2*exp(r[2]*u)*(2-x+y)/pi^(1/2)*c[3,2]*der*iint(exp(-r[2]*u[prime]-1/4/u[prime]*(2-x+y)^2)/u[prime]^(3/2),u[prime]=0..u)),y=0..1),der) % n'est pas différentié :-(
% donc quelque part, diff refuse de différentier les choses trop longues :-(
function notes

%{
% commentaires par blocs pour matlab.
'; EQ41 := (((((1)/(q * e(3)((q))),"=",(1)/(q) * (s/B + 1)/(s^(2)/B^(2) + s * (2/B - K(0)^(2)) + 1) - (K(0))/(s^(2)/B^(2) + s * (2/B - K(0)^(2)) + 1),"=",(1)/(q) * (1)/(s - r(2)) * (B * r(2) + B^(2))/(r(2) - r(1)) - (1)/(q) * (1)/(s - r(1)) * (B * r(1) + B^(2))/(r(2) - r(1)) + (B^(2) * K(0))/(r(1) - r(2)) * ((1)/(s - r(2)) - (1)/(s - r(1))),"
",Sum((1)/(s - r(n(5))) * (B^(2) * K(0))/(r(3 - n(5)) - r(n(5))), n(5)= 1..2) + (1)/(q) * Sum((1)/(s - r(n(5))) * (B * r(n(5)) + B^(2))/(r(n(5)) - r(3 - n(5))), n(5)= 1..2),"
",Sum((c(1, n(5)))/(s - r(n(5))), n(5)= 1..2) + (1)/(q) * Sum((c(2, n(5)))/(s - r(n(5))), n(5)= 1..2)))) )');
'; EQ43 := (((((e(3)((-q)))/(q * e(3)((q))),"=",(1)/(q) * (s^(2)/B^(2) + s * (2/B + K(0)^(2)) + 1)/(s^(2)/B^(2) + s * (2/B -K(0)^(2)) + 1) - 2 * (s * K(0)/B + K(0))/(s^(2)/B^(2) + s * (2/B -K(0)^(2)) + 1),"=",(1)/(q) + (1)/(q) * ((r(1))/(s - r(1)) - (r(2))/(s - r(2))) * (2*K(0)^(2) * B^(2))/(r(1) - r(2)) - (2)/(s - r(1)) * (r(1) * K(0) * B + B^(2) * K(0))/(r(1) - r(2)) - (2)/(s - r(2)) * (r(2) * K(0) * B + B^(2) * K(0))/(r(2) - r(1)),"=",(1)/(q) + Sum((2)/(s - r(n(5))) * (r(n(5)) * K(0) * B + B^(2) * K(0))/(r(3 - n(5)) - r(n(5))), n(5)= 1..2) + (1)/(q) * Sum((r(n(5)))/(s - r(n(5))) * (2*K(0)^(2) * B^(2))/(r(n(5)) - r(3 - n(5))), n(5)= 1..2),"=",(1)/(q) + Sum((c(3, n(5)))/(s - r(n(5))), n(5)= 1..2) + (1)/(q) * Sum((c(4, n(5)))/(s - r(n(5))), n(5)= 1..2)))) )');
B^2*K0/(r[3-n[5]]-r[n[5]])
>> maple c2n5
ans =
-B*(r[n[5]]+B)/(r[3-n[5]]-r[n[5]])
>> maple c3n5
ans =
2*K0*B*(r[n[5]]+B)/(r[3-n[5]]-r[n[5]])
>> maple c4n5
ans =
-2*r[n[5]]*B^2*K0^2/(r[3-n[5]]-r[n[5]])


maple r1
(-1/2*K0*B+1/2*(B^2*K0^2-4*B)^(1/2))^2
>> maple r2
ans =
(-1/2*K0*B-1/2*(B^2*K0^2-4*B)^(1/2))^2


'<,'>s:\[:(:g
'<,'>s:\]:):g
'<,'>s:K\[0]:K0:g
'<,'>s:r[n[5]]:r:g
'<,'>s:r\[n\[5]]:r:ge
'<,'>s:r\[3-n\[5]]:r:ge
'<,'>s:r\[3-n\[5]]:r([2 1]):ge
'<,'>s:c(\([1-4]\), n(5)):c\1n5:g
'<,'>s:
'<,'>s:[^*/]:.&:g
s:[^*/]:.&:g
s:[*^/]:.&:g
'<,'>s:sqrt(r:sqrtr(:g
'<,'>s:sqrt(r(n(5))):sqrtr:g
w
%}
% commentaires par blocs pour matlab.

function res=Int_exp(funca,yspan)
  % calcule Int(exp(funca(y)),y=yspan(1)..yspan(2));
  % suppose que funca est de la forme -(a*y+b)^2/4+c
% EQ47
% (Int(exp(-1./4.*(a.*y+b).^2)./exp(p.*y),y = 0 .. y(prime)), "="
yspan=yspan(1:2);
y=[yspan(1);mean(yspan);yspan(2)];
a=sqrt(-[4 -8 4]*funca(y)*2/diff(yspan)^2);
b=-[-1 1]*funca(yspan(:))/diff(yspan)*2/a-a*y(2);
ys=[yspan(1):diff(yspan)/10:yspan(2)];
c=funca(ys)+(a*ys+b).^2/4;
if std(c)>eps^.4*(1+abs(mean(c))); warning('Int_exp: funca is not quadratic, NaN issued'); a=NaN; end;
c=mean(c);
% mapletry=['evalf(int(exp(',mat2str(c),'-(',mat2str(a),'*y+(',mat2str(b),'))^2/4),y=',mat2str(1),'..',mat2str(2),'));']
% maple(mapletry)
% maple(strrep(mapletry,'(int(','(Int('))
p=0;% se retrouve dans b.
res=(-exp_a_fois_erfc_b(p.*b./a+p.^2./a.^2+c, 1./2.*a.*yspan+1./2.*b+p./a))./a.*pi.^(1./2);
% res=(erfc_complex(1./2.*b+p./a)-erfc_complex(1./2.*a.*yspan+1./2.*b+p./a))./a./exp(-p.*b./a-p.^2./a.^2-c).*pi.^(1./2) faisait des NaN
res=diff(res);
% funca % (control value)
%[1.003,res(:)',a,b,c,p,p.*b./a+p.^2./a.^2+c, 1./2.*a.*yspan+1./2.*b+p./a,exp_a_fois_erfc_b(p.*b./a+p.^2./a.^2+c, 1./2.*a.*yspan+1./2.*b+p./a),exp(p.*b./a+p.^2./a.^2+c),erfc_complex( 1./2.*a.*yspan+1./2.*b+p./a)]

function res=Int_exp_a_fois_erfc_b(funca,funcb,yspan)
  % calcule Int(exp(funca(y))*erfc(funcb(y)),y=yspan(1)..yspan(2));
  % suppose que funca est de la forme d-p*y
  % suppose que funcb est de la forme (a*y+b)/2
% EQ48
%(((Int(erfc(1./2.*a.*y+1./2.*b)./exp(p.*y),y = 0 .. y(prime)), "=", ((
yspan=yspan(1:2);
y=[yspan(1);mean(yspan);yspan(2)];
p=-[-1 1]*funca(yspan(:))/diff(yspan);
ys=[yspan(1):diff(yspan)/10:yspan(2)];
d=funca(ys)+p*ys;
if std(d)>eps^.4*(1+abs(mean(d))); warning('Int_exp_a_fois_erfc_b: funca is not linear, NaN issued'); p=NaN; end;
d=mean(d);
a=[-1 1]*funcb(yspan(:))/diff(yspan)*2;
b=funcb(ys)*2-a*ys;
if std(b)>eps^.4*(1+abs(mean(b))); warning('Int_exp_a_fois_erfc_b: funcb is not linear, NaN issued'); a=NaN; end;
b=mean(b);
% mapletry=['evalf(int(exp(',mat2str(d),'-(',mat2str(p),')*y)*erfc((',mat2str(a),'*y+(',mat2str(b),'))/2),y=',mat2str(1),'..',mat2str(2),'));']
% maple(mapletry)
% maple(strrep(mapletry,'(int(','(Int('))
if abs(p)>eps^.8;
  res=-erfc_complex(1./2.*a.*yspan+1./2.*b)./p./exp(p.*yspan)+exp_a_fois_erfc_b(p.*b./a+p.^2./a.^2, 1./2.*a.*yspan+1./2.*b+p./a)./p;% pour éviter des NaN
  % res=erfc_complex(1./2.*b)./p-erfc_complex(1./2.*a.*yspan+1./2.*b)./p./exp(p.*yspan)+exp(p.*b./a+p.^2./a.^2).*erfc_complex(1./2.*a.*yspan+1./2.*b+p./a)-exp(p.*b./a+p.^2./a.^2).*erfc_complex(1./2.*b+p./a))./p; fait des NaN de temps en temps.
else
  res=(a.*yspan+b)./a.*erfc_complex(1./2.*a.*yspan+1./2.*b)+2.*(-exp(-1./4.*(a.*yspan+b).^2))./a./pi.^(1./2);
end
res=diff(res).*exp(d);
%[1.004,res(:)']

function res=exp_a_fois_erfc_b(a,b)
[erfcb,erfcxb]=erfc_complex(b);
[mini,lequel]=min(abs(log([abs(erfcb),abs(erfcxb)]))+[abs(a),abs(a-b.^2)],[],2);
res=exp(a-b.^2).*erfcxb;
res(lequel==1)=exp(a(lequel==1)).*erfcb(lequel==1);
%mat2str([1.005,res(:).'])

function differfcb=diff_erfc_plus_diff_diff_erfc(b,h,h1,h2)
differfcb=diff_erfc(b)*h.'+h1*diff_diff_erfc(b)*h2.';

function differfcb=diff_diff_erfc(b)% vérifié par maple evalf(D(D(erfc))(2)) ; ans-diff_diff_erfc(2): qui donne 0 !!
differfcb=-2/sqrt(pi).*exp(-b.^2)*(-2*b);
function res=moins_2_h2_exp_a_moins_carre_b_ab_sur_sqrt_pi(a,b,h)
  res=-2*h(2).*exp(a-b.^2)/sqrt(pi);

function res=constante_moins_2_h2_exp_a_moins_carre_b_ab_sur_sqrt_pi(a,b,h)
  res=-2*h(2)/sqrt(pi);


function differfcb=diff_erfc(b)
differfcb=-2/sqrt(pi).*exp(-b.^2);

function differfcb=constante_diff_erfc(b)
differfcb=-2/sqrt(pi);

function res=exp_a_fois_erfc_b_ab_plus_ab_ab(a,b,h,h1,h2)% h=diff(diff([a b],Fo),x), h1=diff([a b],Fo), h2=diff([a b],x)
res=exp_a_fois_erfc_b_ab(a,b,h)+exp_a_fois_erfc_b_ab_ab(a,b,h1,h2);

function res=exp_a_fois_erfc_b_ab_plus_ab_ab_ab(a,b,h,h1,h2,hh1,hh2,hh12)% h=diff(diff([a b],Fo),x), h1=diff([a b],Fo), h2=diff([a b],x), hh1= TODO
res=exp_a_fois_erfc_b_ab(a,b,h)+exp_a_fois_erfc_b_ab_ab(a,b,h1,h2);

function res_ab_ab=exp_a_fois_erfc_b_ab_ab(a,b,h1,h2)% vérifié par maple evalf(D(D(exp*erfc))(2)) ; mat2str(ans-exp_a_fois_erfc_b_ab_ab(2,2,[1 1],[1 1])) ; maple evalf(D(D(a->exp(a)*erfc(2*a-2)))(2)) ; mat2str(ans-exp_a_fois_erfc_b_ab_ab(2,2,[1 2],[1 2])) ; qui donne 0 puis 2.2204e-16
[erfcb,erfcxb]=erfc_complex(b); % >> maple diff(erfc(x),x)       % -2/pi^(1/2)*exp(-x^2)
differfcb=diff_erfc(b);
diffdifferfcb=diff_diff_erfc(b);
[mini,lequel]=min(abs(log([max(max(abs(erfcb))),max(max(abs(erfcxb)))]))+[max(max(abs(a))),max(max(abs(a-b.^2)))]);
if lequel==1
  res=exp(a).*erfcb;
else
  res=exp(a-b.^2).*erfcxb;
end
res_ab=[res,exp(a).*diff_erfc(b)];
res_ab_ab=h1*[res_ab;res_ab(2),exp(a).*diff_diff_erfc(b)]*h2.';
%[1.005,res(:)']

function res_ab=exp_a_fois_erfc_b_ab(a,b,h)
 % bug pour   %      a := 0.24999999999999939681755842359e-5  - 0.9682458365518550e-5   i; b := 0.0025 - 0.00193649167310371 i;h := [2.4999999999999939681755842359 - 9.682458365518550 i, 1250.0000000000000000000000000000 - 968.24583655185500000000000000000 i]
 % bug pour   %      a = 0.24999999999999939681755842359e-5  - 0.9682458365518550e-5i; b = 0.0025 - 0.00193649167310371i;h = [2.4999999999999939681755842359 - 9.682458365518550i, 1250.0000000000000000000000000000 - 968.24583655185500000000000000000i]
 % bug: le resultat est % -1408.00206853651-i*1082.88782304998 alors que la valeur attendue est -1407.9599413950666913846195667893 + 1082.8987002311311997155068782810i
[erfcb,erfcxb]=erfc_complex(b); % >> maple diff(erfc(x),x)       % -2/pi^(1/2)*exp(-x^2) % [0.997179047380717+i*0.00218508593569439,0.997181561441765+i*0.00217543622956666], et maple[erfc(b),erfc(b)*exp(b^2)] donne [0.99717904738071661750784943077620 + 0.0021850859356943913200305063897355i, 0.99718156144176486139621141916746 + 0.0021754362295666607243687359768916i], OK.
[mini,lequel]=min(abs(log([max(max(abs(erfcb))),max(max(abs(erfcxb)))]))+[max(max(abs(a))),max(max(abs(a-b.^2)))]);
if lequel==1
  res=exp(a).*erfcb;
else
  res=exp(a-b.^2).*erfcxb; % 0.997181561441765+i*0.00217543622956666 à comparer à maple exp(a)*erfc(b) = 0.99718156144176486139621141916746 + 0.0021754362295666607243687359768916i, OK !
end
% maple D(erfc)(b) % maple D(erfc)(b)
res_ab=[res,exp(a).*diff_erfc(b)]*h.'; % [0.997181561441765+i*0.00217543622956666,-1.12837916709551]*[2.49999999999999-i*9.68245836551855,1250-i*968.245836551855]', maple [exp(a)*erfc(b),exp(a)*D(erfc)(b)]*h qui donne [0.99718156144176486139621141916746 + 0.0021754362295666607243687359768916i , -1.1283791670955125738961589031215 + 0.64965818319026472977863439621835e-36i]*[ 2.4999999999999939681755842359 - 9.682458365518550i, 1250.0000000000000000000000000000 - 968.24583655185500000000000000000i]' BUG ici, ' pour matlab veut dire transposition ET conjugaison :-(
%[1.005,res(:).']

function res=Sum(func,span)% sum prenant un function_handle à la maple ; je ne sais pas faire ça avec one-liner ... sauf à créer le gros tableau :-(
res=func(1);
%ares=res;% a commenter une fois le debug fini.
%bres=ares;% a commenter une fois le debug fini.
for y=span(2:end)
  res=res+func(y);
  %bres=[bres,res-ares];% a commenter une fois le debug fini.
  %ares=res;% a commenter une fois le debug fini.
end
%func% a commenter une fois le debug fini.
%mat2str([1.006,span(end),res,bres])% a commenter une fois le debug fini.
% res=sum(arrayfun(func,span){:})

function res=majorant_logerfc(z)% majorant de erfc
res=max(imag(z).^2-real(z).^2,log(1-sign(real(z))));

function res=check_Int(varargin)
res=[erfc_complex(1:2),erfc_complex(1+i),erfc(1:2),Int_exp_a_fois_erfc_b(@(y)2*y+1,@(y)3*y-4,1:2),Int_exp(@(y)2*y+1-3*y.^2,1:2),Int_exp_a_fois_erfc_b(@(y)(2+i)*y+1-i,@(y)(3+2*i)*y-4-i,1:2),Int_exp(@(y)(2+i)*y+1-i-(3+i)*y.^2,1:2),Int_exp_a_fois_erfc_b(@(y)-2*y+1,@(y)1+y/10000,1:2),Int_exp(@(y)-y.^2/10000-2*y+1,1:2)*erfc_complex(1)];
mat2str(res)
a=eval(maple('evalf([erfc(1),erfc(2),erfc(1+sqrt(-1)),erfc(1),erfc(2),Int(exp(2*y+1)*erfc(3*y-4),y=1..2),Int(exp(2*y+1-3*y^2),y=1..2),Int(exp((2+sqrt(-1))*y+1-sqrt(-1))*erfc((3+2*sqrt(-1))*y-4-sqrt(-1)),y=1..2),Int(exp((2+sqrt(-1))*y+1-sqrt(-1)-(3+sqrt(-1))*y^2),y=1..2),Int(exp(-2*y+1)*erfc(1+y/10000),y=1..2),Int(exp(y^2/10000-2*y+1),y=1..2)*erfc(1)])'))
disp(mat2str(res-a))
[c_x_t_EQ46,errEQ46]=sensca_shorttime(.9,1,1.01,10,1,1,0,1)
[c_x_t_EQ46,errEQ46]=sensca_shorttime(.1,1,1.01,10,1,1,0,1)
%i;// :noremap à :wa<c-m>:!liste='sensc_compare_long_court.m sensc.m fi.m erfc_complex.m fig1.m fig2.m fig3.m fig4.m fig5.m fig6.m fig7.m cummax.m compare.m senscnonlin_v68.m fig_publiv1.m kronenlogs.m fig9.m sensca_shorttime.m' ; cvs status $liste ; cvs commit -m sauver ; scp -r $liste 10.75.4.18:senscrep ; scp -r $liste 10.75.4.24:senscrep ; scp -r $liste 10.75.4.19:senscrep ; scp -r $liste 10.75.4.18:ComsolDVD/hidden/senscrep <c-m><c-m>
% -rwxr-xr-x 1 goujot goujot   8899 Sep 10 13:59 pdesolve1D.m
% -rwxr-xr-x 1 goujot goujot    301 Sep 10 13:59 pdeloss.m
% -rwxr-xr-x 1 goujot goujot    615 Sep 10 13:59 jetcolormap12.m
% -rw-r--r-- 1 goujot goujot    449 Sep 10 13:59 kronenlogs.m


