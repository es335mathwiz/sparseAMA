%paramsw.m:  Parameter file for linear model used in GLSS AER revision.

%clear all;

%---------------------------Parameter values--------------------------------%

beta  = 0.9987;
pibar = 1.006;
Gz = 1.0041;  %tech growth (gross)
psil = 1;
gamma = 0.858;  %habit persistence
sigmal = 4.49;   %governs labor supply elasticity
phi = 95;
phiw = 8000;  
ep = 6;
epw = 8;
ap = 0.87;    %(1-ap) is the degree of backward indexation of prices
aw = 0.92;    %(1-aw) is the degree to which wages are indexed to price inflation
bw = 0.92;   %(1-bw) is the degree to which wages are indexed to tech shock

intswitch = 1;  % 1 internal investment adjustment cost  0, external costs 
ihabitswitch = 1;  % 1 internal habits; 0 external 
hpswitch = 0;   % 0 one-sided hp-filter  1 two-sided hp-filter
lamhp = 1600;

markup = ep/(ep-1);  %steady state markup
markupw = epw/(epw-1);
alpha = 0.167;  %capital elasticity in C-D production function
delta = 0.025;
phii = 3.14;    %adj. costs on investment (external)
sigmaa = 10000;

%Taylor rule parameters (This is our approximation to JPT)
gam_rs = 0.86;
gam_dp = 1.688;
gamxhp = 0;
gamdy = 0.21/(1-gam_rs);

shrgy = 0.2;

%---------------------------------shock parameters--------------------------------%

%technology shock
rhotech = 0;
sdevtech = 0.01;

%gov't shock
rhog = 0.95;
sdevg = 0.01;

%inv specific shock
rhoinv = 0.77;  
sdevinv = 0.07;

%eta shock
rhoeta = 0.9;  
sdeveta = 0.01;

%monetary shock
rhoint = 0;  
sdevint = 0.01;

%markup shock
rhomark = 0.98;  
sdevmark = 0.01;

%wage markup shock
rhomarkw = 0.98;
sdevmarkw = 0.01;

% MUC shock
rhomuc1 = 1.4;
rhomuc2 = 1 - rhomuc1 - 0.001;
sdevmuc = 0.1;


%--------------------------------------------------%
%Use to overwrite above parameters to compare different parameters of the
%model.
%----------------------------------------------------%
changeparmlist = []; 
changeparmval = [];
nparmchange = size(changeparmlist,1);
for indxi=1:nparmchange
   eval([changeparmlist(indxi,:),'=changeparmval(indxi);'])
end

% %---------------------------------------------------------
% %Steady state and definitions used by linearized model.
% %
% %-----------------------------------------------------------

gamtil = gamma/Gz;
realrate = Gz/beta;
mc = 1/markup;
k2yrat = ((mc*alpha)/(realrate-(1-delta)))*Gz;
shriy = (1-(1-delta)/Gz)*k2yrat;
gg = 1/(1-shrgy);
shrcy = 1-shriy-shrgy;
if (ihabitswitch == 0)
    labss = ((1/markupw)*(1-alpha)*mc*(1/(psil*(1-gamtil)))*(1/shrcy))^(1/(sigmal+1));
else
    labss = ((1/markupw)*(1-alpha)*(1-beta*gamtil)*mc*(1/(psil*(1-gamtil)))*(1/shrcy))^(1/(sigmal+1));
end
kss = labss*(Gz^(alpha/(alpha-1)))*k2yrat^(1/(1-alpha));
gdpss = (kss/Gz)^alpha*labss^(1-alpha);
invss = shriy*gdpss;
css = shrcy*gdpss;
rwss = (1-alpha)*mc*gdpss/labss;
mucss = (1/css)*(1/(1-gamtil));
if (ihabitswitch == 0)
    lamss = mucss;
else
    lamss = mucss*(1-beta*gamtil);
end
kappap = (ep-1)/(phi*(1+beta*(1-ap)));
kappaw = epw*(1-gamtil)*psil*labss^(1+sigmal)/phiw;






