%
% John B. Taylor Multicountry Rational Expectations Model
%	Linearized 12/10/96 by John C. Williams
%
%	Parameter File
%
% ----------------------------------------------------
% Country codes
%	0:	United States
%	1:	Canada
%	2:	France
%	3:	Germany
%	4:	Italy
%	5:	Japan
%	6:	United Kingdom

% ----------------------------------------------------
% tayrulex:	= 1 : Taylor policy rule
%		= 0 : Money supply policy rule

tayrule0=	1;
tayrule1=	0;
tayrule2=	0;
tayrule3=	0;
tayrule4=	0;
tayrule5=	0;
tayrule6=	0;


% ----------------------------------------------------
% nomrulex:	= 1 : price level policy rule
%		= 0 : other policy rule

nomrule0=	0;
nomrule1=	0;
nomrule2=	0;
nomrule3=	0;
nomrule4=	0;
nomrule5=	0;
nomrule6=	0;

% -----------------------------------------------------
% Taylor Rule Policy parameters
%	mrulerx	:	lag of short rate
%	mrulepx	:	current four-quarter inflation rate
%	mrulep1x:	lagged four-quarter inflation rate
%	mruleyx	:	current output gap 
%	mruley1x:	lagged output gap 

mruler0	=	0.84;
mrulep0	=	0.25;
mrulep10=	0.;
mruley0	=	0.5;
mruley10=	0.;
mrulep1	=	0.5;
mruley1	=	0.5;
mrulep2	=	0.5;
mruley2	=	0.5;
mrulep3	=	0.5;
mruley3	=	0.5;
mrulep4	=	0.5;
mruley4	=	0.5;
mrulep5	=	0.5;
mruley5	=	0.5;
mrulep6	=	0.5;
mruley6	=	0.5;

% -----------------------------------------------------
% rholmx: autoregressive parameter in money growth 

rholm0	=	0;
rholm1	=	0;
rholm2	=	0;
rholm3	=	0;
rholm4	=	0;
rholm5	=	0;
rholm6	=	0;

% -----------------------------------------------------
% rholgx: autoregressive parameter in government spending

rhog0	=	0;
qtr4g0	=	0; % = 1 for four-quarter shock (set rhog0=1)
rhog1	=	0;
rhog2	=	0;
rhog3	=	0;
rhog4	=	0;
rhog5	=	0;
rhog6	=	0;

% -----------------------------------------------------
% normalization factors for demand components 
%	set approximately to size of GDP

dummy0=1/3300;
dummy1=1/400;
dummy2=1/1200;
dummy3=1/1500;
dummy4=1/86600;
dummy5=1/267300;
dummy6=1/200;

% -----------------------------------------------------
% linearization parameters
%	computed by `taylor88lindata.m'

cdbar0	=        287.44112715*dummy0;
cnbar0	=        800.58688815*dummy0;
csbar0	=       1066.20013860*dummy0;
cdbar1	=         29.34357465*dummy1;
cnbar1	=         81.82718067*dummy1;
csbar1	=         90.21667857*dummy1;
cdbar2	=         76.29188942*dummy2;
cnbar2	=        380.77972252*dummy2;
csbar2	=        299.32193825*dummy2;
cbar3	=        842.89099138*dummy3;
cbar4	=      55370.13193467*dummy4;
cdbar5	=       9978.56301367*dummy5;
cnbar5	=      67322.78207495*dummy5;
csbar5	=      75625.26508231*dummy5;
cdbar6	=         15.70737838*dummy6;
cnbar6	=         74.32376637*dummy6;
csbar6	=         53.33458143*dummy6;
inebar0	=        256.61155891*dummy0;
insbar0	=        142.23620642*dummy0;
irbar0	=        142.81592778*dummy0;
iibar0	=          5.13028538*dummy0;
ifbar1	=         78.95572631*dummy1;
iibar1	=          1.41985544*dummy1;
inbar2	=        185.05576929*dummy2;
irbar2	=         54.21744542*dummy2;
iibar2	=          3.79803453*dummy2;
ifbar3	=        311.63836051*dummy3;
iibar3	=          1.14713037*dummy3;
ifbar4	=      14619.81992843*dummy4;
iibar4	=        269.74273026*dummy4;
inbar5	=      43864.56074645*dummy5;
irbar5	=      14255.95039353*dummy5;
iibar5	=        860.44831732*dummy5;
inbar6	=         34.02379273*dummy6;
irbar6	=          7.99357557*dummy6;
iibar6	=          0.15344138*dummy6;
exbar0	=        366.57605035*dummy0;
exbar1	=        106.15320313*dummy1;
exbar2	=        283.58817058*dummy2;
exbar3	=        488.57738639*dummy3;
exbar4	=      21935.80439685*dummy4;
exbar5	=      46853.35724714*dummy5;
exbar6	=         66.43780918*dummy6;
imbar0	=        389.81939166*dummy0;
imbar1	=         92.03504823*dummy1;
imbar2	=        286.08852761*dummy2;
imbar3	=        434.80583651*dummy3;
imbar4	=      18784.41171314*dummy4;
imbar5	=      40948.51001004*dummy5;
imbar6	=         63.29733797*dummy6;
ybar0	=       3349.96723345*dummy0;
ybar1	=        364.62165084*dummy1;
ybar2	=       1154.29388235*dummy2;
ybar3	=       1512.17692015*dummy3;
ybar4	=      86626.18331931*dummy4;
ybar5	=     267273.91883793*dummy5;
ybar6	=        238.24320511*dummy6;

% -----------------------------------------------------
% trend term for potential output
%	
% t = 1 for 1971:1

t =	50; % t = 1 for 1971:1

yt0 = exp(0* 7.824554 + 0.006021*t)*dummy0;
yt1 = exp(0* 5.513943 + 0.007953*t)*dummy1;
yt2 = exp(0* 6.755867 + 0.006074*t)*dummy2;
yt3 = exp(0* 7.074856 + 0.005077*t)*dummy3;
yt4 = exp(0*11.094975 + 0.005571*t)*dummy4;
yt5 = exp(0*11.977690 + 0.010259*t)*dummy5;
yt6 = exp(0* 5.299993 + 0.003728*t)*dummy6;
