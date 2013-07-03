%
% John B. Taylor Multicountry Rational Expectations Model
%	Linearized 12/10/96 by John C. Williams
%
%	Impulse response simulation file
%
%	This program generates impulse responses to specified
%	shocks, and computes percent deviations from baseline
%	for demand components
%

nsh=112;		% # of stochastic equations
nbeq=size(cofb,1)-nsh;	% # of model variables
neq=nbeq+nsh;		% # of variables + shocks
nlags=4;		% maximum # of lags
np = 50;		% # of data points
nt=np-nlags;		% length of simulation

%
% Compute VAR representation from AIM solution matrices 
%   (scof and cofb)
%

amat1=cofb(1:nbeq,3*neq+1:3*neq+nbeq);
amat2=cofb(1:nbeq,2*neq+1:2*neq+nbeq);
amat3=cofb(1:nbeq,1*neq+1:1*neq+nbeq);
amat4=cofb(1:nbeq,0*neq+1:0*neq+nbeq);

bmat=inv(scof(1:neq,nlags*neq+1:(nlags+1)*neq));
bmat=bmat(1:nbeq,nbeq+1:neq);

%
% Initialize model variables and shocks to zero 
%    (deviation from baseline)
%

x=zeros(np,nbeq);
e=zeros(np,nsh);

%
% Set shock(s) to system
%

%e(nlags+1,1)	=0.01;		% rs0
%e(nlags+1,2)	=0.01;		% rs1
%e(nlags+1,3)	=0.01;		% rs2
%e(nlags+1,4)	=0.01;		% rs3
%e(nlags+1,5)	=0.01;		% rs4
%e(nlags+1,6)	=0.01;		% rs5
%e(nlags+1,7)	=0.01;		% rs6

%e(nlags+1,8)	=0.01;		% le1
%e(nlags+1,9)	=0.01;		% le2
%e(nlags+1,10)	=0.01;		% le3
%e(nlags+1,11)	=0.01;		% le4
%e(nlags+1,12)	=0.01;		% le5
%e(nlags+1,13)	=0.01;		% le6

%e(nlags+1,14)	=0.01;		% rl0
%e(nlags+1,15)	=0.01;		% rl1
%e(nlags+1,16)	=0.01;		% rl2
%e(nlags+1,17)	=0.01;		% rl3
%e(nlags+1,18)	=0.01;		% rl4
%e(nlags+1,19)	=0.01;		% rl5
%e(nlags+1,20)	=0.01;		% rl6

%e(nlags+1,21)  =0.01;          % cd0_
%e(nlags+1,22)  =0.01;          % cn0_
%e(nlags+1,23)  =0.01;          % cs0_
%e(nlags+1,24)  =0.01;          % cd1_
%e(nlags+1,25)  =0.01;          % cn1_
%e(nlags+1,26)  =0.01;          % cs1_
%e(nlags+1,27)  =0.01;          % cd2_
%e(nlags+1,28)  =0.01;          % cn2_
%e(nlags+1,29)  =0.01;          % cs2_
%e(nlags+1,30)  =0.01;          % c3_ 
%e(nlags+1,31)  =0.01;          % c4_ 
%e(nlags+1,32)  =0.01;          % cd5_
%e(nlags+1,33)  =0.01;          % cn5_
%e(nlags+1,34)  =0.01;          % cs5_
%e(nlags+1,35)  =0.01;          % cd6_
%e(nlags+1,36)  =0.01;          % cn6_
%e(nlags+1,37)  =0.01;          % cs6_

%e(nlags+1,38)  =0.01;          % ine0_
%e(nlags+1,39)  =0.01;          % ins0_
%e(nlags+1,40)  =0.01;          % ir0_ 
%e(nlags+1,41)  =0.01;          % ii0_ 
%e(nlags+1,42)  =0.01;          % if1_ 
%e(nlags+1,43)  =0.01;          % ii1_ 
%e(nlags+1,44)  =0.01;          % in2_ 
%e(nlags+1,45)  =0.01;          % ir2_ 
%e(nlags+1,46)  =0.01;          % ii2_ 
%e(nlags+1,47)  =0.01;          % if3_ 
%e(nlags+1,48)  =0.01;          % ii3_ 
%e(nlags+1,49)  =0.01;          % if4_ 
%e(nlags+1,50)  =0.01;          % ii4_ 
%e(nlags+1,51)  =0.01;          % in5_ 
%e(nlags+1,52)  =0.01;          % ir5_ 
%e(nlags+1,53)  =0.01;          % ii5_ 
%e(nlags+1,54)  =0.01;          % in6_ 
%e(nlags+1,55)  =0.01;          % ir6_ 
%e(nlags+1,56)  =0.01;          % ii6_ 
				  
%e(nlags+1,57)  =0.01;          % lex0_
%e(nlags+1,58)  =0.01;          % lex1_
%e(nlags+1,59)  =0.01;          % lex2_
%e(nlags+1,60)  =0.01;          % lex3_
%e(nlags+1,61)  =0.01;          % lex4_
%e(nlags+1,62)  =0.01;          % lex5_
%e(nlags+1,63)  =0.01;          % lex6_
      				  
%e(nlags+1,64)  =0.01;          % lim0_
%e(nlags+1,65)  =0.01;          % lim1_
%e(nlags+1,66)  =0.01;          % lim2_
%e(nlags+1,67)  =0.01;          % lim3_
%e(nlags+1,68)  =0.01;          % lim4_
%e(nlags+1,69)  =0.01;          % lim5_
%e(nlags+1,70)  =0.01;          % lim6_

%e(nlags+1,71)	=0.01;		% lx0
%e(nlags+1,72)	=0.01;		% lx1
%e(nlags+1,73)	=0.01;		% lx2
%e(nlags+1,74)	=0.01;		% lx3
%e(nlags+1,75)	=0.01;		% lx4
%e(nlags+1,76)	=0.01;		% lx5
%e(nlags+1,77)	=0.01;		% lx6

%e(nlags+1,78)	=0.01;		% lp0
%e(nlags+1,79)	=0.01;		% lp1
%e(nlags+1,80)	=0.01;		% lp2
%e(nlags+1,81)	=0.01;		% lp3
%e(nlags+1,82)	=0.01;		% lp4
%e(nlags+1,83)	=0.01;		% lp5
%e(nlags+1,84)	=0.01;		% lp6

%e(nlags+1,85)	=0.01;		% lpi0
%e(nlags+1,86)	=0.01;		% lpi1
%e(nlags+1,87)	=0.01;		% lpi2
%e(nlags+1,88)	=0.01;		% lpi3
%e(nlags+1,89)	=0.01;		% lpi4
%e(nlags+1,90)	=0.01;		% lpi5
%e(nlags+1,91)	=0.01;		% lpi6

%e(nlags+1,92)	=0.01;		% lpe0
%e(nlags+1,93)	=0.01;		% lpe1
%e(nlags+1,94)	=0.01;		% lpe2
%e(nlags+1,95)	=0.01;		% lpe3
%e(nlags+1,96)	=0.01;		% lpe4
%e(nlags+1,97)	=0.01;		% lpe5
%e(nlags+1,98)	=0.01;		% lpe6

%e(nlags+1,99)	=0.01*ybar0;	% g0
%e(nlags+1,100)	=0.01*ybar1;	% g1
%e(nlags+1,101)	=0.01*ybar2;	% g2
%e(nlags+1,102)	=0.01*ybar3;	% g3
%e(nlags+1,103)	=0.01*ybar4;	% g4
%e(nlags+1,104)	=0.01*ybar5;	% g5
%e(nlags+1,105)	=0.01*ybar6;	% g6

%e(nlags+1,106)	=0.03;		% lm0
%e(nlags+1,107)	=0.03;		% lm1
%e(nlags+1,108)	=0.03;		% lm2
%e(nlags+1,109)	=0.03;		% lm3
%e(nlags+1,110)	=0.03;		% lm4
%e(nlags+1,111)	=0.03;		% lm5
%e(nlags+1,112)	=0.03;		% lm6

% Run simulation

for i = nlags+1:np;
   x1=x(i-1,:)';
   x2=x(i-2,:)';
   x3=x(i-3,:)';
   x4=x(i-4,:)';
   ee=e(i,:)';
   xc=amat1*x1+amat2*x2+amat3*x3+amat4*x4+bmat*ee;
   x(i,:)=xc';
end;

%
% Plot US price level, output gap, funds rate, bond rate
%

t=0:40;
subplot(2,2,1);
plot(t,100*x(4:44,92));
title('Price Level');xlabel('Quarters');ylabel('Percentage Points');   
subplot(2,2,2); plot(t,100*x(4:44,71)/ybar0);
title('Output Gap');xlabel('Quarters');ylabel('Percentage Points'); 
subplot(2,2,3); plot(t,100*x(4:44,1));
title('Fed Funds Rate');xlabel('Quarters');ylabel('Percentage Points');
subplot(2,2,4);  plot(t,100*x(4:44,14));
title('10 Year Bond Rate');xlabel('Quarters');ylabel('Percentage Points');
