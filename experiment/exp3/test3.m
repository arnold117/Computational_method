clear;clc;
coef=[1,0,-1,-1];coef1=[1,0,0,0];
[a,b]=diedai(coef,1.5);c=[a,b];[d,e]=niudun(coef,1.5);f=[d,e];g=ones(2,11);
[k,l]=xiashan(coef,0.4,0.1);m=[k,l];
for i=1:11
    coef1(4)=-(i+100);
    g(1,i)=niudun(coef1,fix((i+100)^(1/3)));
end
for i=1:11
    coef1(4)=-(i+100);
    g(2,i)=jianniu(coef1,fix((i+100)^(1/3)));
end
disp(['�򵥵�������ø���',num2str(c(1)),',����',num2str(c(2)),'��']);
disp(['ţ�ٵ�������ø���',num2str(f(1)),',����',num2str(f(2)),'��']);
disp(num2str(g));
disp(['ţ����ɽ����ø���',num2str(m(1)),',����',num2str(m(2)),'��']);
% disp(sym2poly(diff(poly2sym(coef))));
% coeff=-coef(2:1:end);
% disp(polyval(coeff,x)^(1/(length(coef)-1)));
% disp(2^(1/3));
% disp(coeff);
% disp(poly2sym(coef));
% disp(polyval(coef,x));
% disp(length(coef)-1);
function [fx,a]=diedai(coef,x)
a=0;coeff=-coef(2:1:end);
fx=polyval(coeff,x)^(1/(length(coef)-1));
while true
   fxx=polyval(coeff,fx)^(1/(length(coef)-1));
   if fx==fxx
       break;
   else
       fx=fxx;
       a=a+1;
   end
end
end
function [fx,a]=niudun(coef,x)
coeff=sym2poly(diff(poly2sym(coef)));
a=0;fx=x-(polyval(coef,x))/(polyval(coeff,x));
while true
   fxx=fx-(polyval(coef,fx))/(polyval(coeff,fx));
   if fx==fxx
       break;
   else
       fx=fxx;
       a=a+1;
   end
end
end
function fx=jianniu(coef,x)
coeff=sym2poly(diff(poly2sym(coef)));
c=polyval(coeff,x);
fx=x-polyval(coef,x)/c;
while true
    c=polyval(coeff,fx);
   fxx=fx-polyval(coef,fx)/c;
   if fx==fxx
       break;
   else
       fx=fxx;
   end
end
end
function [c,b]=xiashan(coef,x,B)
coeff=sym2poly(diff(poly2sym(coef)));
fx=x-(polyval(coef,x))/(polyval(coeff,x));A=1;
while true
   fxx=fx-A*(polyval(coef,fx))/(polyval(coeff,fx));
   if abs(fx)>abs(fxx)
       [c,b]=niudun(coef,fxx);
       break;
   elseif A>B
       A=A/2;
   else
       disp('������ѡ��');
       break;
   end
end
end