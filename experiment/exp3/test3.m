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
disp(['简单迭代法求得根：',num2str(c(1)),',迭代',num2str(c(2)),'次']);
disp(['牛顿迭代法求得根：',num2str(f(1)),',迭代',num2str(f(2)),'次']);
disp(num2str(g));
disp(['牛顿下山法求得根：',num2str(m(1)),',迭代',num2str(m(2)),'次']);
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
       disp('请重新选根');
       break;
   end
end
end