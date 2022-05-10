clc;clear;
coef=[2,1,1,4;1,3,2,6;1,2,2,5];
coef1=[0,1,1,4;1,3,2,6;1,2,2,5];
coef2=[2,2,3,3;4,7,7,1;-2,4,5,-7];
coef3=[1,1,-1;2,1,0;1,-1,0];

[x]=gaosi(coef);[y]=liezhu(coef1);[z]=yuedang(coef2);[w]=qiuni(coef3);

disp(['高斯消去法求解方程组的根为：x=',num2str(x(1)),',y=',num2str(x(2)),',z=',num2str(x(3))]);
disp(['列主元高斯消去法求解方程组的根为：x=',num2str(y(1)),',y=',num2str(y(2)),',z=',num2str(y(3))]);
disp(['高斯-约当消去法求解方程组的根为：x=',num2str(z(1)),',y=',num2str(z(2)),',z=',num2str(z(3))]);
disp('用归一化的高斯-约当消去法求得逆矩阵为：');disp(w);

 function [x]=qiuni(coef)
    [r,~]=size(coef);
    x=eye(r);
    for i=1:r
        [~,ti]=max(coef(i:r,i));
        coef([ti+i-1,i],:)=coef([i,ti+i-1],:);
        x([ti+i-1,i],:)=x([i,ti+i-1],:);
       for j=1:r
           if j==i
               a=0;
           else
               a=coef(j,i)/coef(i,i);
           end
           coef(j,:)=coef(j,:)-a*coef(i,:);
           x(j,:)=x(j,:)-a*x(i,:);
       end
       x(i,:)=x(i,:)/coef(i,i);
       coef(i,:)=coef(i,:)/coef(i,i);
    end
 end

function [x]=yuedang(coef)
    [r,c]=size(coef);
    for i=1:r
       for j=1:r
          if j==i
             a=0;
          else
             a=coef(j,i)/coef(i,i);
          end
          coef(j,:)=coef(j,:)-a*coef(i,:);
       end
    end
    x=zeros(r,1);
    for i=1:r
       x(i,1)=coef(i,c)/coef(i,i); 
    end
end

function [x]=liezhu(coef)
    [r,c]=size(coef);
    for i=1:r
        [~,ti]=max(coef(i:r,i));
        coef([ti+i-1,i],:)=coef([i,ti+i-1],:);
        for j=i+1:r
           a=coef(j,i)/coef(i,i);
           coef(j,:)=coef(j,:)-a*coef(i,:);
        end
    end
    x=zeros(r,1);
    for i=r:-1:1
       for j=r:-1:1
           coef(i,c)=coef(i,c)-x(j,1)*coef(i,j);
       end
       x(i,1)=coef(i,c)/coef(i,i);
    end
end

 function [x]=gaosi(coef)
    [r,c]=size(coef);
    for i=1:r-1
       for j=i+1:r
           a=coef(j,i)/coef(i,i);
           coef(j,:)=coef(j,:)-a*coef(i,:);
       end
    end
    x=zeros(r,1);
    for i=r:-1:1
       for j=r:-1:1
           coef(i,c)=coef(i,c)-x(j,1)*coef(i,j);
       end
       x(i,1)=coef(i,c)/coef(i,i);
    end
 end
