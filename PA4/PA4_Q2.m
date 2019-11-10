format long
h=input('Enter h = ');
x = [0: h: 2];
str2=menu('Choose implementation for small end boundary','2nd order Backward Difference','2nd order Central Difference');
if str2==2
  x=[x,2+h];
end 
n=length(x);
p=ones(1,n);
q=-(x+3)./(x+1);
r=((x+3)./(x+1))./(x+1);
s=2.*(x+1)+3.*r;
p=p(2:n-1);
q=q(2:n-1);
r=r(2:n-1);
b=s(2:n-1);
l=(p./h)./h-q./(2*h);
d=-((2.*p)./h)./h+r;
u=(p./h)./h+q./(2*h);
b(1)=b(1)-5*l(1);
l(1)=0;
if str2==1
  l(n-2)=l(n-2)-u(n-2)/3;
  d(n-2)=d(n-2)+4/3*u(n-2);
  u(n-2)=0;
  n=n-2;
else
  l(n-2)=l(n-2)+u(n-2);
  u(n-2)=0;
  n=n-2;
end  
n;
length(l);
alpha=zeros(1,n);
beta=zeros(1,n);
y=zeros(1,n);
i=0;
alpha(1)=d(1);
beta(1)=b(1);
for i=[2:n]
    alpha(i)=d(i)-l(i)/alpha(i-1)*u(i-1);
    beta(i)=b(i)-l(i)/alpha(i-1)*beta(i-1);
end
y(n)=beta(n)/alpha(n);
for i = 2:n
    y(n-i+1)=(beta(n-i+1)-u(n-i+1)*y(n-i+2))/alpha(n-i+1);
end
if str2==1
  y=[5,y,4/3*y(n)-1/3*y(n-1)];
else
  y=[5,y];
end  

if str2==1
   plot(x,y,'-o','MarkerEdgeColor','r');
    title('Discretized Temperature');
    xlabel('Node Number');
    ylabel('Temperature');
    legend;
    hold on;
    fID = fopen('output2.txt','wt');
    fprintf(fID,'x            y\n');
    for i=1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
else
    length(x);
    length(y);
    plot(x(1:length(y)),y,'-o','MarkerEdgeColor','r');
    title('Discretized Temperature');
    xlabel('Node Number');
    ylabel('Temperature');
    legend;
    hold on;
    fID = fopen('output2.txt','wt');
    fprintf(fID,'x            y\n');
    for i=1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
end 