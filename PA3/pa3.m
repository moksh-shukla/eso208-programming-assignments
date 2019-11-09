format long
%x=input('Enter the vector x = ');
%y=input('Enter the vector y = ');
prompt = ('What do you want to do\n A-Fit a least square polynomial\n B-Fit a Lagrange Interpolation Polynomial\n C-Fit a Newton’s Interpolation Polynomial \n D-Fit Cubic splines\n ');
str = input(prompt,'s');
if strcmp(str, 'A')
    N=input('Enter the order of polynomial (an integer <n) \n');
    n=N+1;
    %answer = inputdlg('Input file name');
    file = input('Give me the name of the input:','s');
   M = readmatrix(file);
    x = M(:,1);
    y = M(:,2);
    A=zeros(n);
    b=zeros(n,1);
    for i = 1:n
      b(i)=y'*(x.^(i-1));
      for j = i:n
        A(i,j)=transpose(x.^(j-1))*x.^(i-1);
        A(j,i)=A(i,j);
      end
    end
    
    i=0; 
    j=0;
    k=0;
    m=0;
    A;
    b;
    X=zeros(1,n);
    for k = 1:n-1
      for i = k+1:n
        m=A(i,k)/A(k,k);
        A(i,k)=0;
        b(i)=b(i)-m*b(k);
        for j = k+1:n
           A(i,j)=A(i,j)-m*A(k,j);
        end
        A;
      end
    end
    A;
    b;
    X(n)=b(n)/A(n,n);
    for i = n-1:-1:1
      m=0;
      for j = i+1:n 
        m=m+A(i,j)*X(j);
      end   
      X(i)=(b(i)-m)/A(i,i);
    end
    num_points = A(1,1);
    X=fliplr(X); 
    y_=polyval(X,x);
    mu= sum(y)/num_points;
    sigma=sum((y-mu).^2);
    eps=sum((y-y_).^2);
    r_square=1-eps/sigma;
    x1 = linspace(min(x),max(x));
    y1 = polyval(X,x1);
    p = poly2sym(X);
    coff = coeffs(p);
    plot(x,y,'o')
    hold on
    plot(x1,y1)
    hold off
    matrix = coeffs(p);
%     fileID = fopen('Least_Square.txt','w');
%     %fpri
%     fprintf(fileID,'%f\n', r_square);
%     
%     for j = 1:length(matrix)
%        fprintf(fileID,"%f\t",matrix(1,j));
%     end
%             
%       
%     fclose(fileID);
    n=n+1;
    A=zeros(n);
    b=zeros(n,1);
    for i = 1:n
      b(i)=y'*(x.^(i-1));
      for j = i:n
        A(i,j)=transpose(x.^(j-1))*x.^(i-1);
        A(j,i)=A(i,j);
      end
    end
    
    i=0; 
    j=0;
    k=0;
    m=0;
    A;
    b;
    X=zeros(1,n);
    for k = 1:n-1
      for i = k+1:n
        m=A(i,k)/A(k,k);
        A(i,k)=0;
        b(i)=b(i)-m*b(k);
        for j = k+1:n
           A(i,j)=A(i,j)-m*A(k,j);
        end
        A;
      end
    end
    A;
    b;
    X1(n)=b(n)/A(n,n);
    for i = n-1:-1:1
      m=0;
      for j = i+1:n 
        m=m+A(i,j)*X1(j);
      end   
      X1(i)=(b(i)-m)/A(i,i);
    end
    num_points = A(1,1);
    X1=fliplr(X1); 
    y1_=polyval(X1,x);
    mu= sum(y)/num_points;
    sigma=sum((y-mu).^2);
    eps=sum((y-y1_).^2);
    r_square1=1-eps/sigma;
    x2 = linspace(min(x),max(x));
    y2 = polyval(X,x2);
    p1 = poly2sym(X1);
    coff = coeffs(p1);
    %plot(x,y,'o')
    %hold on
    %plot(x1,y1)
    %hold off
    matrix1 = coeffs(p1);
    
    
    if (n-1>0)
        n=n-1;
        A=zeros(n);
        b=zeros(n,1);
        for i = 1:n
          b(i)=y'*(x.^(i-1));
          for j = i:n
            A(i,j)=transpose(x.^(j-1))*x.^(i-1);
            A(j,i)=A(i,j);
          end
        end

        i=0; 
        j=0;
        k=0;
        m=0;
        A;
        b;
        X=zeros(1,n);
        for k = 1:n-1
          for i = k+1:n
            m=A(i,k)/A(k,k);
            A(i,k)=0;
            b(i)=b(i)-m*b(k);
            for j = k+1:n
               A(i,j)=A(i,j)-m*A(k,j);
            end
            A;
          end
        end
        A;
        b;
        X2(n)=b(n)/A(n,n);
        for i = n-1:-1:1
          m=0;
          for j = i+1:n 
            m=m+A(i,j)*X2(j);
          end   
          X2(i)=(b(i)-m)/A(i,i);
        end
        num_points = A(1,1);
        X2=fliplr(X2); 
        y2_=polyval(X2,x);
        mu= sum(y)/num_points;
        sigma=sum((y-mu).^2);
        eps=sum((y-y2_).^2);
        r_square2=1-eps/sigma;
        x3 = linspace(min(x),max(x));
        y3 = polyval(X,x3);
        p2 = poly2sym(X2);
        coff = coeffs(p2);
        %plot(x,y,'o')
        %hold on
        %plot(x1,y1)
        %hold off       
        
        matrix2 = coeffs(p2);

        fileID = fopen('Least_Square.txt','w');
        fprintf(fileID,'For n\n');
        fprintf(fileID,'%f\n', r_square);

        for j = 1:length(matrix)
           fprintf(fileID,"%f\t",matrix(1,j));
        end
        fprintf("\n");
        fprintf(fileID,'For n+1\n');
        fprintf(fileID,'%f\n', r_square1);    
        for j = 1:length(matrix1)
           fprintf(fileID,"%f\t",matrix1(1,j));
        end 
        
        fprintf('\n');
        fprintf(fileID,'For n-1\n');
        fprintf(fileID,'%f\n', r_square2);    
        for j = 1:length(matrix2)
           fprintf(fileID,"%f\t",matrix2(1,j));
        end

        fclose(fileID);
    end
    
    
        fileID = fopen('Least_Square.txt','w');
        fprintf(fileID,'For n\n');
        fprintf(fileID,'%f\n', r_square);

        for j = 1:length(matrix)
           fprintf(fileID,"%f\t",matrix(1,j));
        end
        fprintf("\n");
        fprintf(fileID,'For n+1\n');
        fprintf(fileID,'%f\n', r_square1);    
        for j = 1:length(matrix1)
           fprintf(fileID,"%f\t",matrix1(1,j));
        end
        
        
    
    
end

if strcmp(str, 'B')
  file = input('Give me the name of the input:','s');
   M = readmatrix(file);
  x = M(:,1);
  y = M(:,2);
  h=[];
  k=0; 
  f=0;
  i=0;
  j=0;
  n = size(M,1);
  x1=x;
  g= sym('x');
  K=[];
  for i = 1:n
     h=[];  
     for j = 1:n-1 
        h = [h,x1(mod(i+j-1,n)+1)];  
     end 
     f=expand(prod(g-h)); 
     f;
     C = coeffs(f);
     K=zeros(1,n);
     b=0;
     for b = n:-1:n-length(C)+1
         K(b)=C(b-n+length(C));
     end    
     o=0;
     oo=0;
     %K;
     for oo = 1:n
         x1(i)^(oo-1);
         %K(oo);
         o=o+(x1(i)^(oo-1))*K(oo);
         o;
     end    
     %x1(i)
     %o=polyval(K,x1(i))
     o;
     y(i)*f/o;
     k = k+y(i)*f/o;
     p=k;
  end 
  disp(p);
  C = coeffs(k);
  K=zeros(1,n);
  b=0;
  for b = n:-1:n-length(C)+1
      K(b)=C(b-n+length(C));
  end 
  x2 = linspace(min(x1),max(x1));
  i=0;
  y2=zeros(1,length(x2));
  for i = 1:length(x2)
    oo=0;
    for oo = 1:n
       y2(i)=y2(i)+(x2(i)^(oo-1))*K(oo);
    end    
  end
  plot(x,y,'o')
  hold on
  plot(x2,y2)
  hold off
  matrix = coeffs(p);
  fileID = fopen('Lagrange.txt','w');
    for i = 1:length(matrix)
        fprintf(fileID,"%f ",matrix(1,i));
    end
   fclose(fileID);
   
end

if strcmp(str,'C')
    
   file = input('Give me the name of the input:','s');
   M = readmatrix(file);
   x = M(:,1);
   y = M(:,2);
   n = length(x);
   D = zeros(n,n);
   D(:,1) = y';
   
   for j=2:n
     for k=j:n
         D(k,j) = (D(k,j-1)-D(k-1,j-1))/(x(k)-x(k-j+1));
     end
   end
   
   C = D(n,n);
   for k=(n-1):-1:1
     C = conv(C,poly(x(k)));
     m = length(C);
     C(m) = C(m) + D(k,k);
   end
   
   m = sym('x');
   f = 0;
   for k = n-1:-1:1
       f = f + C(n-k)*(m^k);
       if(k==1)
           f = f + C(n);
       end
   end
   
   X = sym2poly(f);
   x1 = linspace(min(x),max(x));
   y1 = polyval(X,x1);
   plot(x,y,'o')
   hold on
   plot(x1,y1)
   hold off
   fileID = fopen('Newton.txt','w');
    for i = 1:length(X)
        fprintf(fileID,"%f ",X(1,i));
    end
   fclose(fileID);
   
end

if strcmp(str,'D')
   file = input('Give me the name of the input:','s');
   M = readmatrix(file);
   X = readmatrix("input_spline.txt");
   x = M(:,1);
   y = M(:,2);
   prompt2 = 'Options to choose Linear Spline (1), Quadratic Spline (2), Natural Spline (3), Not-a-knot (4), Periodic (5), Clamped Spline (6)';
  str2 = input(prompt2);
  %linear spline
  if str2 == 1
      l = length(x);
      a0 = zeros(1,l-1);
      a1 = zeros(1,l-1);
 
   for i = 1:l-1
        a = [a0(i);a1(i)];
        A = [1 x(i);1 x(i+1)];
        b = [y(i);y(i+1)];
        a = A\b;
        a0(i) = a(1,1);a1(i) = a(2,1);
   end
 
    file = input('Give me the name of the output file:','s');
    fid = fopen(file,'w');
    for i = 1:l-1
        fprintf(fid,'%f %f in [%f,%f]\n',a1(i),a0(i),x(i),x(i+1));
    end
    plot (x,y,'-*r');
    fclose(fid);

  end
      
  if str2 == 2
      l = length(x);
    u = zeros(1,l);
    u(1) = 0; %chosen condition
    for i = 2:l
        u(i) = 2*((y(i) - y(i-1))/(x(i) - x(i-1)))-u(i-1);
    end
    c = zeros(1,l-1);
    for i = 1:l-1
        c(i) = y(i)+ (u(i)*(x(i+1) - x(i)))/2;
    end
    % coefficients 
    a0 = zeros(1,l-1);a1 = zeros(1,l-1);a2 = zeros(1,l-1);
    for i = 1:l-1
        l2 = 1;
        for j = 1:2
            a = poly(x(i));
            l2 = conv(l2,a);
        end
        term1 = ((u(i+1)/(x(i+1) - x(i)))*l2)/2;
        l1 = 1;
        for j = 1:2
            a_1 = poly(x(i+1));
            l1 = conv(l1,a_1);
        end
        term2 = ((u(i)/(x(i+1) - x(i)))*l1)/2;
        p = term1 - term2;
        p(3) = p(3) + c(i);
        a0(i) = p(3);a1(i) = p(2);a2(i) = p(1);
    end
    file = input('Give me the name of the output:','s');
    fid = fopen(file,'w');
    l2 = l;
    for k = 1:l-1
        fprintf(fid,'%f %f %f in [%f,%f]\n',a2(k),a1(k),a0(k),x(k),x(k+1));
    end
    x1 = x(1):0.01:x(2);
    y1 = polyval([a2(1),a1(1),a0(1)],x1);
    plot (x1,y1);
    hold on
    plot (x,y,'*');
    for i = 2:l-1
        x1 = x(i):0.01:x(i+1);
        y1 = polyval([a2(i),a1(i),a0(i)],x1);
        plot(x1,y1);
    end
    hold off

  end
  
  %natural spline
  if str2 == 3
      n=length(x);
      %X = readMatrix("input_spline.txt");
h=zeros(n-1,1);
for i=1:n-1
    h(i,1)=x(i+1,1)-x(i,1);
end
a=zeros(n,n);
%boundary condition
a(1,1)=1;
a(n,n)=1;
for i=2:n-1
    a(i,i-1)=h(i-1,1);
    a(i,i)=2*(h(i-1,1)+h(i,1));
    a(i,i+1)=h(i,1);
end
b=zeros(n,1);
g=zeros(n,1);
for i=2:n
    g(i,1) = (y(i,1)-y(i-1,1))/h(i-1,1);
end
for i=2:n-1
    b(i,1) = 6*(g(i+1,1)-g(i,1));
end
sigma=inv(a)*b;
A=zeros(n-1,1);
B=zeros(n-1,1);
C=zeros(n-1,1);
D=zeros(n-1,1);
for i=1:n-1
    A(i,1)=sigma(i+1,1)/(6*h(i,1));
    B(i,1)=sigma(i,1)/(6*h(i,1));
    C(i,1)=(y(i+1,1)/h(i,1))-(sigma(i+1,1)*h(i,1)/6);
    D(i,1)=(y(i,1)/h(i,1))-(sigma(i,1)*h(i,1)/6);
end

i = ones(size(X));
for j=1:n
    i(x(j,1) <= X) = j;
end
Y=A(i,1).*((X-x(i,1)).^3)-B(i,1).*((X-x(i+1,1)).^3)+C(i,1).*(X-x(i,1))-D(i,1).*(X-x(i+1,1));
fileid=fopen('natural_spline.txt','w');
fprintf(fileid,'%s','Interpolated values of y* at given x* ');
fprintf(fileid,'\n');
fprintf(fileid,'%s','Natural Cubic Spline: ');
fprintf(fileid,'\n');
for i=1:size(X)
    fprintf(fileid,'%.4f %.4f\n',X(i,1),Y(i,1));
end
    for i = 1:n-1
    l0 = 1;
    for j = 1:3
        p1 = poly(x(i));
        l0 = conv(p1,l0);
    end
    l1 = 1;
    for j = 1:1
        p = poly(x(i));
        l1 = conv(p,l1);
    end
      l2 = 1;
    for j = 1:3
        p1 = poly(x(i+1));
        l2 = conv(p1,l2);
    end
        l3 = 1;
    for j = 1:1
        p1 = poly(x(i+1));
        l3 = conv(p1,l3);
    end
    term1 = A(i,1)*l0;term2 = B(i,1)*l2;term3 = C(i,1)*l1;term4 = D(i,1)*l3;
    term3(3) = 0;term3(4) = 0;term4(3) = 0;term4(4) = 0;
    p = term1 + term2 + term3 + term4;
for i = 1:n-1
    fprintf(fileid,'%.4f %.4f %.4f %.4f in [%.4f,%.4f]\n',p(4),p(3),p(2),p(1),x(i),x(i+1));
end
for i = 1:n
    d = 3*p(4)*x(i)*x(i) + 2*p(3)*x(i) + p(2);
    fprintf(fileid,'u(%d) = %f ',i,d);
end
fprintf(fileid,'\n');
for i = 1:n
    d = 6*p(4)*x(i) + 2*p(3);
    fprintf(fileid,'v(%d) = %f ',i,d);
end
   
fclose(fileid);
type('output3.txt');
plot(x,y,'ro');
hold on;
for i=1:n-1
    x1=x(i,1):0.001:x(i+1,1);
    y1=A(i,1).*((x1-x(i,1)).^3)-B(i,1).*((x1-x(i+1,1)).^3)+C(i,1).*(x1-x(i,1))-D(i,1).*(x1-x(i+1,1));
    k=plot(x1,y1,'g');
end
xlabel('x');
ylabel('y');
legend(k,'Natural Cubic Spline');
    end
  end

  
  if str2 == 4
      n=length(x);
    h=zeros(n-1,1);
    for i=1:n-1
        h(i,1)=x(i+1,1)-x(i,1);
    end

    a=zeros(n,n);
    a(1,1)=h(2,1);
    a(1,2)=-(h(1,1)+h(2,1));
    a(1,3)=h(1,1);
    a(n,n-2)=h(n-1,1);
    a(n,n-1)=-(h(n-1,1)+h(n-2,1));
    a(n,n)=h(n-2,1);
    for i=2:n-1
        a(i,i-1)=h(i-1,1);
        a(i,i)=2*(h(i-1,1)+h(i-1,1));
        a(i,i+1)=h(i-1,1);
    end
    b=zeros(n,1);
    g=zeros(n,1);
    for i=2:n
        g(i,1) = (y(i,1)-y(i-1,1))/h(i-1,1);
    end
    for i=2:n-1
        b(i,1) = 6*(g(i+1,1)-g(i,1));
    end
    sigma=inv(a)*b;
    A=zeros(n-1,1);
    B=zeros(n-1,1);
    C=zeros(n-1,1);
    D=zeros(n-1,1);
    for i=1:n-1
        A(i,1)=sigma(i+1,1)/(6*h(i,1));
        B(i,1)=sigma(i,1)/(6*h(i,1));
        C(i,1)=(y(i+1,1)/h(i,1))-(sigma(i+1,1)*h(i,1)/6);
        D(i,1)=(y(i,1)/h(i,1))-(sigma(i,1)*h(i,1)/6);
    end
    i = ones(size(X));
    for j=1:n
        i(x(j,1) <= X) = j;
    end
    Y=A(i,1).*((X-x(i,1)).^3)-B(i,1).*((X-x(i+1,1)).^3)+C(i,1).*(X-x(i,1))-D(i,1).*(X-x(i+1,1));
    fileid=fopen('not_a_knot.txt','w');
    fprintf(fileid,'%s','Interpolated values of y* at given x* ');
    fprintf(fileid,'\n');
    fprintf(fileid,'%s','NotaKnot Cubic Spline: ');
    fprintf(fileid,'\n');
    for i=1:size(X)
        fprintf(fileid,'%.4f %.4f\n',X(i,1),Y(i,1));
    end
        for i = 1:n-1
    l0 = 1;
    for j = 1:3
        p1 = poly(x(i));
        l0 = conv(p1,l0);
    end
    l1 = 1;
    for j = 1:1
        p = poly(x(i));
        l1 = conv(p,l1);
    end
      l2 = 1;
    for j = 1:3
        p1 = poly(x(i+1));
        l2 = conv(p1,l2);
    end
        l3 = 1;
    for j = 1:1
        p1 = poly(x(i+1));
        l3 = conv(p1,l3);
    end
    term1 = A(i,1)*l0;term2 = B(i,1)*l2;term3 = C(i,1)*l1;term4 = D(i,1)*l3;
    term3(3) = 0;term3(4) = 0;term4(3) = 0;term4(4) = 0;
    p = term1 + term2 + term3 + term4;
for i = 1:n-1
    fprintf(fileid,'%.4f %.4f %.4f %.4f in [%.4f,%.4f]\n',p(4),p(3),p(2),p(1),x(i),x(i+1));
end
for i = 1:n
    d = 3*p(4)*x(i)*x(i) + 2*p(3)*x(i) + p(2);
    fprintf(fileid,'u(%d) = %f ',i,d);
end
fprintf(fileid,'\n');
for i = 1:n
    d = 6*p(4)*x(i) + 2*p(3);
    fprintf(fileid,'v(%d) = %f ',i,d);
end
    fclose(fileid);
    type('output4.txt');
    plot(x,y,'ro');
    hold on;
    for i=1:n-1
        x1=x(i,1):0.001:x(i+1,1);
        y1=A(i,1).*((x1-x(i,1)).^3)-B(i,1).*((x1-x(i+1,1)).^3)+C(i,1).*(x1-x(i,1))-D(i,1).*(x1-x(i+1,1));
        k=plot(x1,y1,'c');
    end
    xlabel('x');
    ylabel('y');
    legend(k,'Not-a-Knot Cubic Spline');
        end
  end
    
  if str2 == 5
      n=length(x);
h=zeros(n-1,1);
for i=1:n-1
    h(i,1)=x(i+1,1)-x(i,1);
end

a=zeros(n,n);
a(1,1)=2*h(1,1);
a(1,2)=h(1,1);
a(1,n-1)=h(n-1,1);
a(1,n)=2*h(n-1,1);
a(n,1)=1;
a(n,n)=-1;
for i=2:n-1
    a(i,i-1)=h(i-1,1);
    a(i,i)=2*(h(i-1,1)+h(i,1));
    a(i,i+1)=h(i,1);
end
b=zeros(n,1);
g=zeros(n,1);
for i=2:n
    g(i,1) = (y(i,1)-y(i-1,1))/h(i-1,1);
end
for i=2:n-1
    b(i,1) = 6*(g(i+1,1)-g(i,1));
end
b(1)=-6*((y(n,1)-y(n-1,1))/h(n-1,1)) + (6*(y(2,1) - (y(1,1)))/h(1,1));
sigma=inv(a)*b;
A=zeros(n-1,1);
B=zeros(n-1,1);
C=zeros(n-1,1);
D=zeros(n-1,1);
for i=1:n-1
    A(i,1)=sigma(i+1,1)/(6*h(i,1));
    B(i,1)=sigma(i,1)/(6*h(i,1));
    C(i,1)=(y(i+1,1)/h(i,1))-(sigma(i+1,1)*h(i,1)/6);
    D(i,1)=(y(i,1)/h(i,1))-(sigma(i,1)*h(i,1)/6);
end
i = ones(size(X));
for j=1:n
    i(x(j,1) <= X) = j;
end
Y=A(i,1).*((X-x(i,1)).^3)-B(i,1).*((X-x(i+1,1)).^3)+C(i,1).*(X-x(i,1))-D(i,1).*(X-x(i+1,1));
fileid=fopen('periodic_spline.txt','w');
fprintf(fileid,'%s','Interpolated values of y* at given x* ');
fprintf(fileid,'\n');
fprintf(fileid,'%s','Periodic Cubic Spline: ');
fprintf(fileid,'\n');
for i=1:size(X)
    fprintf(fileid,'%.4f %.4f\n',X(i,1),Y(i,1));
end
    for i = 1:n-1
    l0 = 1;
    for j = 1:3
        p1 = poly(x(i));
        l0 = conv(p1,l0);
    end
    l1 = 1;
    for j = 1:1
        p = poly(x(i));
        l1 = conv(p,l1);
    end
      l2 = 1;
    for j = 1:3
        p1 = poly(x(i+1));
        l2 = conv(p1,l2);
    end
        l3 = 1;
    for j = 1:1
        p1 = poly(x(i+1));
        l3 = conv(p1,l3);
    end
    term1 = A(i,1)*l0;term2 = B(i,1)*l2;term3 = C(i,1)*l1;term4 = D(i,1)*l3;
    term3(3) = 0;term3(4) = 0;term4(3) = 0;term4(4) = 0;
    p = term1 + term2 + term3 + term4;
for i = 1:n-1
    fprintf(fileid,'%.4f %.4f %.4f %.4f in [%.4f,%.4f]\n',p(4),p(3),p(2),p(1),x(i),x(i+1));
end
for i = 1:n
    d = 3*p(4)*x(i)*x(i) + 2*p(3)*x(i) + p(2);
    fprintf(fileid,'u(%d) = %f ',i,d);
end
fprintf(fileid,'\n');
for i = 1:n
    d = 6*p(4)*x(i) + 2*p(3);
    fprintf(fileid,'v(%d) = %f ',i,d);
end
fclose(fileid);
type('output5.txt');
plot(x,y,'ro');
hold on;
for i=1:n-1
    x1=x(i,1):0.001:x(i+1,1);
    y1=A(i,1).*((x1-x(i,1)).^3)-B(i,1).*((x1-x(i+1,1)).^3)+C(i,1).*(x1-x(i,1))-D(i,1).*(x1-x(i+1,1));
    k=plot(x1,y1,'m');
end
xlabel('x');
ylabel('y');
legend(k,'Periodic Cubic Spline');
    end
  end
  
  if str2 == 6
      n=length(x);
      s = readmatrix('S_value.txt');
    h=zeros(n-1,1);
    for i=1:n-1
        h(i,1)=x(i+1,1)-x(i,1);
    end
    a=zeros(n,n);
    a(1,1)=2*h(1,1);
    a(1,2)=h(1,1);
    a(n,n)=2*h(n-1,1);
    a(n,n-1)=h(n-1,1);
    for i=2:n-1
        a(i,i-1)=h(i-1,1);
        a(i,i)=2*(h(i-1,1)+h(i,1));
        a(i,i+1)=h(i,1);
    end
    b=zeros(n,1);
    g=zeros(n,1);
    for i=2:n
        g(i,1) = (y(i,1)-y(i-1,1))/h(i-1,1);
    end
    for i=2:n-1
        b(i,1) = 6*(g(i+1,1)-g(i,1));
    end
    b(1,1)=6*(((y(2)-y(1))/h(1,1))-s(1,1));
    b(n,1)=6*(((y(n-1)-y(n))/h(n-1,1)) + s(2,1));
    sigma=inv(a)*b;
    A=zeros(n-1,1);
    B=zeros(n-1,1);
    C=zeros(n-1,1);
    D=zeros(n-1,1);
    for i=1:n-1
        A(i,1)=sigma(i+1,1)/(6*h(i,1));
        B(i,1)=sigma(i,1)/(6*h(i,1));
        C(i,1)=(y(i+1,1)/h(i,1))-(sigma(i+1,1)*h(i,1)/6);
        D(i,1)=(y(i,1)/h(i,1))-(sigma(i,1)*h(i,1)/6);
    end
    i = ones(size(X));
    for j=1:n
        i(x(j,1) <= X) = j;
    end
    Y=A(i,1).*((X-x(i,1)).^3)-B(i,1).*((X-x(i+1,1)).^3)+C(i,1).*(X-x(i,1))-D(i,1).*(X-x(i+1,1));
    fileid=fopen('clamped.txt','w');
    fprintf(fileid,'%s','Interpolated values of y* at given x* ');
    fprintf(fileid,'\n');
    fprintf(fileid,'%s','Clamped Cubic Spline: ');
    fprintf(fileid,'\n');
    for i=1:size(X)
        fprintf(fileid,'%.4f %.4f\n',X(i,1),Y(i,1));
    end
        for i = 1:n-1
    l0 = 1;
    for j = 1:3
        p1 = poly(x(i));
        l0 = conv(p1,l0);
    end
    l1 = 1;
    for j = 1:1
        p = poly(x(i));
        l1 = conv(p,l1);
    end
      l2 = 1;
    for j = 1:3
        p1 = poly(x(i+1));
        l2 = conv(p1,l2);
    end
        l3 = 1;
    for j = 1:1
        p1 = poly(x(i+1));
        l3 = conv(p1,l3);
    end
    term1 = A(i,1)*l0;term2 = B(i,1)*l2;term3 = C(i,1)*l1;term4 = D(i,1)*l3;
    term3(3) = 0;term3(4) = 0;term4(3) = 0;term4(4) = 0;
    p = term1 + term2 + term3 + term4;
for i = 1:n-1
    fprintf(fileid,'%.4f %.4f %.4f %.4f in [%.4f,%.4f]\n',p(4),p(3),p(2),p(1),x(i),x(i+1));
end
for i = 1:n
    d = 3*p(4)*x(i)*x(i) + 2*p(3)*x(i) + p(2);
    fprintf(fileid,'u(%d) = %f ',i,d);
end
fprintf(fileid,'\n');
for i = 1:n
    d = 6*p(4)*x(i) + 2*p(3);
    fprintf(fileid,'v(%d) = %f ',i,d);
end
    %for i = 1:(size)
    fclose(fileid);
    type('output6.txt');
    plot(x,y,'ro');
    hold on;
    for i=1:n-1
        x1=x(i,1):0.001:x(i+1,1);
        y1=A(i,1).*((x1-x(i,1)).^3)-B(i,1).*((x1-x(i+1,1)).^3)+C(i,1).*(x1-x(i,1))-D(i,1).*(x1-x(i+1,1));
        k=plot(x1,y1,'y');
    end
    xlabel('x');
    ylabel('y');
    legend(k,'Clamped Cubic Spline');
    
        end
  end
  end
