func = input('Enter the function (dy/dx)\n','s');
func = char(func);
f = inline(func,'x','y');
file = input('Enter name of file to run\n','s');
%file = char(file);
fileID = fopen(file);
inputs = textscan(fileID,'%f');
fclose(fileID);
inputs = cell2mat(inputs);

x0 = inputs(1);
y0 = inputs(2);
xf = inputs(3);
h = inputs(4);
hmax = inputs(5);
alpha = inputs(6);
tolerance = inputs(7);

n = (xf-x0)/h;
x = zeros(n,1);
y = zeros(n,1);
x(1) = x0;
y(1) = y0;

val = menu('Enter the option to solve ODE','Euler Forward','Euler Backward','Trapezoidal','4th order Adams-Bashforth','4th order Adams-Moulton','4th order Backward Difference Formulation','4th order Runge Kutta');
%y=zeros(xf-x0+1);

%Euler Forward
if val==1
    %disp("Hello!")
    for i = 2:n+1
        x(i) = x(i-1)+h;
        y(i) = y(i-1)+h*f(x(i-1),y(i-1));
    end 
        
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','Euler Forward Method');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fID = fopen('output.txt','wt');
    fprintf(fID,'x                 y\n');
    for i = 1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
    
end

%Euler Backward
if val==2
    %Iterative Predictor Corrector
    for i=2:n+1
        x(i)=x(i-1)+h;
        err = 100;
        y2 = y(i-1) + h*f(x(i-1),y(i-1));
        while err>tolerance
            t = y2;
            y2 = y(i-1)+h*(f(x(i),t));
            err = abs((y2-t)/y2*100);
        end
        y(i) = y2;
    end
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName',' Euler Backward');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fid = fopen('output.txt','wt');
    fprintf(fid,'x      y\n');
    for i=1:n+1
        fprintf(fid,'%f  %f\n',x(i),y(i));
    end
end

%Trapezoidal
if val==3
    for i=2:n+1
        x(i)=x(i-1)+h;
        err = 100;
        y2 = y(i-1) + h*f(x(i-1),y(i-1));
        while err>tolerance
            t = y2;
            y2 = y(i-1)+h*0.5*(f(x(i),t)+f(x(i-1),y(i-1)));
            err = abs((y2-t)/y2*100);
        end
        y(i) = y2;
    end
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','Trapezoidal Method');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fid = fopen('output.txt','wt');
    fprintf(fid,'x           y\n');
    for i=1:n+1
        fprintf(fid,'%f  %f\n',x(i),y(i));
    end
end

%Adam's Bashforth
if val == 4
    
    x(2) = x(1)+h;
    y(2) = y(1)+h*f(x(1),y(1));
    x(3) = x(2)+h;
    y(3) = y(2)+h*(3*f(x(2),y(2)) - f(x(1),y(1)))/2;
    x(4) = x(3)+h;
    y(4) = y(3)+h*((23/12)*f(x(3),y(3))-(4/3)*f(x(2),y(2))+(5/12)*f(x(1),y(1)));
    
    for i = 5:n+1
        x(i) = x(i-1)+h;
        y(i) = y(i-1)+h*((55/24)*f(x(i-1),y(i-1))-(59/24)*f(x(i-2),y(i-2))+(37/24)*f(x(i-2),y(i-2))-(3/8)*f(x(i-3),y(i-3)));
    end
    
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','4th order Adams Bashforth');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fID = fopen('output.txt','wt');
    fprintf(fID,'x                 y\n');
    for i = 1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
    
end

%4th order Adams Moulton
if val == 5
    for i=2:4
        x(i)=x(i-1)+h;
        y(i)=y(i-1)+h*f(x(i-1),y(i-1));
    end
    for i=5:n+1
        x(i)=x(i-1)+h;
        err = 100;
        y2=y(i-1)+h*((55/24)*f(x(i-1),y(i-1))-(59/24)*f(x(i-2),y(i-2))+(37/24)*f(x(i-3),y(i-3))-(3/8)*f(x(i-4),y(i-4)));
        while err>tolerance
            t = y2;
            y2=y(i-1)+h*((9/24)*f(x(i),t)+(19/24)*f(x(i-1),y(i-1))-(5/24)*f(x(i-2),y(i-2))+(1/24)*f(x(i-3),y(i-3)));
            err = abs((y2-t)/y2*100);
        end
        y(i) = y2;
    end
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','4th-order Adams-Moulton');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fid = fopen('output.txt','wt');
    fprintf(fid,'x     y\n');
    for i=1:n+1
        fprintf(fid,'%f  %f\n',x(i),y(i));
    end
end   

%4th order BDF
if val==6
     for i=2:4
        x(i)=x(i-1)+h;
        y(i)=y(i-1)+h*f(x(i-1),y(i-1));
    end
    for i=5:n+1
        x(i)=x(i-1)+h;
        err = 100;
        y2 = y(i-1) + h*f(x(i-1),y(i-1));
        while err>tolerance
            t = y2;
            y2 = (12/25)*(h*f(x(i),t)+4*y(i-1)-3*y(i-2)+(4/3)*y(i-3)-.25*y(i-4));
            err = abs((y2-t)/y2*100);
        end
        y(i) = y2;
    end
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','4th-order(BDF)');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fid = fopen('output.txt','wt');
    fprintf(fid,'x      y\n');
    for i=1:n+1
        fprintf(fid,'%f  %f\n',x(i),y(i));
    end    
end
%4th order Runge Kutta
if val==7
    for i = 2:n+1
        x(i) = x(i-1)+h;
        phi0 = f(x(i-1),y(i-1));
        phi1 = f(x(i-1)+0.5*h,y(i-1)+0.5*h*phi0);
        phi2 = f(x(i-1)+0.5*h,y(i-1)+0.5*h*phi1);
        phi3 = f(x(i-1)+h,y(i-1)+h*phi2);
        y(i) = y(i-1)+h*(phi0/6+(phi1+phi2)/3+phi3/6);
    end
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','4th order Runge Kutta');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fID = fopen('output.txt','wt');
    fprintf(fID,'x                 y\n');
    for i = 1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
end