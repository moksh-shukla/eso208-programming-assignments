promptans = input('Is the input equation a polynomial?(Y/N)\n','s');
if(promptans == 'N')
    promptans1 = input('Choose the method of solution by selecting number shown with method:\n Bisection-1,\n False Position-2,\n Fixed Point Method-3,\n Newton-Raphson-4,\n Secant-5 \n','s');
    
    
    % bisection method
    if(promptans1 == '1')
        clear
        format long 
        func_ask1 = input('Enter your function\n','s');
        str1 = strcat('@(x)', func_ask1);
        f_bisection = str2func(str1);
        a_bisection = input('Enter first starting point\n');
        b_bisection = input('Enter second starting point\n');
        disp('Now subsequent enter stopping criteria:\n');
        rel_err = input('Enter the percentage relative error allowed in solution\n');
        func_value = input('Enter Convergence criteria for the function value, i.e., how close f(x) is to zero\n');
        max_iter = input('Enter allowed maximum number of iterations\n');
                      
        c1 = 0;
        iter = 0;
        c2 = (a_bisection+b_bisection)/2;
        err = abs(100*(c2-c1)/c2);
        stop_err = 0;
        stop_iter=0;
        stop_value=0;
        i=1;

        while ((err > rel_err) && (iter<max_iter) && ((abs(f_bisection(c2))>func_value)))
            if(f_bisection(a_bisection)*f_bisection(b_bisection)<0)
                c2 = (a_bisection+b_bisection)/2;
                if(f_bisection(a_bisection)*f_bisection(c2)<0) b_bisection = c2;
                else a_bisection = c2;
                end
            end
            err = abs(100*(c2-c1)/c2);
            c1 = c2;
            p=c2;
            errg(i)=err;
            iter = iter + 1;
            i=i+1;
            if(err<rel_err)
                stop_err=1;
            end
    
            if(abs(f_bisection(c2))<func_value)
                stop_value=1;
            end
            if(iter==max_iter)
                stop_iter=1;
            end
        end

        fprintf('Root is %f\n',p);

        if(stop_err)
            disp('Iterations stopped as relative error stopping criteria was met');
        end

        if(stop_iter)
            disp('Iterations stopped as maximum number of iterations were executed');
        end

        if(stop_value)
            disp('Iterations stopped as function value was close to zero as required');
        end

        figure(1)
        plot(errg);
        title('Plot of error');
        xlabel('Iterations');
        ylabel('Error');
        figure(2)
     
        ezplot(f_bisection);
        grid on;
    end
    
    % regula-falsi
    if(promptans1 == '2')
        clear
        func_ask2 = input('Enter your function\n','s');
        str1 = strcat('@(x)',func_ask2);
        f_falsepos = str2func(str1);
        a_falsepos = input('Enter first starting point\n');
        b_falsepos = input('Enter second starting point\n');
        disp('Now subsequent enter stopping criteria:');
        rel_err = input('Enter the relative error allowed in solution\n');
        func_value = input('Enter Convergence criteria for the function value, i.e., how close f(x) is to zero\n');
        max_iter = input('Enter allowed maximum number of iterations\n');

        iter = 0;
        x0=a_falsepos;
        x1=b_falsepos;
        err = 100;
        x2 = x1 - f_falsepos(x1)*(x1 - x0)/(f_falsepos(x1) - f_falsepos(x0));
        stop_err = 0;
        stop_iter=0;
        stop_value=0;
        i=1;

        while ((err > rel_err) && (iter<max_iter) && (abs(f_falsepos(x2))>func_value)) 
            x2 = x1 - f_falsepos(x1)*(x1 - x0)/(f_falsepos(x1) - f_falsepos(x0));
            if (f_falsepos(x2)*f_falsepos(x0)<0)
                x1=x2;
                %x1=x1;
            else
                x0=x2;
                %x0=x0;
            end
    
            p=x2;

            err = abs(100*(x2-x1)/x2);
            errg(i) = err;
            iter = iter+1;
            i=i+1;
             if(err<rel_err)
                stop_err=1;
            end
    
            if(abs(f_falsepos(x2))<func_value)
                stop_value=1;
            end
            if(iter==max_iter)
                stop_iter=1;
            end
        end

        fprintf('answer = %f\n',p);

        if(stop_err)
            disp('Iterations stopped as relative error stopping criteria was met');
        end

        if(stop_iter)
            disp('Iterations stopped as maximum number of iterations were executed');
        end

        if(stop_value)
            disp('Iterations stopped as function value was close to zero as required');
        end

        figure(1)
        plot(errg);
        title('Plot of error');
        xlabel('Iterations');
        ylabel('Error');
        figure(2)
        ezplot(f_falsepos);

    end
   
   % fixed point method 
    
    if(promptans1 == '3')
       clear
        format long
        func_ask31 = input('Enter the function\n','s');
        func_ask3 = input('Enter your function g(x) such that your function f(x) is expressed as x=g(x)\n','s');
        str1 = strcat('@(x)',func_ask3);
        str2 = strcat('@(x)',func_ask31);
        f_fixedpoint = str2func(str1);
        f_fixedpoint1=str2func(str2);
        a_fixedpoint = input('Enter starting point\n');
        disp('Now subsequent enter stopping criteria:\n');
        rel_err = input('Enter the relative error allowed in solution\n');
        func_value = input('Enter Convergence criteria for the function value, i.e., how close f(x) is to zero\n');
        max_iter = input('Enter allowed maximum number of iterations\n');

        iter=0;
        err = 100;
        x1 = a_fixedpoint;
        x2 = f_fixedpoint(x1);
        stop_err = 0;
        stop_iter=0;
        stop_value=0;
        i=1;

        while ((err > rel_err) && (iter<max_iter) && (abs(f_fixedpoint(x2))>func_value))
            x2 = f_fixedpoint(x1);
            err = abs(100*(x2-x1)/x2);
            x1 = x2;
            p = x2;
            errg(i) = err;
            iter = iter+1;
            i=i+1;
            if(err<rel_err)
                stop_err=1;
            end

            if(abs(f_fixedpoint(x2))<func_value)
                stop_value=1;
            end
            if(iter==max_iter)
                stop_iter=1;
            end
        end
        
        fprintf('Root is  = %f\n',p);

        if(stop_err)
            disp('Iterations stopped as relative error stopping criteria was met');
        end

        if(stop_iter)
            disp('Iterations stopped as maximum number of iterations were executed');
        end

        if(stop_value)
            disp('Iterations stopped as function value was close to zero as required');
        end

        figure(1)
        plot(errg);
        title('Plot of error');
        xlabel('Iterations');
        ylabel('Error');
        figure(2)
        ezplot(f_fixedpoint1);

    end
    
    % Newton Raphson method
    
    if(promptans1 == '4')
       clear
       format longg
       func_ask4 = input('Enter your function f(x)\n','s');
       str1 = strcat('@(x)',func_ask4);
       func_ask41 = input('Enter first derivative of function f(x)\n','s');
       str2 = strcat('@(x)',func_ask41);
       f_newton1 = str2func(str1);
       f_newton2 = str2func(str2);
       a_newton = input('Enter starting point\n');
       disp('Now subsequent enter stopping criteria:');
       rel_err = input('Enter the relative error allowed in solution\n');
       func_value = input('Enter Convergence criteria for the function value, i.e., how close f(x) is to zero\n');
       max_iter = input('Enter allowed maximum number of iterations\n');

       x1=a_newton;
       x2 = x1 - f_newton1(x1)/f_newton2(x1);
       err = 100;
       iter = 0; 
       stop_err = 0;
       stop_iter=0;
       stop_value=0;
       i=1;

       while ((err > rel_err) && (iter<max_iter) && (abs(f_newton1(x2))>func_value))
           x2 = x1 - f_newton1(x1)/f_newton2(x1);
           err = abs(100*(x2-x1)/x2);
   
    
           if(err<rel_err)
               stop_err=1;
           end

           if(abs(f_newton1(x2))<func_value)
              stop_value=1;
           end
           x1 = x2;
           p=x2;
           errg(i) = err;
           iter = iter+1;
           i=i+1;
           if(iter==max_iter)
               stop_iter=1;
           end
           error(iter)=err;
       end

       fprintf('Root is %f\n',p);

       if(stop_err)
           disp('Iterations stopped as relative error stopping criteria was met');
       end

       if(stop_iter)
           disp('Iterations stopped as maximum number of iterations were executed');
       end

       if(stop_value)
           disp('Iterations stopped as function value was close to zero as required');
       end

       figure(1)
       plot(errg);
       title('Plot of error');
       xlabel('Iterations');
       ylabel('Error');
       figure(2)
       ezplot(f_newton1)


    end
    
    %secant method
    if(promptans1 == '5')
        format longg
        func_ask5 = input('Enter your function\n','s');
        str1 = strcat('@(x)',func_ask5);
        f_secant = str2func(str1);
        a_secant = input('Enter first starting point\n');
        b_secant = input('Enter second starting point\n');
        disp('Now subsequent enter stopping criteria:\n');
        rel_err = input('Enter the relative error allowed in solution\n');
        func_value = input('Enter absolute value of Convergence criteria for the function value, i.e., how close f(x) is to zero\n');
        max_iter = input('Enter allowed maximum number of iterations\n');

        x0 = a_secant;
        x1 = b_secant;  
        x2 = x1 - f_secant(x1)*(x1 - x0)/(f_secant(x1) - f_secant(x0));
        iter=0;
        err = 100;
        stop_err=0;
        stop_val=0;
        stop_iter=0;
        i=1;

        while (((err > rel_err) && (iter<=max_iter) && (abs(f_secant(x2))>func_value)))
            x2 = x1 - f_secant(x1)*(x1 - x0)/(f_secant(x1) - f_secant(x0));
            err = abs(100*(x2-x1)/x2);
            x0 = x1;
            x1 = x2;
            errg(i)=err;
            p=x2;
            i=i+1;
            if(err < rel_err)
                stop_err=1;
            end
            if(iter>max_iter)
                stop_iter=1;
            end
            if(abs(f_secant(x2))<func_value)
                stop_val=1;
            end 
            iter = iter + 1;
        end

        fprintf('Root is %f\n',p);

        if(stop_err)
            disp('Iterations stopped as relative error stopping criteria was met');
        end

        if(stop_iter)
            disp('Iterations stopped as maximum number of iterations were executed');
        end

        if(stop_val)
            disp('Iterations stopped as function value was close to zero as required');
        end

        figure(1)
        plot(errg);
        title('Plot of error');
        xlabel('Iterations');
        ylabel('Error');
        figure(2)
        ezplot(f_secant)

    end 
    
elseif(promptans == 'Y')
    clear
    algo = input('Choose One of the following methods: Muller-1, Bairstow-2\n','s');
    if(algo=='1')
        degree = input('Input Degree of Polynomial\n');
        coef = input('Enter all degree+1 number of coefficients with each input followed by space\n','s');
        coef1 = str2num(coef);
    
        i=0;
        str=num2str(coef1(1));
        while(i<degree)
            i=i+1;
            str = strcat(str,'+(',num2str(coef1(i+1)),')*x^',num2str(i));
        end

        str = strcat('@(x)',str);
        fun = str2func(str);
        i=0;
    
        a_muller = input('Enter first starting point\n');
        b_muller = input('Enter 2nd starting point\n');
        c_muller = input('Enter third starting point\n');
    
        disp('Now subsequent enter stopping criteria:');
        rel_err = input('Enter the relative error allowed in solution\n');
        func_value = input('Enter Convergence criteria for the function value, i.e., how close f(x) is to zero\n');
        max_iter = input('Enter allowed maximum number of iterations\n');

        x1=a_muller;
        x2=b_muller;
        x3=c_muller;
        err=100;
        iter=0;
        stop_err = 0;
        stop_iter=0;
        stop_value=0;
        i=1;

        while((err>rel_err) && (iter<max_iter) && (abs(fun(x3))>func_value))
            c = fun(x3);
            a = ((fun(x3)-fun(x2))/(x3-x2)-(fun(x3)-fun(x1))/(x3-x1))/(x2-x1);
            b = (fun(x3)-fun(x1))/(x3-x1)-a*(x1-x3);
            x1=x2;
            x2=x3;
            if(b>0)
                x3=x3-(2*c/(b+sqrt(b*b-4*a*c)));
            else
                x3=x3-(2*c/(b-sqrt(b*b-4*a*c)));
                p=x3;
            end
            err=abs((x3-x2)/x3);
            errg(i)=err;
            iter=iter+1;
            i=i+1;

            if(err<rel_err)
                stop_err=1;
            end

            if(abs(fun(x3))<func_value)
                stop_value=1;
            end
        
            if(iter==max_iter)
                stop_iter=1;
            end

        end
        answer = x3;
        fprintf('Root is %f+%fi\n',real(p),imag(p));

        if(stop_err)
            disp('Iterations stopped as relative error stopping criteria was met');
        end

        if(stop_iter)
            disp('Iterations stopped as maximum number of iterations were executed');
        end

        if(stop_value)
            disp('Iterations stopped as function value was close to zero as required');
        end
        
        figure(1)
        plot(errg);
        title('Plot of error');
        xlabel('Iterations');
        ylabel('Error');
        figure(2)
        ezplot(fun);
        grid on;
    end 
    
     if(algo=='2')
        func = input('Input your polynomial: ','s');
        f = str2sym(func);
        r = input('Starting value of r: ');
        s = input('Starting value of s: ');
        relApproxErr = input('Allowed Value of relative error: ');
        max_iter = input('Allowed maximum iteration: ');
        y = matlabFunction(f);
        n = polynomialDegree(f) + 1;
        m = n - 1;
        a = sym2poly(f);
        a = fliplr(a);
        M = zeros;
        N = zeros;
        while n > 0
            for j = 1: max_iter
                b = zeros;
                c = zeros;
                b(n) = a(n);
                b(n-1) = a(n-1) + r*b(n);
                c(n) = b(n);
                c(n-1) = b(n-1) + r*c(n);

                for i = n-2: -1: 1 
                    b(i) = a(i) + r*b(i+1) + s*b(i+2);
                    c(i) = b(i) + r*c(i+1) + s*c(i+2);
                end

                syms x y 
                eqn1 = c(2)*x + c(3)*y == -b(1);
                eqn2 = c(3)*x + c(4)*y == -b(2);              
                sol = solve([eqn1, eqn2], [x, y]);
                delr = sol.x;
                dels = sol.y;       
                errPercentr = abs((delr/(r+delr))*100);       
                errPercents = abs((dels/(r+dels))*100);
                N(j) = errPercents;
                M(j) = errPercentr;                
                if errPercents > relApproxErr || errPercentr > relApproxErr 
                    r = r + delr;
                    s = s + dels;
                end 
                if errPercents < relApproxErr && errPercentr < relApproxErr
                    r = r + delr;
                    s = s + dels;
                break;
                end
            end
            m = m - 2;
            if (r^2 + 4*s) < 0
                x1 = r/2;
                x2 = (-(r^2 + 4*s))/2;
                fprintf('Roots of the functions are: %f + %fi, %f - %fi\n', x1, x2, x1, x2);
            else    
                x1 = (r + (r^2 + 4*s)^(1/2))/2;
                x2 = (r - (r^2 + 4*s)^(1/2))/2;    
                fprintf('Roots of the functions are: %f, %f\n', x1, x2);
            end 
            if m > 2
                a = ones;
                for k = n : -1: 3
                   a(k-2) = b(k);
                end  
            end
            n = n - 2;
            if m == 2
                if (b(4)^2 - 4*b(5)*b(3)) < 0
                    x3 = -b(4)/2*b(5);
                    x4 = (-(b(4)^2 - 4*b(5)*b(3)))^(1/2)/2*b(5);
                    fprintf('Roots of the functions are: %f + %fi, %f - %fi\n', x3, x4, x3, x4);
                else
                    x3 = (-b(4) + (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                    x4 = (-b(4) - (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                    fprintf('Roots of the functions are: %f, %f\n', x3, x4);
                end
            break;
            end
            if m == 1
                x1 = -b(3)/b(4);
                fprintf('Roots of the function are: %f\n', x1);
            break;
            end
        end
        figure(1);       
        ezplot(f)
        grid on
        figure(2);
        plot(M)
        title('Error plot of r');
        figure(3);
        plot(N)
        title('Error plot of s');
        
        
    end

            
    
end    
    

             
   
  
    