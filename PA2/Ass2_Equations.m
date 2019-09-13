prompt = 'What would you like to do?\nA.Solve the system of equation.\nB.Perform an LU Decomposition.\nC.Perform a Matrix Inversion.\n';
t = input(prompt,'s');
%thomas algorithm
if strcmp(t, 'A')
    prompt = 'Is the system is tridiagonal?(Y/N): ';
    reply = input(prompt, 's');
    if strcmp(reply, 'Y')
        fID = fopen('thomas.txt', 'rt');
        sizen = 1;
        n = fscanf(fID, '%f', sizen);
        %disp(n);
        
        sized = n;
        d = fscanf(fID, '%f', sized);
        %disp(d);
        
        sizeb = n;
        b = fscanf(fID, '%f', sizeb);
        %disp(b);
        
        sizeu = n;
        u = fscanf(fID, '%f', sizeu);
        %disp(u);
      
        sizel = n;
        l = fscanf(fID, '%f', sizel);
        %disp(l);
        fclose(fID);
        
        alpha=zeros(1,n);
        beta=zeros(1,n);
        x=zeros(1,n);
        i=0;
        alpha(1)=d(1);
        beta(1)=b(1);
        for i=(2:n)
           alpha(i)=d(i)-l(i)/alpha(i-1)*u(i-1);
           beta(i)=b(i)-l(i)/alpha(i-1)*beta(i-1);
        end
        x(n)=beta(n)/alpha(n);
        for i=(2:n)
           x(n-i+1)=(beta(n-i+1)-u(n-i+1)*x(n-i+2))/alpha(n-i+1);
        end
        disp(x);
        fileID = fopen('thomas_out.txt','w');
        fprintf(fileID,'Roots are:\n');
        fprintf(fileID,'%f\n', x);
        fclose(fileID);
    end
    
    %Gauss-Elimination with partial pivoting
    if strcmp(reply,'N')
     fID = fopen('GaussElimination_partialpivoting.txt', 'rt');
     sizen = 1;
     n = fscanf(fID, '%f', sizen);
        
     sizeA1 = [n n];
     A = fscanf(fID, '%f', sizeA1);
     A = A';
        
     sizeb = n;
     b = fscanf(fID, '%f', sizeb);
     fclose(fID);
     %A1 = reshape(A1,[n,n]);
     %A = [A1 b];
     disp(A);
     % x = zeros;
     i=0; 
     b_=b;
     j=0;
     k=0;
     m=0;
     x=zeros(1,n);
     for k= 1:n-1
       [M,index]=max(abs(A(:,k)));
       A([index k],:) = A([k index],:);
       b([index k]) = b([k index]);
       %A;
       %b;
       for i= k+1:n
         m=A(i,k)/A(k,k);
         A(i,k)=0;
         b(i)=b(i)-m*b(k);
         for j= k+1:n
            A(i,j)=A(i,j)-m*A(k,j);
         end
       end
     end
     x(n)=b(n)/A(n,n);
     for i= n-1:-1:1
       m=0;
       for j= i+1:n 
         m=m+A(i,j)*x(j);
       end   
       x(i)=(b(i)-m)/A(i,i);
       
     end
     disp(x');
     fileID = fopen('Gaussian-Elimination_partial-pivoting_out.txt','w');
     fprintf(fileID,'Roots are:\n');
     fprintf(fileID,'%f\n', x');
     fclose(fileID);
    end
   
end    

if strcmp(t, 'B')
    prompt = 'Is the matrix is symmetric and positive definite?(Y/N): ';
    reply = input(prompt, 's');
    
    %cholesky
    %cholesky
    %cholesky
    if strcmp(reply,'Y')
        fID = fopen('cholesky2.txt', 'rt');
        sizen = 1;
        n = fscanf(fID, '%f', sizen);
        
        sizeA = [n n];
        A = fscanf(fID, '%f', sizeA);
        A = A';
        
        sizeb = n;
        b = fscanf(fID, '%f', sizeb);
        fclose(fID);
        
        
        
         
     exchange_i=[];
     exchange_j=[];
     exchange_with = [];
     
     i=0;
     for i= 1:n-1
       [M,index]=max(abs(A(i:n,i:n)));
       [M,J]=max(M);
       index=index(J);
       A([index+i-1 i],:) = A([i index+i-1],:);
       b([index+i-1 i]) = b([i index+i-1]);
       A(:,[J+i-1 i]) = A(:,[i J+i-1]);
       exchange_i=[exchange_i,index+i-1];
       exchange_j=[exchange_j,J+i-1];
       exchange_with=[exchange_with,i];
     end
    
     exchange_i=[exchange_i;exchange_with];
     exchange_j=[exchange_j;exchange_with];
     [M,N] = size(exchange_i);
     [P,Q] = size(exchange_j);
     
     l=zeros(n,n);
     for j= 1:n
       for i= j:n
         m=0;
         for k= 1:j-1
           m=m+l(i,k)*l(j,k);            
         end
         if i==j
           l(j,j)=sqrt(A(j,j)-m);
         else
           l(i,j)=(A(i,j)-m)/l(j,j);
         end
       end
     end
     disp(l);
     %disp(l*transpose(l));
     disp(exchange_i);
     disp(exchange_j);
     fileID = fopen('cholesky2_out.txt','w');
        fprintf(fileID,"Row Pivoting\n");
        for i = 1:M
            for j = 1:N
                fprintf(fileID,"%d\n",exchange_i(i,j));
            end
            fprintf(fileID,"\n");
        end
        
        fprintf(fileID,"\nColumn Pivoting\n");
        for i = 1:P
            for j = 1:Q
                fprintf(fileID,"%d\n",exchange_j(i,j));
            end
            fprintf(fileID,"\n");
        end
        fprintf(fileID,'\nL Matrix:\n');
        for i = 1:n
            for j = 1:n
                fprintf(fileID,"%f ",l(i,j));
            end
            fprintf(fileID,"\n");
        end
        fprintf(fileID,'\nL Matrix:\n');
        for i = 1:n
            for j = 1:n
                fprintf(fileID,"%f ",l(i,j));
            end
            fprintf(fileID,"\n");
        end
         
        fclose(fileID);
        
    end
    
    
    %lu decomposition
    %lu decomposition
    if strcmp(reply,'N')
        fID = fopen('LUD.txt', 'rt');
        sizen = 1;
        n = fscanf(fID, '%f', sizen);

        sizeA = [n n];
        A = fscanf(fID, '%f', sizeA);
        A = A';

        sizeb = n;
        b = fscanf(fID, '%f', sizeb);
        fclose(fID); 
        disp(A);
 
     exchange_i=[];
     exchange_j=[];
     exchange_with=[];
     i=0;
     for i= 1:n-1
       [M,idx]=max(abs(A(i:n,i:n)));
       [M,J]=max(M);
       idx=idx(J);
       A([idx+i-1 i],:) = A([i idx+i-1],:);
       b([idx+i-1 i]) = b([i idx+i-1]);
       A(:,[J+i-1 i]) = A(:,[i J+i-1]);
       exchange_i=[exchange_i,idx+i-1];
       exchange_j=[exchange_j,J+i-1];
       exchange_with=[exchange_with,i];
     end
     exchange_i=[exchange_i;exchange_with];     
     exchange_j=[exchange_j;exchange_with];
     [M,N] = size(exchange_i);
     [P,Q] = size(exchange_j);
     
     %disp("**");
     %disp(A);
     %disp("**");
     %disp(exchange_i);
     %disp("**");
     %disp(exchange_j);
     
     prompt3 = 'Do you want Crout Type (Y) or Dolittle type (N) ?';
     reply = input(prompt3,'s');
     
     %Crout Algo
     %Crout Algo
     if strcmp(reply,'Y')        
        l=zeros(n);
        u=zeros(n);
        for j= 1:n 
           for i= j:n
              m1=0;
              m2=0;
              for k= 1:i-1 
                m1=m1+l(i,k)*u(k,j);
              end
              for k= 1:i-1 
                m2=m2+l(j,k)*u(k,i);
              end
              l(i,j)=A(i,j)-m1;
              u(j,i)=(A(j,i)-m2)/l(j,j);  
           end
        end
        disp("L");
        disp(l);
        disp("U");
        disp(u);  
        %disp(l*u);
        disp(exchange_i);
        disp(exchange_j);
        fileID = fopen('LUD_crout.txt','w');
        fprintf(fileID,"Row Pivoting\n");
        for i = 1:M
            for j = 1:N
                fprintf(fileID,"%d\n",exchange_i(i,j));
            end
            fprintf(fileID,"\n");
        end
        
        fprintf(fileID,"\nColumn Pivoting\n");
        for i = 1:P
            for j = 1:Q
                fprintf(fileID,"%d\n",exchange_j(i,j));
            end
            fprintf(fileID,"\n");
        end
        	
        fprintf(fileID,'\nL Matrix:\n');
        for i = 1:n
            for j = 1:n
                fprintf(fileID,"%f ",l(i,j));
            end
            fprintf(fileID,"\n");
        end
        fprintf(fileID,'\nU Matrix:\n');
        for i = 1:n
            for j = 1:n
                fprintf(fileID,"%f ",u(i,j));
            end
            fprintf(fileID,"\n");
        end 
        fclose(fileID);
     end
        
     
  %Doolittle Algorithm   
  %Doolittle Algorithm
        if strcmp(reply,'N')
            l=zeros(n);
            u=zeros(n);
            for i= 1:n
               for j= i:n
                  m1=0;
                  for k= 1:i-1  
                    m1=m1+l(i,k)*u(k,j);
                  end
                  m2=0;
                  for k= 1:i-1 
                    m2=m2+l(j,k)*u(k,i);
                  end
                  u(i,j)=A(i,j)-m1;
                  l(j,i)=(A(j,i)-m2)/u(i,i);  
               end
            end
            disp("L");
            disp(l);
            disp("U");
            disp(u);
            %disp(l*u);
            fileID = fopen('LUD_doolittle.txt','w');
        fprintf(fileID,"Row Pivoting\n");
        for i = 1:M
            for j = 1:N
                fprintf(fileID,"%d\n",exchange_i(i,j));
            end
            fprintf(fileID,"\n");
        end
        
        fprintf(fileID,"\nColumn Pivoting\n");
        for i = 1:P
            for j = 1:Q
                fprintf(fileID,"%d\n",exchange_j(i,j));
            end
            fprintf(fileID,"\n");
        end
        fprintf(fileID,'\nL Matrix:\n');
        for i = 1:n
            for j = 1:n
                fprintf(fileID,"%f ",l(i,j));
            end
            fprintf(fileID,"\n");
        end
        fprintf(fileID,'\nU Matrix:\n');
        for i = 1:n
            for j = 1:n
                fprintf(fileID,"%f ",u(i,j));
            end
            fprintf(fileID,"\n");
        end 
        fclose(fileID);
      
        end
    end
end

%inverse of matrix
if strcmp(t,'C') 
    fID = fopen('Matrix_inverse.txt', 'rt');
    sizen = 1;
    n = fscanf(fID, '%f', sizen);
    
    sizeA1 = [n n];
    A1 = fscanf(fID, '%f', sizeA1);
    A1 = A1';
    fclose(fID);
    index = eye(n,n);
   % B = eye(n);
    A = [A1 index];
    %disp(A);
   % I=eye(n,n);
    %B=A;
    A=[A index];
    for k= 1:n
       for j= k+1:2*n
           A(k,j)=A(k,j)/A(k,k);
       end
       A(k,k)=1;
       for i= 1:n
           if i==k
              continue;
           end
           for j= k+1:2*n
              A(i,j)=A(i,j)-A(i,k)*A(k,j);
           end
           A(i,k)=0;
       end 
    end
     A=A(:,n+1:n*2);
     disp("Matrix Inverse:");
     disp(A);
     fileID = fopen('matrix-inverse_out.txt','w');
     fprintf(fileID,"Matrix Inverse is:\n");
     for i = 1:n
         for j = 1:n
             fprintf(fileID,"%f\t",A(i,j));
         end
         fprintf(fileID,"\n");
     end 
     fclose(fileID);
end  