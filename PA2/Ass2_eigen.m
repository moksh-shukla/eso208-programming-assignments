format long
fID = fopen('eigen3.txt', 'rt');
sizen = 1;
n = fscanf(fID, '%f', sizen);
%disp(n);
    
sizeA = [n n];
A = fscanf(fID, '%f', sizeA);
A = A';

%disp(A);

sizeb = 1;
rel_err = fscanf(fID, '%f', sizeb);
fclose(fID);
%disp(rel_err);

prompt1 = 'Whether largest eigenvalue type (L) or all eigenvalues of A type(A) ';
reply = input(prompt1,'s');

%power method
if strcmp(reply,'L')
  y=zeros(n,1);
  z=zeros(n,1);
  z(1)=1;
  i=0;
  value1=0;
  value2=0;
  for i= 1:100
    y=A*z;
    y=y/sqrt(sum(y.*y));
    z=y;
    if i>1
       value1=(y')*(A*y);
       if abs((value1-value2)/value2) < rel_err/100
         i;
         value2=value1;
         break;
       end
       value2=value1;
       i;
    else
       value2=(y')*(A*y);
    end
  end 
  disp("Maximum Eigen Value is: ");
  disp(value2);
  fprintf("Iterations");
  disp(i);
  fileID = fopen('eigen-power3_out.txt','w');
  fprintf(fileID,'Eigen Vector:\n%.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n %.3f\t\n',y.');
  fprintf(fileID,'Maximum Eigen Value is: %f\nIterations: %d\n', value2,i);
  fclose(fileID);
end

%QR decomposition
if strcmp(reply,'A')
  q=zeros(n);
  r=zeros(n);
  m=zeros(n,1);
  value2=zeros(1,n);
  value1=[];
  i=0;
  j=0;
  k=0;
  for i= 1:100
    for j= 1:n
      m=zeros(n,1);
      value1= [];
      for k= 1:j-1
         m=m+(transpose(A(:,j))*q(:,k))*q(:,k); %%calculation of sigma in gram schimdt formula
      end
      q(:,j)=A(:,j)-m;    %calculation of orthonormal vectors
      q(:,j)=q(:,j)/sqrt(sum(q(:,j).*q(:,j))); %normalizing vector
    end
    for i= 1:n
       for j= i:n
          r(i,j)=transpose(q(:,i))*A(:,j); %calculation of R in QR
       end
       value1=[value1,r(i,i)];
    end
    A=r*q;
    if i>1
       if max(abs((value1-value2)./value2)) < rel_err/100
         break;
       end
    end
    value2=value1;
  end
  disp("Eigen Values are: ");
  disp(value2);
  fileID = fopen('eigen-qr3_out.txt','w');
  fprintf(fileID,'Eigen Values: %f\n', value2);
  fprintf(fileID,'Iterations: %d\n',i);
  fclose(fileID);
     
end