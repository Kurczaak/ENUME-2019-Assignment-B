clear all
close all


step=100; %number of x's to generate the function on

%Data containers preallocation
approximation=zeros(1,step); %values of approximation
x=linspace(-1,1,step); %vactor of linearly spread x's


id=1;%to keep track of the figures' numbers
pair_generator=[0.1,0.5,2];
%---------Task 1 & 2------------
for N=10:10:30
    %Data containers preallocation
    y=zeros(1,N); %vector of y coordinates of generated points
    x_n=zeros(1,N); %vector of x coordinates of generated points
    %Task 1
    
    %Generating different pairs of N and K
    for i=1:3
        K=N*pair_generator(i);
        K=round(K)+1;
        for s=1:step %Generating values of approximation on linearly spread x's
            approximation(s)=Approximate(x(s),N,K);
        end
        figure(id)
        id=id+1;
        hold on;
        plot(x,FirstFunction(x)); %Main function 
        x_n=GenerateXn(N);
        y=GenerateY(N);
        plot(x_n,y,'og') %Generated points
        plot(x,approximation,'r') %Approximation of the main function
        ttle=strcat('Graph made for N=', num2str(N),' And for K= ', num2str(K));
        title(ttle);
        legend('original function','chosen points', 'function Approximate');
        xlabel('-1 < x < 1') 
        ylabel('y') 
        hold off
    end
end


%----------Task3----------
rms=zeros(45,45);
mxe=zeros(45,45);
[N,K]=meshgrid(5:50,5:50);
for n=5:50
    for k=5:50
        if(k<n)
    rms(n-4,k-4)=RMS(n,k); %root-mean-square error
    mxe(n-4,k-4)=MxError(n,k); %maximum error
        else
            rms(n-4,k-4)=NaN;
            mxe(n-4,k-4)=NaN;
        end
    end
end
%RMS
figure(id) 
id=id+1;
surf(K,N,rms);
title('Root-mean-square error dependency on N and K');
xlabel('N') 
ylabel('K') 
zlabel('RMS')

%Maximum error
figure(id)
id=id+1;
surf(K,N,mxe);
title('Maximum error dependency on N and K');
xlabel('N') 
ylabel('K') 
zlabel('Max Error')


%-------Task4----------
N_const=20; %Number of points of the function upon the investigation is done
step=10; %Number of values of sigma upon the investigation is done
STD=logspace(-5,-1,step); %vector containing logaritmically spread values of standard deviation
rms=zeros(N_const,N_const);
rms_min=ones(1,step);
mxe=zeros(N_const,N_const);
mxe_min=ones(1,step);
n_rms=zeros(1,step); %N minimising the RMS
k_rms=zeros(1,step); %K minimising the RMS
n_mxe=zeros(1,step); %N minimising the maximum error
k_mxe=zeros(1,step); %K minimising the maximum error
for i=1:step
    for n=1:N_const
        for k=1:N_const
            if(k<n)
                rms(n,k)=CorruptedRMS(n,k,STD(i));
                mxe(n,k)=CorruptedMxError(n,k,STD(i));
                %RMS investigation
                if(rms(n,k)<rms_min(i))
                    rms_min(i)=rms(n,k);
                    n_rms(i)=n;
                    k_rms(i)=k;
                end
                %Max Error investigation
                if(mxe(n,k)<mxe_min(i))
                    mxe_min(i)=rms(n,k);
                    n_mxe(i)=n;
                    k_mxe(i)=k;
                end
            else
                rms(n,k)=NaN;
                mxe(n,k)=NaN;
            end
        end
    end
end
%RMS
disp('N miminum for RMS:');
n_rms
disp('K miminum for RMS:');
k_rms
figure(id)
id=id+1;
p_rms = polyfit(STD, rms_min,3);
polynomial_rms = polyval(p_rms, STD);
loglog(STD,rms_min, "ro");
hold on
loglog(STD,polynomial_rms,'m');
hold off
title('RMS dependency on the STD');
legend('RMS min','polyfit approximation');
xlabel('STD') 
ylabel('RMS') 


%Maximum error
disp('N miminum for maximum error:');
n_mxe
disp('K miminum for maximum error:');
k_mxe
figure(id)
id=id+1;
p_mxe = polyfit(STD, mxe_min,3);
polynomial_mxe = polyval(p_mxe, STD);
loglog(STD,mxe_min, "ro");
hold on
loglog(STD,polynomial_mxe,'m');
title('Maximum error dependency on the STD');
legend('Maximum error min','polyfit approximation');
hold off
xlabel('STD') 
ylabel('Max error') 


%Function calculating root-mean-square error 
%On parameters N and K
function y=RMS(N,K)
nominator=zeros(1,N);
denominator=zeros(1,N);
x_n=GenerateXn(N);
for i=1:N
nominator(1,i)=Approximate(x_n(i),N,K)-FirstFunction(x_n(i));
denominator(1,i)=FirstFunction(x_n(i));
end
y=norm(nominator,2)/norm(denominator,2);
end

%Function calculating maximum error 
%On parameters N and K
function y=MxError(N,K)
nominator=zeros(1,N);
denominator=zeros(1,N);
x_n=GenerateXn(N);
for i=1:N
nominator(1,i)=Approximate(x_n(i),N,K)-FirstFunction(x_n(i));
denominator(1,i)=FirstFunction(x_n(i));
end
y=norm(nominator,inf)/norm(denominator,inf);
end

%Function generating Matrix of Phi parameters
function y=GeneratePhi(N,K)
Fi=zeros(N,K);
for n=1:N
    x_n=-1+2*(n-1)/(N-1);
    for k=1:K
        Fi(n,k)=Bsk(x_n,k,K);
    end
end
y=Fi;
end

%Generating values of the main function
%For N points
function y=GenerateY(N)
yn=zeros(1,N);
 for n=1:N
   x_n=-1+2*(n-1)/(N-1);
   yn(n)=FirstFunction(x_n);
 end
y=yn;
end

%Generating x values of the N points
function y=GenerateXn(N)
x_n=zeros(1,N);
for n=1:N
   x_n(n)=-1+2*(n-1)/(N-1);
end
y=x_n;
end

%Spline Function, used for Phi matrix generation
function y=Bsk(x,k,K)
xk=x_k(k,K);
y=Bs(2*(x-xk)+2);
end

%Function used for approximation the main function
function y=Approximate(x,N,K)
  Fi=GeneratePhi(N,K);
        y=GenerateY(N);
        p=Fi.'*Fi\Fi.'*y.';
y=0;
K=length(p);
for i=1:K
    y=y+p(i)*Bsk(x,i,K);
end
end

%auxilary function used in B-spline function
function x=x_k(k,K)
x=2*((k-1)/(K-1))-1;
end

%The main function to approximate
function y=FirstFunction(x)
y=-cos(pi.*x).*exp(1).^(-1/3-x);
end

%B-sline function
function y=Bs(x)
if (x>=0 && x<1)
    y=x^3;
elseif (x>=1 && x<2)
    y= -3*(x-1)^3+3*(x-1)^2+3*(x-1)+1;
elseif (x>=2 && x<3)
    y= 3*(x-2)^3-6*(x-2)^2+4;
elseif (x>=3 && x<=4)
    y= -(x-3)^3+3*(x-3)^2-3*(x-3)+1;
else
    y=0;
end
end

%---------Task 4 error corrupted functions-----
function y=CorruptedRMS(N,K,sigma)
nominator=zeros(1,N);
denominator=zeros(1,N);
x_n=GenerateXn(N);
for i=1:N
nominator(1,i)=CorruptedApproximate(x_n(i),N,K,sigma)-FirstFunction(x_n(i));
denominator(1,i)=FirstFunction(x_n(i));
end
y=norm(nominator,2)/norm(denominator,2);
end

function y=CorruptedMxError(N,K,sigma)
nominator=zeros(1,N);
denominator=zeros(1,N);
x_n=GenerateXn(N);
for i=1:N
nominator(1,i)=CorruptedApproximate(x_n(i),N,K,sigma)-FirstFunction(x_n(i));
denominator(1,i)=FirstFunction(x_n(i));
end
y=norm(nominator,inf)/norm(denominator,inf);
end


function y=GenerateCorruptedY(x, sigma)
N=length(x);
yn=zeros(1,N);
 for n=1:N
   x_n=-1+2*(n-1)/(N-1);
   yn(n)=FirstFunction(x_n)+randn()*sigma^2;
 end
y=yn;
end


function y=CorruptedApproximate(x,N,K, sigma)
  Fi=GeneratePhi(N,K);
        x_n=GenerateXn(N);
        y=GenerateCorruptedY(x_n,sigma);
        p=Fi.'*Fi\Fi.'*y.';
y=0;
K=length(p);
for i=1:K
    y=y+p(i)*Bsk(x,i,K);
end
end