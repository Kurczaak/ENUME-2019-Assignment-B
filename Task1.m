clear
close

%Data containers preallocation
step=100;
App=zeros(1,step);
x=linspace(-1,1,step);
y=zeros(1,10);
x_n=zeros(1,10);
y_t=zeros(1,10);

%to keep track of the figures' numbers
id=1;
for N=10:10:30
    y=zeros(1,N);
    x_n=zeros(1,N);
    y_t=zeros(1,N);
    %Task 1
    index=N/10;
    for i=8:9
        K=N*0.1*i;
        K=round(K);
        for s=1:step
            App(s)=Approximation(x(s),N,K);
        end
        figure(id)
        id=id+1;
        hold on;
        plot(x,FirstFunction(x));
        x_n=GenerateXn(N);
        y=GenerateY(x_n);
        plot(x_n,y,'og')
        plot(x,App,'r')
        ttle=strcat('Graph made for N=', num2str(N),' And for K= ', num2str(K));
        title(ttle);
        legend('original function','chosen points', 'function approximation');
        hold off
    end
end

%Task3
[N,K]=meshgrid(5:50,5:50);
for n=5:50
    for k=5:50
        if(k<n)
    roms(n-4,k-4)=RMS(-0.5,n,k);
    mxe(n-4,k-4)=MxError(-0.5,n,k);
        else
            roms(n-4,k-4)=NaN;
            mxe(n-4,k-4)=NaN;
        end
    end
end
figure(10)
surf(K,N,roms);
figure(11)
surf(K,N,mxe-roms);

function y=RMS(x,N,K)
nominator=norm(Approximation(x,N,K)-FirstFunction(x),2);
denominator=norm(FirstFunction(x),2);
y=nominator/denominator;
end

function y=MxError(x,N,K)
nominator=norm(Approximation(x,N,K)-FirstFunction(x),inf);
denominator=norm(FirstFunction(x),inf);
y=nominator/denominator;
end

function y=GenerateFi(N,K)
Fi=zeros(N,K);
for n=1:N
    x_n=-1+2*(n-1)/(N-1);
    for k=1:K
        Fi(n,k)=Bsk(x_n,k,K);
    end
end
y=Fi;
end

function y=GenerateY(x)
N=length(x);
yn=zeros(1,N);
 for n=1:N
   x_n=-1+2*(n-1)/(N-1);
   yn(n)=FirstFunction(x_n);
 end
y=yn;
end

function y=GenerateXn(N)
x_n=zeros(1,N);
for n=1:N
   x_n(n)=-1+2*(n-1)/(N-1);
end
y=x_n;
end

function y=Bsk(x,k,K)
xk=-1+2*((k-1)/(K-1));
y=Bs(2*(x-xk)+2);
end
function y=Approximation(x,N,K)
  Fi=GenerateFi(N,K);
        x_n=GenerateXn(N);
        y=GenerateY(x_n);
        p=Fi.'*Fi\Fi.'*y.';
y=0;
K=length(p);
for i=1:K
    y=y+p(i)*Bsk(x,i,K);
end
end

function x=x_k(k,K)
x=2*((k-1)/(K-1))-1;
end

function y=FirstFunction(x)
y=-cos(pi.*x).*exp(1).^(-1/3-x);
end

function y=Bs(x)
if (x>=0 && x<1)
    y=x^3;
elseif (x>=1 && x<2)
    y= 3*(x-1)*(x-(x-1)^2)+1;
elseif (x>=2 && x<3)
    y=3*((x-2)^2)*(x-2-2*(x-2))+4;
elseif (x>=3 && x<=4 )
    y=(x-3)*(-3+3*(x-3)-(x-3)^2)+1;
else
    y=0;
end
end


