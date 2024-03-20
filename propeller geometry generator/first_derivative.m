function d = first_derivative(x,y)
    n=length(x);    
    d=zeros(n,1);
    for i=2:n-1
        h1=x(i)-x(i-1);
        h2=x(i+1)-x(i);
        f0=y(i-1);
        f1=y(i);
        f2=y(i+1);
        d(i)=-h2*f0/h1/(h1+h2) - (h1-h2)*f1/h1/h2 + h1*f2/h2/(h1+h2);
    end
    
    h1=x(2)-x(1);
    h2=x(3)-x(2);
    f0=y(1);
    f1=y(2);
    f2=y(3);
    d(1)=-f0*(2*h1+h2)/h1/(h1+h2) + (h1+h2)*f1/h1/h2 -h1*f2/h2/(h1+h2);
    
    h1=x(n-1)-x(n-2);
    h2=x(n)-x(n-1);
    f0=y(n-2);
    f1=y(n-1);
    f2=y(n);
    d(n)=h2*f0/h1/(h1+h2) - f1*(h1+h2)/h1/h2 + f2*(h1+2*h2)/h2/(h1+h2);
end

