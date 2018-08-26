function [f,G]=objfungrad(x)

f=exp(x(1))*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);

t=exp(x(1))*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);

G=[t+exp(x(1))*(8*x(1)+4*x(2)),exp(x(1))*(4*x(1)+4*x(2)+2)];

end 

