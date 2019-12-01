close all
clear
clc

x = linspace(0,2*pi,100);
f = @(x) 3./(5-4*cos(x));
fN2 = uf(x,4);
fN3 = uf(x,6);
fN4 = uf(x,10);

figure(1)
subplot(2,1,1)
hold on
plot(x,f(x),'--','LineWidth',3)
plot(x,fN2,'LineWidth',1)
plot(x,fN3,'LineWidth',1)
plot(x,fN4,'LineWidth',1)
legend("f(x)=3/(5-4cos(x))","N=4","N=6","N=10")
grid on

kint = 2:2:120;
error = 0*kint;
cntr = 1;
for k = kint
    fN = uf(x,k);
    error(cntr) = max(abs(f(x)-fN));
    cntr = cntr + 1;
end
    
figure(1)
subplot(2,1,2)
semilogy(kint,error,'LineWidth',3)
grid on


function res = uf(x,N)
    res = 0;
    for k = -N/2:(N/2-1)
        res = res + 2^(-abs(k))*cos(k*x);
    end
end