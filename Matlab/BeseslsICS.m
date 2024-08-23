% Script to play with bessel ICS on parabolic problem on disc 




x=linspace(0,10,1001);

plot(x,besselj(0,x))




syms q_n
alpha = 1.0;
eqn = besselj(0, q_n) == 0
guess = 2.5;
z01=vpasolve(eqn, q_n, guess);


alpha=1
beta=0.28
L=5;

r=linspace(0,5,1001)'
u0=besselj(0, z01./L.*r);

figure
plot(r,u0);
hold on 

for t=1:10
    u=u0*exp(-((z01*alpha/L)^2)*t);
    plot(r,u);
end


figure
plot(r,u0);
hold on 

for t=1:10
    u=u0*exp(-((z01*alpha/L)^2-beta)*t);
    plot(r,u);
end