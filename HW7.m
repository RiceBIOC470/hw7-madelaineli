%HW7
clear all

%GB comments
1a 90 An initial population of what? Also, by the graphs its clear that some maximum occurs…the more interesting aspect is to explain it. 
1b 100
1c 95 missing axis labels
1d 50 should observe points of bifurcations as a function of a. 
2a 50 wrong equations used. It should be [V/(1+x(2)^4)-x(1); V/(1+x(1)^4)-x(2)];
2b. 100 correctly generated plots, but the equations used are wrong. I will give full credit
2c  100 same as 2b. 
overall: 84


% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).
% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.
% at x = 0: initial population
% at x = 1: point of maximum growth

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 
syms X;
a1 = 1; 
a2 = -1;
figure(1)
fplot(a1*X*(1-X))
figure(2)
fplot(a2*X*(1-X))
% These fixed points does not depends on the value of a. a determines the
% speed at which population reaches the maximum at x = 1

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 
figure(3)
[t1, treq1] = timeforward(0.1,10)

% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

figure(4)
for a = 0:0.5:5
    x0 = rand(1,200);
    for i = 1:200
        x1 = a.*x0(i).*(1-x0(i));
        for j = 1:200
            x2 = a.*x1.*(1-x1);
            x3 = a.*x2.*(1-x2);
            plot(a,x3)
            hold on
        end
    end
end
% the maximum value of x(1-x) = 0.25, as a increases, the maximum value
% increases. and thus the Xf increases. 

% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 

% da/dt = (V+(V-k)*b^4)/(1+b^4)-b
% db/dt = (V+(V-k)*a^4)/(1+a^4)-a

% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 
V = 5;
k = 1;
da = @(t,b) (V+(V-k)*b^4)/(1+b^4)-b;
db = @(t,a) (V+(V-k)*a^4)/(1+a^4)-a;
% A0>B0
A0 = 20;
B0 = 3;
sol_a = ode23(da,[0 10],B0);
sol_b = ode23(db,[0 10],A0);
figure(4)
xlabel('time')
ylabel('expressions')
plot(sol_a.x,sol_a.y,'b')
hold on
plot(sol_b.x,sol_b.y,'g')
% A0<B0
A0 = 3;
B0 = 20;
sol_a = ode23(da,[0 10],B0);
sol_b = ode23(db,[0 10],A0);
figure(5)
xlabel('time')
ylabel('expressions')
plot(sol_a.x,sol_a.y,'b')
hold on
plot(sol_b.x,sol_b.y,'g')

% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 
figure(6)
hold on;
x=1;
for y=0:0.05:5
    polycoeff=[5 5-y -x 0];
    rts=roots(polycoeff);
    rts=rts(imag(rts)==0);
    plot(y*ones(length(rts),1),rts,'r-');
end

