function [timecourse,timereq] = timeforward(x0,a)
rhs = @(t,x) a*x*(1-x);
sol = ode23(rhs,[0 10], x0);
plot(sol.x,sol.y)

timecourse = sol.y;
lg = timecourse>=0.99;
array = sol.x(lg);
timereq = array(1);
end
