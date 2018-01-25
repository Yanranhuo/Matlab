% wolfe condition
% yanran Huo

close all;
% Define objective function
f = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Define gradient of objective function
gradf = @(x) [ -400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)) ; 200*(x(2)-x(1)^2) ];

% Define hessian matrix
hessf = @(x) [-400*(x(2)-3*x(1)^2)+2 -400*x(1); -400*x(1) 200];

% Set initial point
xk = [0;1];
fk = f(xk);
k = 0;
Tol = 1e-4;

% Surf and contour plot of the objection function
figure('Position',[30 100 1200 500])
N = 100;
x = linspace(-2,2,N);
[X,Y] = meshgrid(x,x);
fplot = @(x,y) 100*(y-x.^2).^2 + (1-x).^2;
Z = fplot(X,Y);
subplot(1,2,1)
surf(X,Y,Z);
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
hold on;
subplot(1,2,2)
contour(X,Y,Z,60)
hold on
xlabel('x'); ylabel('y');
hold on
% Plot initial point
plot(xk(1),xk(2),'o','MarkerFaceColor','r')

fprintf('iter                   xk                          fk                      alphak\n')
fprintf('-----------------------------------------------------------------------------------------\n')
fprintf('%3d        %3.8e    %3.8e       %3.8e        %3.8e\n',k,xk(1),xk(2),fk, 1)


KeepIterate = true;

while KeepIterate == true
    % Choose pk
    gradfk = gradf(xk);
    hessfk = hessf(xk);
    [~, p]=chol(hessfk);
    if p == 0
        p;
        pk = -(hessfk\gradfk);
        pk = pk/norm(pk);
    else
        p;
        pk = -gradfk/norm(gradfk);
    end

    %set start point
    alpha_0=0;
    alpha_max = 5;
    c1 = 1e-4;
    c2 = 0.9;
    alpha = Wolfecondition(f,gradf,alpha_0,alpha_max,xk,pk,c1,c2);
    xk;
    fk = f(xk);
    grad_last = gradfk;
    xknew = xk+alpha*pk;
    fktest = f(xknew);
    if fktest < fk
        % Plot the xknew
        plot(xknew(1),xknew(2),'o','MarkerFaceColor','r')
        pause(0.6)
        Err_1 = norm(abs(xknew-xk));
        Err_2 = norm(abs(gradf(xknew)-grad_last));
    if Err_1<Tol && Err_2<Tol
        KeepIterate = false;
    end
    k = k+1; xk = xknew; fk = fktest;
    end
    % Print the iteration
plot(xknew(1),xknew(2),'o','MarkerFaceColor','r')
quiver(xk(1),xk(2),pk(1),pk(2),alpha,'k','Linewidth',1);
drawnow;
fprintf('%3d        %3.8e    %3.8e       %3.8e        %3.8e\n',k,xk(1),xk(2),fk,1)
keepIterate = false;
end
disp(xk);