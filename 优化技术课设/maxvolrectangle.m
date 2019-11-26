function [low,up,volume]=maxvolrectangle(A,b)
% Calculate the maximum inner rectangle of the polygon,
% where A and b are the semi-space constraints forming convex polygons or
% polyhedrons.

[row,col] = size(A);
combination = nchoosek(1:row,col); % Find the combinations of rows in matrix A
[m,n] = size(combination);
points =[];
for i = 1:m
    A_star = A(combination(i,:),:); % Extract the number of rows corresponding to the combination label
    if rank(A_star)~= n % Check there is a solution or not
        continue; % If no solution,next combination
    else
        x = A_star\b(combination(i,:)); % If there is solution,give the answer
        if all(A*x<=b) % If the point is inside the polygon, it is the boundary point of the polygon
            points =[points x];
        end
    end
end
points =(unique(points','rows'))'; % Remove duplicate boundary points
[row_s,col_s]=size(points);
if row_s >= col_s
    error('The provided constraints cannot form a polygon or polyhedron')
end

% Compute A- and A+
Ap = zeros(row,col);
Am = zeros(row,col);
for i = 1:row
    for j = 1:col
        Ap(i,j)=max(A(i,j),0);
        Am(i,j)=max(-A(i,j),0);
    end
end


% Define general variables
K=10000;
delta = 1e-10;% Inner loop closure condition
epsilon =1e-8;% Outer loop closure condition

% Significant digits is 10
digits(10)
% Define symbol variable
L=sym('L',[row_s 1]);
U=sym('U',[row_s 1]);

% Give the initial value
center = mean(points,2); % Calculate the center of the polygon to determine the position of the initial point
init = 1; % The initial point parameter ensures that the initial point is inside the polygon
X0 = [center-init;center+init];
while ~(all(A*X0(1:row_s)< b) && all(A*X0(row_s+1:2*row_s)< b))
    init = 0.9*init;
    X0 = [center-init;center+init];
end

% Construct penalty function
G = Ap*U-Am*L-b;
phi = sum(-log(-G));

% Function to be optimized
f = sum(-log(U-L));
F = K*f + phi;

% Solve for the gradient and the hessian matrix
F1 = jacobian(F,[L.',U.']);
F2 = jacobian(F1,[L.',U.']);


% Newton method and interior point method iterative optimal value X0	
while(1)
    t=1;
    while(1)
        % The values of the gradient and the hessian matrix are computed numerically
        F1_X = vpa(subs(F1,[L.',U.'],X0'));
        F2_X = vpa(subs(F2,[L.',U.'],X0'));
        % Calculate the next point
        X1 = X0 - t*(F2_X\(F1_X'));
        % Determine whether the next point is a feasible solution, 
        % if not, update step size by exponential attenuation and the attenuation rate of 0.9
        while ~(all(A*X1(1:row_s)< b) && all(A*X1(row_s+1:2*row_s)< b))
            t=t*0.95;
            X1 = X0 - t*(F2_X\(F1_X'));
        end
        % Calculate the function values of the last feasible point and the current feasible point
        f0 = vpa(subs(f,[L.',U.'],X0'));
        f1 = vpa(subs(f,[L.',U.'],X1'));
        % Compare the value of Euclidean norm with last value
        n = norm(f1-f0);
        % If the condition is met to end the current loop, otherwise the current optimal solution is set to a new X0
        if n < delta
            break;
        else
            X0=X1;
        end
    end
    % Determine whether the outer loop meets the end condition, otherwise update the value of K and then calculate
    if row/K < epsilon
        break;
    else
        K=K*10;
    end
end
low = double(X0(1:row_s));
up = double(X0(row_s+1:2*row_s));

% plotting
hold on
axis equal
if row_s == 2
    patch(points(1,:),points(2,:),'c')
    pos =double([X0(1:2)' (X0(3:4)-X0(1:2))']);
    rectangle('Position',pos)
    volume = prod(double(X0(row_s+1:2*row_s)-X0(1:row_s)));
elseif row_s == 3
    view(3) % Open 3D view
    P1 = points';
    K1 = convhulln(P1); % Calculate the convex hull
    patch('Vertice',P1,'Faces',K1,'FaceColor','c','FaceAlpha',.3) % Construct the polytope
    % Build therectangular vertexes
    P2=[ X0(1),X0(1),X0(1),X0(1),X0(4),X0(4),X0(4),X0(4);
        X0(2),X0(2),X0(5),X0(5),X0(2),X0(2),X0(5),X0(5);
        X0(3),X0(6),X0(3),X0(6),X0(3),X0(6),X0(3),X0(6)]';
    P2 = double(P2);
    K2 = convhulln(P2); % Calculate the convex hull
    patch('Vertice',P2,'Faces',K2,'FaceColor','r') % Building cuboid	
    volume = prod(double(X0(row_s+1:2*row_s)-X0(1:row_s)));
else
    disp("High-dimensional data cannot draw images")
end



