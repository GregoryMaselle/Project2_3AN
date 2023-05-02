%Coupled schrodinger equations
%Gregory Maselle MSLGRE001
%26/04/2023
%Solving coupled schrodinger equations using explicit finite difference
%method. periodic boundary conditions.

%Firstly lets choose an interval to work over. Lets start with [-10,10]
N = 100;% How many intervals we are having
a = -10; % Right side of interval
b = 10; % Left side of interval
h = (b-a)/N; % Step Size
j = 0:N; % Array of n integers
x = a + j*h; % The array of x points, xj.
c = 1; % This si the velocity c for now.
x1 = 0;
x2 = 0; % these are the values used in the initial condition for u. 
alpha =1;

epsilon = 1e-6;
%Lets solve for u first ONLY with a dummy variable for v.

%%Initial conditions:
u = (sech(x-x1)).*exp(1i*c*(x-x1)) + (sech(x-x2)).*exp(1i*c*(x-x2)); % this is for t = 0.
u = u'; % we now have the initial conditions set for u. It is the first column of our solution matrix since it is t = 0.
v = (sech(x-1)).*exp(1i*c*(x-1)) + (sech(x-2)).*exp(1i*c*(x-2));
v= v';

%% Let us set the values for the times.
t = 0.001; % this is the size of the time steps.
M = 100; % this is the amount of time iterations it will do.
k = 0:M;
T = t*k;
n =1;
%% Lets set up the transformation matrices.
s = t/(2*h^2);
main1 = (1i -2*s)*ones(1,N);
side = s*ones(1,N-1);
A = diag(main1,0) + diag(side,1) + diag(side,-1);
A(1,end) = s;
A(end,1) = s;

main2 = (1i +2*s)*ones(1,N);
side2 = -1*side;

B = diag(main2,0) + diag(side2,1) + diag(side2,-1);
B(1,end) = -s;
A(end,1) = -s;

%% Let us loop through.
while (n <= M)
    NLTermu =zeros(N,1);
    NLTermv = zeros(N,1);
    for c= 1:N
        if c == 1 
            NLTermu(c) = (abs(u(end-1,n))^2 + (abs(v(end-1,n)))^2/(1+ alpha*(abs(v(end-1,n)))^2))*u(end-1,n) ...
                        + (abs(u(2,n))^2 + (abs(v(2,n)))^2/(1+ alpha*(abs(v(2,n)))^2))*u(2,n);

            NLTermv(c) = (abs(v(end-1,n))^2 + (abs(u(end-1,n)))^2/(1+ alpha*(abs(u(end-1,n)))^2))*v(end-1,n) ...
                        + (abs(v(2,n))^2 + (abs(u(2,n)))^2/(1+ alpha*(abs(u(2,n)))^2))*v(2,n);
        elseif c == N
            NLTermu(c) = (abs(u(end-2,n))^2 + (abs(v(end-2,n)))^2/(1+ alpha*(abs(v(end-2,n)))^2))*u(end-2,n) ...
                        + (abs(u(1,n))^2 + (abs(v(1,n)))^2/(1+ alpha*(abs(v(1,n)))^2))*u(1,n);

            NLTermv(c) = (abs(v(end-2,n))^2 + (abs(u(end-2,n)))^2/(1+ alpha*(abs(u(end-2,n)))^2))*v(end-2,n) ...
                        + (abs(v(1,n))^2 + (abs(u(1,n)))^2/(1+ alpha*(abs(u(1,n)))^2))*v(1,n);
        else
            NLTermu(c) = (abs(u(c-1,n))^2 + (abs(v(c-1,n)))^2/(1+ alpha*(abs(v(c-1,n)))^2))*u(c-1,n) ...
                        + (abs(u(c+1,n))^2 + (abs(v(c+1,n)))^2/(1+ alpha*(abs(v(c+1,n)))^2))*u(c+1,n);

            NLTermv(c) = (abs(v(c-1,n))^2 + (abs(u(c-1,n)))^2/(1+ alpha*(abs(u(c-1,n)))^2))*v(c-1,n) ...
                        + (abs(v(c+1,n))^2 + (abs(u(c+1,n)))^2/(1+ alpha*(abs(u(c+1,n)))^2))*v(c+1,n);
        end
    end
    u(1:end-1,n+1) = A\( B*u(1:end-1,n) - t*NLTermu );
    u(end,n+1) = u(1,n+1);
    v(1:end-1,n+1) = A\( B*v(1:end-1,n) - t*NLTermv );
    v(end,n+1) = v(1,n+1);

    n= n+1;
end

figure(1), surf(x,T,real(u))
figure(2), surf(x,T,real(v))
figure(3), surf(x,T,imag(u))
figure(4), surf(x,T,imag(v))
% for h = 1:M
%     plot(real(u(:,h)))
%     pause
% end


            
    
