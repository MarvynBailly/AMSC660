function lazy(d)

x = [2;0];

if d == 1
    b = -2;
    a = [1 -2];
elseif d == 2
    b = -6;
    a = [-1 -2];
elseif d == 3
    b = -2;
    a = [-1 2];
elseif d == 4
    b = 0;
    a = [1 0];
elseif d == 5
    b = 0;
    a = [0 1];
end

p = [-1  5/2];

if(a*p' < 0)
    alpha = (b-a*x)/(a*p');
    xnew = x + alpha*p'
    fprintf("Accepted with alpha = %d\n",alpha);
else
    fprintf("a'p = %d\n",a*p')
end    
end

