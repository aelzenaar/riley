k = QQ[i]/(i^2+1);
R = k[x,y];
S = QQ[x,y];
use R;

a = 1;
b = 1;
XX = matrix (R, {{a,1},{0,a^(-1)}});
xx = matrix (R, {{a^(-1),-1},{0,a}});
YY = matrix (R, {{b,0},{x+y*i,b^(-1)}});
yy = matrix (R, {{b^(-1),0},{-x-y*i,b}});

lookupTable = {{xx, XX}, {YY, yy}};

modifiedCeiling = n -> ceiling(if ceiling n == n then n + 1/2 else n);

word = (p,q) -> product toList apply(1 .. 2*q, n-> lookupTable#(if n%2 == 0 then 0 else 1)#(if modifiedCeiling(n*p/q)%2 == 0 then 0 else 1));
farey = (p,q) -> trace word (p,q);

pleatingComponents = (p,q) -> factor substitute(-i*farey(p,q),S)

coprime = n -> select(1..n, m->(gcd(m,n)==1))

