k =  toField (QQ[i]/(i^2+1))
R = k[x,y]


a = 1
b = 1
XX = matrix {{a,1},{0,a^(-1)}}
xx = inverse XX
YY = matrix {{b,0},{x+yi,b^(-1)}}
yy = inverse YY

lookupTable = {{xx, XX}, {YY, yy}};

modifiedCeiling = n -> ceiling(if ceiling n == n then n + 1/2 else n);

word = (p,q) -> product(apply(1 .. 2*q, i-> lookupTable#(if i%2 == 0 then 0 else 1)#(if modifiedCeiling(i*p/q)%2 == 0 then 0 else 1)))
farey = (p,q) -> trace word (p,q)

farey(0,1)
substitute(farey(0,1),{i=>0})
