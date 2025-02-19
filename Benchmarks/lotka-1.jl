# number of variables
@polyvar x1 x2 x3
vars = [x1, x2, x3]

# vector field
f = [ x1*(1-x3),
      x2*(1-2*x3),
      x3*(x1+x2-1)
     ];

# description polynomial of the initial/unsafe set, >=0 by default
gi = [ 0.09 - (x1-1)^2 - (x2-1)^2 - (x3-1.5)^2 ]  
gu = [ 0.09 - (x1+1)^2 - (x2+1)^2 - (x3+1.5)^2]
