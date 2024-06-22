# number of variables
@polyvar x1 x2
vars = [x1, x2]
# vector field
f = [ x2,
     -x1 + 0.333*x1^3 - x2]

# description polynomial of the initial/unsafe set, >=0 by default
gi = [ 0.25 - (x1-1.5)^2 - x2^2 ]  
# gu = [ x1+1.5, -0.5-x1, x2+3, -2-x2]
gu = [ 0.16 - (x1+1)^2 - (x2+1)^2 ] #all these methods fail to produce a valid bc
