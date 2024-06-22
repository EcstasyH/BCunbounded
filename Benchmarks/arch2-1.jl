# lti
# number of variables
@polyvar x1 x2
vars = [x1, x2]

# vector field
f = [ x1^2+x2^2-1,
      5*(x1*x2-1)
     ]
     
gi = [ 2.25 - (x1)^2 - (x2-2)^2]
gu = [-1-x1, x1+3, -2-x2,x2+3]
