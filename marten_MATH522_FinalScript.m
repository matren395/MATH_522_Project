# Daniel Marten
# 11/22/21
# Error Analysis and Manual Calculation

# Copying values over, could be more optimized
t = [0, 7.143, 14.286, 21.429, 28.571, 35.714, 42.857, 50, 57.143, 64.286, 71.429, 78.571, 85.714, 92.857, 100];
x1 = [0.623, 0.374, 0.249, 0.183, 0.145, 0.12, 0.103, 0.089, 0.078, 0.068, 0.06, 0.053, 0.047, 0.041, 0.037];
x2 = [0, 0.113, 0.151, 0.157, 0.15, 0.137, 0.124, 0.11, 0.098, 0.087, 0.077, 0.068, 0.06, 0.053, 0.047];

# values where the log is linear, figured after 9th timepoint
tLin = t(9:15);
x1Lin = log(x1)(9:15);
x2Lin = log(x2)(9:15);

oneVec = ones(7,1); # just a vector of ones
tLinTransp = transpose(tLin); # dimensions are tricky

# matrix as: [ av11, av12 ; lambda1, lambda1 ]
ansMat = [oneVec , tLinTransp] \ [transpose(x1Lin) , transpose(x2Lin)]

av11 = exp(ansMat(1,1)) # names should be explanatory 
av21 = exp(ansMat(1,2))
lambda11 = ansMat(2,1);
lambda12 = ansMat(2,2);
lambda1 = (lambda11 + lambda12)/2 # lambda1 as the average of the two
lambda1Big = min(lambda11,lambda12);
lambdaErr = -1*(lambda1 - lambda1Big)/(lambda1Big) * 100 

# Finding av2 and lambda2 from the residual

tRes = t(1:6);
x1Res = x1(1:6);
x2Res = x2(1:6);

x1EstErr = (x1 - av11*exp(lambda1*t))(1:6); # Difference between actual and our
x2EstErr = (x2 - av21*exp(lambda1*t))(1:6); # ... predicted values from 
# ... parameters above
# Only working with 1-6

av12 = x1(1) - av11*1 # Finding av2 VERY simple, so that it matches at t=0
av22 = x2(1) - av21*1

x1EElog = log((x1EstErr)/av12); # difference OVER av2
x2EElog = log((x2EstErr)/av22);

# ansMat2 as: two separate values for Lambda 2, averaged later on
ansMat2 = transpose(tRes) \ transpose([x1EElog ; x2EElog])

lambda21 = ansMat2(1);
lambda22 = ansMat2(2);
lambda2 = (lambda21 + lambda22)/2
lambda2Big = min(lambda21,lambda22);
lambdaErr = -1*(lambda2 - lambda2Big)/(lambda2Big) * 100


x1Ans = av11*exp(lambda1*t) + av12*exp(lambda2*t); # our 'answer' for predicted
x2Ans = av21*exp(lambda1*t) + av22*exp(lambda2*t); # ... values

# MSE Computer Code
# 'Check' used, as we are 'checking' accuracy
# x1Ans 

x1Check01 = x1Ans - x1; # difference
x1Check02 = sum(x1Check01.^2); # sum of squared errors
x1mse = sqrt((x1Check02)/15) # root mean squared error
disp("RMSE for better method, parameters, 1.1") # see table S.1


x2Check01 = x2Ans - x2;
x2Check02 = sum(x2Check01.^2);
x2mse = sqrt((x2Check02)/15)
disp("RMSE for better method, parameters, 1.2") # see table S.1

# NOTE: OTHER THAN RMSE ANALYSIS, ALL CODE BELOW THIS POINT IS COPIED FROM 
# ... FROM THE ORIGINAL REPORT
# NO ORIGINAL CODE BELOW HERE

adiv = av11/av21; # v1 as the ratio, av21 is the larger one
bdiv = bv22/bv12; # v2 as the ration

# in AK=b, A is left, b is right 

left = [adiv, 1, 0, 0 ; 1, bdiv, 0, 0, ; 0, 0, adiv, 1 ; 0, 0, 1, bdiv];
right = [adiv*lambda1; lambda2; lambda1 ; bdiv*lambda2];

K = left\right;

k11 = K(1) # assigning individual values from K
k12 = K(2)
k21 = K(3)
k22 = K(4) # should be negative k12 ? 
pDifK12K22 = ((k12-abs(k22))/(k22)) * 100

k01 = lambda1*lambda2/k12
lSum = lambda1+lambda2
kSum = -(k01+k12+k21)
sumDiff = (lSum-kSum)/(lSum) * 100 # ~0.6%, I'll take it!
k01k12 = k01*k12 # identical to below
lProd = lambda1*lambda2 # identical to above!

tSim = 0:1:100; # running ODE45 with our K values to check if they are correct
x10 = 0.623; # x1(0)
x20 = 0.00; # x2(0)

[xsoln,ysoln] = ode45(@odeSim02redo,t,[x10,x20]); # I'm having problems with ODE45
# so I had to hardcode the values in :/ 
# Not best practices

x1Sim = ysoln(:,1);
x2Sim = ysoln(:,2);

x1Check01 = transpose(x1Sim) - x1; # difference
x1Check02 = sum(x1Check01.^2); # sum of squared errors
x1mse = sqrt((x1Check02)/15) # root mean squared error
disp("RMSE for better method, ODE45, 2.1") # see table S.1

x2Check01 = transpose(x2Sim) - x2;
x2Check02 = sum(x2Check01.^2);
x2mse = sqrt((x2Check02)/15)
disp("RMSE for better method, ODE45, 2.2") # see table S.1

kOld = [-0.081889,   0.051496,   0.041126,  -0.050160]
kNew = transpose(K)

kDiff = (kOld - kNew)./(kOld) * 100;
kDiff2 = (kOld-kNew)./(kNew) * 100;






