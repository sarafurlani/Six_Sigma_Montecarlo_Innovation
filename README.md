https://pastebin.com/raw/K5rqVMCE

# Montecarlo Simulation Functions

# Mode Calculation
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Lognormal - Return MeanLog and SDLog for a given Mean and SD
GetLogMeanSD <- function(NormMean, NormSD) {
  LogMean = log(NormMean^2 / sqrt((NormMean^2)+(NormSD^2)))
  LogSD = sqrt(log(1+(NormSD^2)/(NormMean^2)))
  return(c(LogMean, LogSD))
}

# Triangular Distribution - Quantile Function
qtri <- function(p, a, b, c) {
  if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
  
  ifelse(p > 0 & p < 1,
         ifelse(p <= ptri(c, a, b, c),
                a+sqrt((a^2-a*b-(a-b)*c)*p),
                b-sqrt(b^2-a*b+(a-b)*c+(a*b-b^2-(a-b)*c)*p)),
         NA)
}

# Triangular Distribution - Random Generation
rtri <- function(n, a,b,c) {
  if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
  qtri(runif(n, min = 0, max = 1), a, b, c)
}

# Triangular Distribution - Density Function
dtri <- function(x, a, b, c) {
  if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
  
  ifelse(x >= a & x <= b,
         ifelse(x <= c,
                2*(x-a)/((b-a)*(c-a)),
                2*(b-x)/((b-a)*(b-c))),
         0)
}

# Triangular Distribution - Probability Function
ptri <- function(x, a, b, c) {
  if (!all(a <= c && c <= b && a < b)) stop("It must be that a ≤ c ≤ b and a < b!");
  
  ifelse(x > a & x < b,
         ifelse(x <= c,
                ((x-a)^2)/((b-a)*(c-a)),
                1-((b-x)^2)/((b-a)*(b-c))),
         ifelse(x <= a,
                0,
                1))
}

# Normalise a vector
Normalise <- function(x){
  range = max(x)-min(x)
  if(range == 0) range = 1
  return((x-mean(x))/range)
}

# MonteCarlo Sensitivity from Function
# And X Values; X must be a List
MonteCarloSens <- function(FUN, X) {
  n = length(X)
  xnames = names(X)
  xvals = unlist(X, use.names = F)
  m = matrix(rep(1,n^2),nrow = n, ncol = n)
  m = m + diag(0.01,n)
  y = do.call(FUN,as.list(xvals))
  sens = c()
  for(i in 1:n) {
    sens[i] = 100*(do.call(FUN,as.list(t(m *xvals)[i,])) - y)/y
  }
  xcol = rep("",n)
  xcol[which(sens >= 0)] = "darkolivegreen3"
  xcol[which(sens < 0)] = "brown1"
  index = order(abs(sens), decreasing = F)
  graphics.off()
  par(las=2); m = par("mar"); par(mar = c(m[1],2*m[2],m[3],m[4]))
  barplot(abs(sens[index]), names.arg = xnames[index], 
          horiz = TRUE, col = xcol[index])
  cat("\n")
  return(data.frame(Variable = xnames[rev(index)],
                    PercentSens = sens[rev(index)]))
}


# PerformMonteCarlo Simulation
# Output Simulated X-Y Data &Histogram
SimulateX <- function(nsim, Function, df) {
  
  if(length(formals(Function)) != nrow(df)) {
    cat("\n","The no of rows in the data frame is not equal", "\n",
        "to the number of arguments in the function.","\n", "\n")
  }
  # Create Data Frame of Simulated Xs
  nvar = length(formals(Function))
  X = data.frame(matrix(NA, ncol=nvar, nrow=nsim))
  colnames(X) = df[,1]
  for(i in 1:nvar) {
    parm = list(df[i,3],df[i,4],df[i,5])
    parm = parm[!is.na(parm)]
    X[,i] = do.call(as.character(df[i,2]),c(nsim,parm))
  }
  return(X)
}

# PerformMonteCarlo Simulation
# Output Simulated X-Y Data &Histogram
SimulateY <- function(Function, X) {
  
  if(length(formals(Function)) != ncol(X)) {
    cat("\n","The no of columns in X dataframe is not equal", "\n",
        "to the number of arguments in the function.","\n", "\n")
  }
  
  # Calculate Simulated Y from Function
  out = c()
  Xlist <- as.list(as.data.frame(t(X)))
  for(i in 1:nrow(X)) {
    out[i] = do.call(GetProfit,as.list(Xlist[[i]]))
  }
  hist(out, main = "Histogram of Simulated Output", xlab = NULL)
  return(out)
}

# MonteCarlo Sensitivity for all Simulated Xs
# Input - nsim, Transfer Function, df
# Produces Boxplot of Sensitivities
MonteCarloSens2 <- function(Function, X) {
  
  if(length(formals(Function)) != ncol(X)) {
    cat("\n","The no of columns in X dataframe is not equal", "\n",
        "to the number of arguments in the function.","\n", "\n")
  }
  
  nvar = ncol(X)
  nsim = nrow(X)
  
  # Calculate Sensitivity Matrix
  S = data.frame(matrix(0, ncol=nvar, nrow=nsim))
  colnames(S) = names(X)
  for(i in 1:nsim) {
    xvals = unlist(X[i,], use.names = F)
    m = matrix(rep(1,nvar^2),nrow = nvar, ncol = nvar)
    m = m + diag(0.01,nvar)
    y = do.call(Function,as.list(xvals))
    for(j in 1:nvar) {
      S[i,j] = 100*(do.call(Function,as.list(t(m *xvals)[j,])) - y)/y
    }
  }
  
  med = as.vector(mapply(median,S))
  xcol = rep("",nvar)
  xcol[which(med >= 0)] = "darkolivegreen3"
  xcol[which(med < 0)] = "brown1"
  index = order(abs(med), decreasing = F)
  
  graphics.off()
  par(las=2); mg = par("mar"); par(mar = c(mg[1],2*mg[2],mg[3],mg[4]))
  boxplot(S[index],horizontal=TRUE,outline = FALSE, col = xcol[index])
  
  q1 = as.vector(mapply(quantile,S)[2,])
  q3 = as.vector(mapply(quantile,S)[4,])
  
  w1 = q1 - 1.5*(q3-q1)
  w2 = q3 + 1.5*(q3-q1)
  
  Y = data.frame(Variable = colnames(S)[rev(index)],
                 WLeft = w1[rev(index)],
                 Q1 = q1[rev(index)],
                 Median = med[rev(index)],
                 Q3 = q3[rev(index)],
                 WRight = w2[rev(index)])
  
  return(Y)
}

# Throughput Equation
# Wo = Critical WIP
# To = Raw Process Time
# Rb = Bottleneck Rate
# TH = Throughput
# FT = Flow Time

Throughput <- function(Rb, To, WIP) {
  Wo = Rb * To
  TH = Rb * WIP / (Wo + WIP - 1)
  return(TH)
}

# Mortgage Calculator
MonthlyPayment <- function(Amount, Months, Rate) {
  Rate = Rate/1200
  m = (Amount * Rate * (1+Rate)^Months)/((1+Rate)^Months - 1)
  return(round(m,2))
}


cat("\014")

#from rstudio exercise


rbinom(n, size, prob)
# prob = 0.9; probability of good products
# size = trials (daily quantity produced = 100)
# n = number of observations (days in a year = 365)
# x = number of good products per day in 365 days

x = rbinom(365, 100, 0.9)       
hist(x); print(x[1:10])

how many good products are produced every day

# Generate Integers - Poisson Distribution
# Number of times an event can occur 
# Number of arrivals/hour at petrol station
# Not possible to count non-events

x = rpois(100,1)             # rpois(n, lambda)
hist(x); print(x[1:10])      # lambda = expected rate


# Generate Integers - User Defined Distribution

x = sample(c(1, 2, 3), size = 100, replace = TRUE, 
           prob = c(0.5, 0.1, 0.4))
print(x); hist(x)

# Generate Discrete Categories - Uniform Distribution
x = sample(c("Eastern","Central","Mountain","Pacific"), 
           size = 50, replace = TRUE)
print(x); table(x); barplot(table(x))

# Generate Discrete Categories - User Defined Distribution
x = sample(c("pass","fail"), size = 20, replace = TRUE, 
           prob = c(0.85, 0.15))
print(x); table(x); barplot(table(x))

# Generate Discrete Categories - User Defined Distribution
x = sample(c("pass","fail"), size = 20, replace = TRUE, 
           prob = c(0.85, 0.15))
print(x); table(x); barplot(table(x))

# Exercise - Generate These Categories
#            A, E, I, O, U - 40% prob, uniform dist
#            X, Y, Z       - 60% prob, uniform dist
#            Total Number  - 1000

x1=c("A", "E", "I", "O", "U")
x2=c("X","Y","Z")

p1=rep(0.4/5, 5)
p2=rep(0.6/3, 3)
print(p1); print(p2)

#Task 2 - Random Generation of Continuous Variables
#          Uniform, Normal, Exponential Distribution
#          Gamma, Triangular, Log-Normal

# Generate Continuous Data - Uniform Distribution
x = runif(1000, 10,20)     # n, min, max
hist(x); print(x[1:6])
# Generate Continuous Data - Normal Distribution
x = rnorm(1000, mean = 10, sd = 2)
hist(x); print(x[1:6])

# Generate Continuous Data - Exponential Distribution
x = rexp(1000, rate = 0.25)  # mean = 1/rate
hist(x); print(x[1:6])

# Generate Continuous Data - Gamma Distribution

x = rgamma(1000, shape = 1, rate = 0.25)   # mean = shape/rate
hist(x); print(x[1:6]); mean(x)           

x = rgamma(1000, shape = 2, rate = 0.5)    # mean = shape/rate
hist(x); print(x[1:6]); mean(x); sd(x)     # sd = shape/rate^2    

x = rgamma(1000, shape = 4, rate = 1)      # mean = shape/rate
hist(x); print(x[1:6]); mean(x); sd(x)     # sd = shape/rate^2 

x = rgamma(1000, shape = 16, rate = 4)     # mean = shape/rate
hist(x); print(x[1:6]); mean(x); sd(x)     # sd = shape/rate^2     

x = rgamma(1000, shape = 64, rate = 16)    # mean = shape/rate
hist(x); print(x[1:6]); mean(x); sd(x)    # sd = shape/rate^2 

x = rgamma(1000, shape = 640000, rate = 160000)     # mean = shape/rate
hist(x); print(x[1:6]); mean(x); sd(x)

# Triangular Distribution
x = rtri(1000, 5, 15, 10)               # n, min, max, mid
hist(x); print(x[1:6]); mean(x)         

# Generate Continuous Data - Lognormal Distribution
MeanLog = GetLogMeanSD(50,20)[1]  # provide mean and std dev
SDLog   = GetLogMeanSD(50,20)[2]  # provide mean and std dev

x = rlnorm(1000, MeanLog, SDLog)
hist(x)

# Exercise - Generate 1000 Continuous Values
#      Picked Randomly, 30% from A, 70% from B
#      A - Normal Dist, Mean = 100, Stdev = 12
#      B - Exponential Dist, Rate = 0.01

x = rnorm(300, 100, 12)
y = rexp(700, 0.01)

z = sample(c(x,y), size = 1000)
hist(z)


#Task 3 - Transfer Function; Y = function(X1, X2, etc)
  #          Triange made of 3 bars, AB, BC, CA 
  #          See https://pasteboard.co/JMPvQPaY.png
  #          Each bar is 100 units in length
  #          Angle C is 60 degrees as designed
  #          Mfg variation in bars +/- 1.0 
  #          Calculate Angle Distribution
  
AB = 99.8; BC = 100.12; CA = 99.94

AngleC = acos( (BC^2 + CA^2 - AB^2) / (2*BC*CA) )  # cosine law
AngleC = AngleC * 180/pi                           # radians to degrees

cat("\n", "AngleC :", AngleC, "\n")
    
# Simulating Distribution of Angle-C
    
AB = rnorm(1000, 100, 1/3)   # 99.73 % fall within 3 std devs
BC = rnorm(1000, 100, 1/3)
CA = rnorm(1000, 100, 1/3)

AngleC = acos( (BC^2 + CA^2 - AB^2) / (2*BC*CA) )   # cosine law
AngleC = AngleC * 180/pi                            # CONVERSION:radians to degrees

AngleC[1:10]; mean(AngleC); sd(AngleC)
hist(AngleC)

# Task 4 - Transfer Function; Y = function(X1, X2, etc)
#          Florist Business 
#          Small, Medium, Large (Profit: $5, $7, $10)
#          Avg. Orders (20,15,10)
#          Cake Orders (Profit: $10) - 5% orders
#          Chocolates  (Profit: $3)  - 10% orders
#          If > 50 bouquets are in a day
#               handling costs 10% of the profit
#          Calculate Annual Profit Distribution
#-----

VS = 20; VM = 15; VL = 10     
PS = 5;  PM = 7;  PL = 10  
XCake = 0.05; XChoc = 0.10
PCake =  10;  PChoc = 3; Max = 50

VCake = XCake * (VS+VM+VL) 
VChoc = XChoc * (VS+VM+VL)   

# Step 1 - Using R Function As A Transfer Function 

GetProfit <- function(VS,VM,VL, PS,PM,PL, VCake,VChoc, PCake,PChoc, Max) {
  Profit = VS*PS + VM*PM + VL*PL + 
    VCake*PCake + VChoc*PChoc
  if(VS+VM+VL > Max) Profit = 0.9 * Profit
  return(Profit)
}

GetProfit(VS,VM,VL, PS,PM,PL, VCake,VChoc, PCake,PChoc, Max)

GetProfit(VS+5,VM,VL, PS,PM,PL, VCake,VChoc, PCake,PChoc+3, Max)

# Step 2 - Generate Simulated Xs
#          From Statistical Distributions
#          Vectors prefixed with "x"

nsim = 365                # Days in the Year

xVS = rpois(nsim, VS)     # Daily Qty of Small Boquets
xVM = rpois(nsim, VM)     # Daily Qty of Medium Boquets
xVL = rpois(nsim, VL)     # Daily Qty of Large Boquets

xVCake = rbinom(nsim, VS+VM+VL, XCake) # Daily Qty of Cake Orders
xVChoc = rbinom(nsim, VS+VM+VL, XChoc) # Daily Qty of Choc Orders

#put all the outcomes in a table
Xtable = data.frame(Qty_Small    = xVS, 
                    Qty_Medium   = xVM,
                    Qty_Large    = xVL,
                    Price_Small  = 5, 
                    Price_Medium = 7,
                    Price_Large  = 10,
                    Qty_Cake     = xVCake,
                    Qty_Choc     = xVChoc,
                    Price_Cake   = 10,
                    Price_Choc   = 3,
                    Qty_Penalty  = 50)

head(Xtable)

#explanaition
The author is using the Poisson distribution and then the binomial distribution for different parts of the simulation because they are modeling different types of events with 
different characteristics.
Poisson Distribution:
The Poisson distribution is commonly used to model the number of events that occur in a fixed 
interval of time or space. It is appropriate when events happen independently and at a constant average rate. 
In this case, the author is using the Poisson distribution to simulate the daily quantities of small, medium, 
and large bouquets (xVS, xVM, and xVL). These quantities are continuous, non-negative values representing the 
number of bouquets sold each day, and they are likely influenced by factors that follow a Poisson distribution 
(e.g., the number of customers, demand fluctuations, etc.).
 Binomial Distribution:
The binomial distribution is used to model the number of successes in a fixed number of independent 
Bernoulli trials, where each trial has the same probability of success. In this context, the author is 
using the binomial distribution to simulate the daily quantities of cake and chocolate orders (xVCake and xVChoc). 
The variables XCake and XChoc are likely probabilities representing the likelihood that a customer will order 
cake or chocolate on a given day. The binomial distribution is appropriate here since each day order of
cake or chocolate is like a binary outcome (order or no order), and it is based on a probability of occurrence.

# Step 3 - Generate Simulated Y & Analyse
# Apply Function To Simulated Xs

Profit = SimulateY(GetProfit, Xtable)
will take every row in the x table, apply the get profit function and get an output for all the rows.

hist(Profit)

cat("\n", "Avg Daily Profit    :", mean(Profit), "\n",
    "Total Annual Profit :", sum(Profit), "\n","\n")

# Step 3 - Generate Simulated Y & Analyse
#          Apply Function To Simulated Xs

# What is the chance the Florist makes atleast $300/day
length(Profit)
length(Profit[Profit >= 300]) / length(Profit)

# Export Data Into A Spreadsheet
cbind(Xtable,Profit)[1:5,]
write.csv(cbind(Xtable,Profit), "out.csv")

# Exercise - Plot Change In Profit Histogram
#      if the price of Small Boquets (PS)
#      is increased by 1% 

Xtable1 = Xtable;
Xtable1$Price_Small = Xtable$Price_Small * 1.01

Profit0 = SimulateY(GetProfit, Xtable)
Profit1 = SimulateY(GetProfit, Xtable1)

# Task 5 - Sensitivity Analysis
#          When Xs are changed by 1%
#          Find % change in Y 
#          Find the most sensitive Xs 
#----------------------------------------------------------

# METHOD A - Sensitivity of Transfer Function
#            Determined at a given X values

# Step1  - Define the Transfer Function
GetProfit <- function(VS,VM,VL, PS,PM,PL, VCake,VChoc, PCake,PChoc, Max) {
  Profit = VS*PS + VM*PM + VL*PL + 
    VCake*PCake + VChoc*PChoc
  if(VS+VM+VL > Max) Profit = 0.9 * Profit
  return(Profit)
}

# Step 2 - Define a List of Variables with Values
Xvals = list(QtySmall = 20, 
             QtyMedium   = 15,
             QtyLarge    = 10,
             PriceSmall  = 5, 
             PriceMedium = 7,
             PriceLarge  = 10,
             QtyCake     = 0.05*(20+15+10),
             QtyChoc     = 0.10*(20+15+10),
             PriceCake   = 10,
             PriceChoc   = 3,
             QtyPentalty = 50)

# Step 3 - Run the Sensitivity Analysis Function
#          Calculated At The Current X values

S = MonteCarloSens(GetProfit, Xvals)
print(S)
 it will show the %changes when yu change that variable in the other variables.

# METHOD B - Sensitivity of Transfer Function
#            Over the range of Simulated Xs 

# Step 1 - Define the Transfer Function

# Step 2 - Read Data File with X Distribution Data
Xparm = read.csv("Monte Carlo Analysis.csv")
print(Xparm)

# Step 3 - Convert Xparm data into Xtable (Simulated Xs)
#          Xtable can be manually created as in Task 4
Xtable = SimulateX(1000, GetProfit, Xparm)
head(Xtable)


# Step 4 - Simulated Y Values (Optional)
SimulateY(GetProfit,Xtable)

# Step 5 - Run the Sensitivity Analysis Function
#          Calculated At All Simulated Xs
#          Boxplots give a range of sensitivity values
#          Dark line in the centre is Median
#          Coloured box -  Inter Quartile Range (IQR)
#          Whiskers (Q1 - 1.5*IQR) to (Q3 + 1.5*IQR)

S = MonteCarloSens2(GetProfit, Xtable)
print(S)

# Sensitivity of Loan/ Mortgage Formula
# Input - Loan Amount, Term (Months), Rate (Annual)
# Output - Monthly Repayment

MonthlyPayment(10000,36,5)

Xval = list(Loan = 10000,
            Term = 36,
            Rate = 5.0)

MonteCarloSens(MonthlyPayment, Xval)


hist(Profit1 - Profit0)
100*mean(Profit1 - Profit0)/mean(Profit0)
