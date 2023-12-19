#https://mi21.vsb.cz/sites/mi21.vsb.cz/files/unit/numericke_metody_interaktivne.pdf

#1. Seminar
jacobi <- function(A, b, tol=1e-8, max_iter=1000) {
  n <- length(b)
  xold <- numeric(n)
  xnew <- numeric(n)
  
  for (k in 1:max_iter) {
    for (i in 1:n) {
      xnew[i] <- (b[i] - sum(A[i, -i] * xold[-i])) / A[i, i]
    }
    
    if (sum(abs(xnew - xold)) < tol) {
      return(xnew)
    }
    
    xold <- xnew
  }
  
  cat("Jacobi method did not converge within the specified tolerance and maximum iterations.\n")
  return(NULL)
}

A <- matrix(c(1, 0.5, 0.2, 0.3), nrow=2)
b <- c(1, 2)
result <- jacobi(A, b)

if (!is.null(result)) {
  cat("Approximate solution:", result, "\n")
}

solve(A,b)



Lagrange <- function(t,x,y){
  suma <- 0
  for(j in 1:n){
    pitko <- 1
    for (i in 1:n){
      if(i != j) pitko <- pitko *(t - x[i] / x[j] - x[i])
    }
    suma <- suma + y[j] * pitko
  }
  return(suma)
}

n <- 50
x <- 1:n
y <- runif(n)
plot(x,y)
t <- seq(1,10,0.001)

ft <- Lagrange(t,x,y)
lines(t, ft, col="red")

#2. Seminar

MidPointRule <- function(f, a,b, n=1) {
  h <- (b - a) / n
  return(h * sum(f(a + (1 : n) * h / 2)))
}

MidPointRule(sin, 0, pi/2, 1000)


MCIntergralI <- function(f, a, b, ymax, n=1){
  x <- runif(n, a, b)
  y <- runif(n, 0, ymax)
  ymax * (b - a) *sum(y < f(x)) / n
}

MCIntergralI(sin, 0, pi / 2, 1, 10000)


MCIntergralII <- function(f, a, b, n=1){
  return((b - a)*mean(runif(n,a,b)))
}

MCIntergralII(sin, 0, pi / 2, 100000)
#----------------------------------------------

n <- 100
x <- pi * (1:n) / n
y <- 4 * x + runif(n, -1, 1)

plot(x,y)
a <- sum(x * y) / sum(x*x)
abline(a=0, b=a, col='red', lw=2)


EulerStep <- function(f, x, y, h){
  return(y + h * f(x,y))
}

f <- function(t, N, k){return(-k*N)}

N0 <- 1000
dt <- 0.1
t <- seq(0,10, dt)
n <- length(t)
N <- numeric(n)
N[1] <- N0
k <- 0.74
for( i in 2:n) N[i] <- 
  EulerStep(function(t,N) f(t,N,k), t[i-1], N[i-1], dt)

plot(t,N)
  plot(function(N) N0*exp(-k*N), xlim=c(0,10), col='red', add=TRUE, lw=2)
  
#3. Seminar
NewtonInterpol <- function(x,y)
{
  n <- length(x)
  a <- y
  if(n > 1)
  {
    for (i in 2:n)
    {
      sumicka <- 0
      nasobic <- 1
      for(j in 1:(i-1))
      {
        sumicka <- sumicka + a[j] * nasobic
        nasobic <- nasobic * (x[i] - x[j])
      }
      a[i] <- (y[i] - sumicka) / nasobic
    }
  }
  return(a)
}

NewtonPolynomEvaluate <- function(t, a, x)
{
  y <- a[1]
  nasobic <- 1
  for (i in 2:length(a))
  {
    nasobic <- nasobic * (t - x[i-1])
    y <- y + a[i] * nasobic
  }
  return(y)
}
  
n <- 3
x <- 1:n
y <- c(1,5,17)

a <- NewtonInterpol(x,y)
y <- NewtonPolynomEvaluate(x[1], a, x)
y
y <- NewtonPolynomEvaluate(x[2], a, x)


#4. Seminar

SimpsonRule <- function(f,a,b,n=1)
{
  h <- (b-a)/(2*n)
  res <- f(a)+f(b)+4*sum(f(a+2*h*1:n-h))
  if(n>1) res <- res + 2*sum(f(a+2*h*1:(n-1)))
  return(res*h/3)
}


print(sin(1)-sin(0))
print(SimpsonRule(cos,0,1, 1000))


#5. Seminar

f <- function(x){
  return (2*x^3/2*x/5  )
}

g1 <- function(x){
  return ((((2*x+5)/2)^1/3.))
}

g2 <- function(x){
  return (x^3-2.5)
}

g <- function(x){
  return (cos(x))
}

plot(function(x) x, xlim=c(0,2))
plot(g1, xlim=c(0,2), add=TRUE, col='blue', lw=2)
x0 <- 0

for (i in 1:15){
  print(x0)
  x1 <- g(x0)
  segments(x0,x0,x0,x1, col='red')
  segments(x0,x1,x1,x1, col='red')
  x0 <- x1
}

# 6. Seminar

# Gausova eliminace

gauss_jordan_elimination <- function(A, b) {
  # P�ipoj�me b k A, aby se na n�j vztahovaly stejn� ��dkov� operace
  Ab <- cbind(A, b)
  
  n <- length(b)
  # P�eveden� na horn� troj�heln�kovou formu
  for (i in 1:n) {
    # Normalizace ��dk�
    a <- Ab[i, ]
    Ab[i, ] <- a / a[i]
    # Eliminace v�ech ostatn�ch prvk� ve sloupci
    for (j in c(1:n)[-i]) {
      Ab[j, ] <- Ab[j, ] - Ab[j, i] * Ab[i, ]
    }
  }
  
  # Zp�tn� chod - z�sk�n� �e�en�
  for (i in n:1) {
    for (j in (1:n)[-i]) {
      Ab[j, ] <- Ab[j, ] - Ab[j, i] * Ab[i, ]
    }
  }
  
  # V�sledek je v posledn�m sloupci po transformaci
  return(Ab[, n + 1])
}

# P��klad pou�it�
A <- matrix(c(2, 1, -1, -3, -1, 2, -2, 1, 2), nrow = 3, byrow = TRUE)
b <- c(8, -11, -3)

result <- gauss_jordan_elimination(A, b)
print(result)

# 7. Seminar

Hermite <- function(x0,x1,y0,y1,d0,d1){
  h <- x1-x0
  hrec <- 1/h
  delta <- (y1-y0)*hrec-d0
  Delta <- d1-d0
  a3 <- (Delta - 2*delta)
  return(c(y0,d0,(delta-a3)*hrec,a3*hrec*hrec))
}
Horner <- function(x,x0,coef){
  n <- length(coef)
  y <- coef[n]
  t <- x-x0
  for(i in (n-1):1) y <- y*t+coef[i]
  return(y)
}
fd <- function(f,x,h=0.0000001){
  return((f(x+h)-f(x-h))/(2*h))
}

x0 <- -1
x1 <- 3
y0 <- 7
y1 <- -2
d0 <- -2
d1 <- 3

coef <- Hermite(x0,x1,y0,y1,d0,d1)
plot(function(x) Horner(x,x0,coef), xlim=c(x0-0.5,x1+0.5))
lines(c(x0,x1), c(y0,y1), col="red")

print(fd(function(x) Horner(x,x0,coef), x0)) #d0
print(fd(function(x) Horner(x,x0,coef), x1)) #d1

#8. Seminar

ODEstep <- function(f, y, t, dt) {
  return(y + dt * f(y, t))
}

SIR <- function(y, t, beta = 0.4, nu = 0.2) {
  dS <- -beta * y[1] * y[2]
  dR <- nu * y[3]
  dI <- -dS - dR
  return(c(dS, dI, dR))
}

n <- 1000
Y <- matrix(0, n, 3)
Y[1, ] <- c(0.999, 0.001, 0)

for (i in 2:n) {
  Y[i, ] <- ODEstep(SIR, Y[i - 1, ], i, 0.1)
}

plot(1:n, Y[, 1], type='l', col='blue', ylab='Population', xlab='Time', main='SIR Model')
lines(1:n, Y[, 2], col='red')
lines(1:n, Y[, 3], col='green')
ylim <- c(0, 1)
axis(2, ylim=ylim)



ODEstep <- function(f, y, t, dt) {
  return(y + dt * f(y, t))
}

ODEstepMiddle <- function(f, y, t, dt) {
  dtpul <- 0.5 * dt
  y_mid <- ODEstep(f, y, t, dtpul)
  return(y + dt * f(y_mid, t + dtpul))
}

ODEstepA <- function(f, y, t, dt) {
  y_mid <- ODEstep(f, y, t, dt)
  return(y + dt * (f(y, t) + f(y_mid, t + dt)) / 2)
}

LV <- function(y, t, r = 0.67, a = 1.33, s = 1, b = 1) {
  return(c(y[1] * (r - a * y[2]), y[2] * (-s + b * y[1])))
}

n <- 5000
Y <- matrix(0, n, 2)
Y[1, ] <- c(0.9, 0.9)

for (i in 2:n) {
  Y[i, ] <- ODEstepA(LV, Y[i - 1, ], i, 0.02)
}

plot(1:n, Y[, 1], type='l', col='blue', ylab='Population', xlab='Time', main='Lotka-Volterra Predator-Prey Model')
lines(1:n, Y[, 2], col='red')
legend("topright", legend=c("Rabbits", "Foxes"), col=c("blue", "red"), lty=1)



TK <- function(y, t, o = 4, z = 0.05) {
  return(c(y[2], -(2 * z * o * y[2] + o^2 * y[1])))
}

n <- 1000
Y <- matrix(0, n, 2)
Y[1, ] <- c(10, -10)

for (i in 2:n) {
  Y[i, ] <- ODEstepA(TK, Y[i - 1, ], i, 0.02)
}

plot(1:n, Y[, 1], type='l', col='blue', ylab='Position', xlab='Time', main='Damped Harmonic Oscillator')


#9. Seminar
### Metoda Runge-Kutta 4. Řádu

ODEstepRK4 <- function(f,t,y, dt)
{
  dtpul <- 0.5 * dt
  tp <- t * dtpul
  k1 <- f(t,y)
  k2 <- f(t+dtpul, y+k1*dtpul)
  k3 <- f(tp, y+k2*dtpul)
  k4 <- f(t+dt, yd+dt*k3)
  return(y+dt*(k1+2*k2+2*k3+k4)/6)
}

f <- function(x)
{
  return(exp(-x*x)*2/sqrt(pi))
}

trapezodial <- function(f,a,b,n=1)
{
  h <- (b-a)/n
  suma <- f(a)+f(b)
  if (n>1)
  {
    suma <- suma + 2*sum(f(h*(1:(n-1))))
  }
  return(0.5*h*suma)
}

options(digits = 8)
m <- 4
a <- sapply(2**(0:m), function(n) trapezodial(f,0,1,n))
for (j in 1:m)
{
  a <- (4**j*a[2:(m+2-j)]-a[1:(m+1-j)])/(4**j-1)
}
print(a)