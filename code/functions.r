library(pracma)

a <- rep(0, 250)
b <- rep(0, 250)
rv <- 51:300

for (r in rv){
	s <- function(x){
		n <- exp(-x^2)-exp(-(r-x)^2)+sqrt(pi)*x*(erf(r-x)+erf(x))
		n
	}

	x <- 1 + 1:r * 0.3

	y <- s(x)

	l <- lm(y ~ x)

	b[r-50] <- l$coefficients[1]
	a[r-50] <- l$coefficients[2]
}

r <- 100
x <- 1 + 1:r * 0.3
s <- function(x){
		n <- exp(-x^2)-exp(-(r-x)^2)+sqrt(pi)*x*(erf(r-x)+erf(x))
		n
}
s1 <- function(x){
		n <- 3.54*x + 0.187/(r-3.156)
		n
}

plot(x, s(x), type="l", col="blue", lwd=2)
lines(x, s1(x), col="red", lwd=2, lty=2)

f <- function(r){r^3*(3*log(r)-1)}
r0 <- 0.2
r <- 0:100/100*r0
t <- f(r)-f(r0)

#plot(t, r, type="l", col="blue", lwd=2)

