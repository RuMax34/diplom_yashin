bet <- function(theta){
	cot <- 1/tan(theta)
	b <- sqrt(1+cot^2)-cot
	1.5*b*(1+b^2/3)
}
# Angle between 0 and pi/2, no more!

theta <- pi/4
beta <- bet(theta)
kT <- 0.026*550/300
FAs <- 0.1
DGa <- 1
DAs <- 1
tAs <- exp(0.1/kT)
C0Ga <- 2
C0As <- FAs*tAs
OGaAs <- 1
OGa <- 0.1
kr <- 0.1
Rd0 <- 20
Rinf <- 70
dr <- 0.4
dt <- 0.05
Nr <- floor(Rinf/dr)
AGa <- DGa*dt/dr^2
AAs <- DAs*dt/dr^2
K <- kr*dt
F <- FAs*dt
S <- dt/tAs
P <- 4*dt/beta*OGa*DGa/dr
Q <- OGaAs*K

qd <- floor(Rd0/dr)

CGa1 <- rep(0,Nr)
CGa1[1:qd] <- C0Ga
CGa2 <- CGa1
CAs1 <- rep(C0As,Nr)
CAs1[1:qd] <- 0
CAs2 <- CAs1

h <- rep(0, Nr)
plot((1:(Nr-1))*dr, h[-Nr], pch = NA, xlab = "Расстояние, ячейки", ylab = "Высота кольца, ячейки", ylim = c(0, 40))

Rd1 <- Rd0
Rd2 <- Rd0
Rd <- Rd0

while (Rd2 > 3) {
	for (j in qd:(Nr-1)){
		CGa2[j] <- CGa1[j] + AGa*((1-0.5/j)*CGa1[j-1] - 2*CGa1[j] + (1+0.5/j)*CGa1[j+1]) - K*CGa1[j]*CAs1[j]
		CAs2[j] <- CAs1[j] + AAs*((1-0.5/j)*CAs1[j-1] - 2*CAs1[j] + (1+0.5/j)*CAs1[j+1]) - K*CGa1[j]*CAs1[j] + F - S*CAs1[j]
	}
	h <- h + Q*CGa1*CAs1
	Rd2 <- sqrt(Rd2^2 + P*(CGa1[qd]-CGa1[qd-1]))
	if (Rd1 - Rd2 > 1) {
		Rd1 <- Rd2
		lines((1:(Nr-1))*dr, h[-Nr], col = "blue", lwd = 2, lty = sample(1:6, 1))
	}
	Rd <- c(Rd, Rd2)
	qd <- floor(Rd2/dr)
	CGa1 <- CGa2
	CAs1 <- CAs2
}

