data <- read.table("bessel_ratio_values.dat", col.names=c("x", "y"))

fit <- lm(y ~ poly(x, 9, raw=TRUE), data=data)

# Step 3: Show summary
summary(fit)

# Predict and plot
new_data <- data.frame(x = seq(min(data$x), max(data$x), length.out = 100))
new_data$y_pred <- predict(fit, newdata = new_data)

f <- function(x){
    4.657125e-06 + -7.136973e-04*x + 6.389987e-01*x^2 + 1.260996e+00*x^3 + -1.788035e+00*x^4 + 9.097378e+00*x^5 + -2.079507e+01*x^6 + 3.516017e+01*x^7 + -3.297203e+01*x^8 + 1.510697e+01*x^9
}

plot(data$x, data$y-f(data$x), col = "blue", pch = 19)
lines(data$x, data$y-f(data$x), col = "red", lwd = 2)