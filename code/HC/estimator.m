function alpha = estimator(X, Y)

alpha = pinv(X' * X) * (X' * Y)

end