data {
  int<lower=1> J;
}
parameters {
  vector[J] x;
}
model {
  x ~ normal(0, 1);
}
