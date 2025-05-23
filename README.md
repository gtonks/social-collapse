# social-collapse
Stochastic modeling of social collapse

## Contents

### `dhanalysis`
Analyzes system behavior under different equations for `dh/dt`.

Equation 1: `dh/dt = 0`.

Equation 2: `dh/dt = mu_h*x*h*(1 - h/k_h)`.

Equation 3: `dh/dt = mu_h*h - c_h*h/x`.

Equation 4:
```
dh/dt = h*(r - gamma(x)*h)
gamma(x) = beta - (beta - alpha)*x/k_gamma.
```

#### Execution
From this directory,
```
python -m socialcollapse.dhanalysis.eq1
python -m socialcollapse.dhanalysis.eq2
python -m socialcollapse.dhanalysis.eq3
python -m socialcollapse.dhanalysis.eq4
```

### `eulermaruyama`
Contains Euler-Maruyama simulations with a fixed `dH/dt = 0`.

#### Execution
From this directory,
```
python -m socialcollapse.eulermaruyama.H_var
python -m socialcollapse.eulermaruyama.mu_var
```

### `gillespie`
Contains old attempts at Gillespie simulations.