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
Euler-Maruyama simulations.

#### `H_var`
Variation in `H` with an average `dH/dt = 0`.

#### `mu_var`
Variation in `mu` with a fixed `dH/dt = 0`.

#### `mu_var_linear_dH`
Variation in `mu` with `dH/dt = rH - gamma(x)H^2`.

#### `mu_var_nonlinear_dH`
Variation in `mu` with `dH/dt = BxH / (rho + x) - dH - gamma(x)H^2`.

#### `final_search_over_a_mu_var_linear_dH`
Computes data for final points after `n_steps` iterations. Repeats for several `a` and `sigma` over `n_trials`.

#### `plot_prob_mesh_over_a`
Plot the probability of resource availability over the parameter search data from `final_search_over_a_mu_var_linear_dH`.

#### `final_search_over_b_mu_var_linear_dH`
Computes data for final points after `n_steps` iterations. Repeats for several `b` and `sigma` over `n_trials`.

#### `plot_prob_mesh_over_b`
Plot the probability of resource availability over the parameter search data from `final_search_over_b_mu_var_linear_dH`.

#### Execution
From this directory,
```
python -m socialcollapse.eulermaruyama.H_var
python -m socialcollapse.eulermaruyama.mu_var
python -m socialcollapse.eulermaruyama.mu_var_linear_dH
python -m socialcollapse.eulermaruyama.mu_var_nonlinear_dH
python -m socialcollapse.eulermaruyama.final_search_over_a_mu_var_linear_dH
python -m socialcollapse.eulermaruyama.plot_prob_mesh_over_a
python -m socialcollapse.eulermaruyama.final_search_over_b_mu_var_linear_dH
python -m socialcollapse.eulermaruyama.plot_prob_mesh_over_b
```

### `gillespie`
Contains old attempts at Gillespie simulations.
