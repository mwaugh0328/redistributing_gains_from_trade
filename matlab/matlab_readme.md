This file describes the structure of the matlab code used to compute the results from ["Redistributing the Gains From Trade Through Progressive Taxation"](http://www.waugheconomics.com/uploads/2/2/5/6/22563786/lw_tax.pdf)

Contents and description of code...

- ``island_solve_calibrate_NAG.m`` which takes a progressivity parameter and then calibrates and solves for an equilibrium to the economy.

- ``clustergrid`` which creats an asset grid with points clusterd around places of interest.

- ``open_economy_markovchain`` creats the markov chain which characterizes the joint dynamics of local productivity shocks and world prices. Uses ``rouwenhorst.m``

-

---

### Calibration

Simply run ``island_solve_calibrate_NAG.m`` with the argument 0.18. To adjust the calibration targets, lines 16-20 contain the other data moments.

---

### Compute Optimal Policy

---

### Plotting
