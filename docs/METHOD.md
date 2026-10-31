# Method

## System

A dual-rail logical qubit is one excitation shared by two modes, `|a>` and `|b>`; `|0>` (no excitation)
is the erased state. The Lindblad master equation (time in units of the mean lifetime 1/gamma_bar) has

* loss `L_a = sqrt(1 + eta) |0><a|`, `L_b = sqrt(1 - eta) |0><b|` (asymmetry eta),
* differential dephasing at rate gamma_phi (default 0.01),
* optional quasi-static mode detuning delta drawn from a Gaussian of standard deviation sigma_delta
  (averaged exactly with a 7-node Gauss-Hermite rule),
* a swap pulse: beam-splitter coupling for a time t_pi (default 0.003) with fractional over-rotation
  eps, loss and dephasing acting during the pulse.

A hold of total duration T with n pulses uses Carr-Purcell-Meiboom-Gill (CPMG) spacing
`tau/2, pi, tau, pi, ..., pi, tau/2` and keeps the total elapsed time fixed at T, so all methods
see the same decay. n must be even so the net logical map is the identity.

## Methods compared

1. **Standard quantum baseline: bare erasure postselection.** Hold, then discard flagged shots.
2. **Classical baseline: calibrated no-jump inversion.** Estimate the survival of `|a>` and `|b>`
   from 20000 calibration shots each, and undo `K = diag(exp(-x_a/2), exp(-x_b/2))` on the
   estimated logical Bloch vector (clipped to the Bloch ball). An oracle variant uses the exact rates
   at run time and isolates the shot-noise floor.
3. **Proposed: swap echo** with n = 2, 4, 8 pulses (Floquet average of the two loss rates).

## Metrics

* Conditional infidelity (exact): mean over the six Pauli eigenstate inputs of one minus the overlap
  of the renormalised accepted state with the input; acceptance probability.
* Bloch-vector mean squared error (MSE, shot-limited): 20000 shots per Pauli basis and input,
  flagged shots discarded, squared distance of the estimated Bloch vector to the ideal input,
  averaged over inputs and 200 repetitions per seed.

## Seeds and randomness

Five seeds. Each seed draws a device pulse over-rotation eps ~ Normal(0, 0.01); each repetition
draws a drift of the loss asymmetry (standard deviation 0.1) between calibration and run, plus
binomial shot and calibration noise. Phase maps and ablations are deterministic (exact).

## Experiments

* **E1 (H1)** infidelity versus eta for three hold times against (2 eta T)^2 / 24.
* **E2 (H2)** two phase maps, log10 ratio echo/bare, over (eta, T) and (eps, T), no detuning so that
  only the non-Hermitian term is probed.
* **E3 (H3)** methods table at T = 0.1, 0.3, 0.6; a 10-point hold sweep for the cost-quality Pareto
  plot; a drift-versus-over-rotation map of echo versus calibrated inversion at T = 0.3.
* **Ablation** ideal pulses, over-rotation only, finite pulse duration only, both; with and without
  quasi-static detuning. Pulse-count scaling n = 0...32.
* **Noise sweep** dephasing rate gamma_phi = 0...0.2 (share of the bare error that is bias),
  over-rotation and drift axes of the maps above.
