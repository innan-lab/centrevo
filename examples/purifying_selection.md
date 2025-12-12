# Scenario: Purifying Selection (GC Content)

Centromeres play a critical structural role, and their sequence composition (e.g., GC content) can be functionally important. This scenario simulates **purifying selection**, where deviating from a specific GC content reduces an individual's fitness.

This acts as a constraint, preventing the sequence from drifting too far towards AT-rich or GC-rich extremes, which can happen under neutral drift or biased gene conversion.

## 1. Initialize Simulation

- **Optimal GC (`--fit-gc-opt`)**: 0.42 (42% GC content).
- **Concentration (`--fit-gc-conc`)**: 50.0 (Strength of selection; higher = stricter).

```bash
centrevo init \
  -N purifying \
  -n 50 \
  -g 200 \
  --fit-gc-opt 0.42 \
  --fit-gc-conc 50.0 \
  -o purifying.db
```

## 2. Run Simulation

```bash
centrevo run -N purifying -d purifying.db
```

## 3. Expected Outcome

- **Composition Stability**: The mean GC content of the population should remain tightly clustered around 42%.
- **Variance Reduction**: Compared to neutral drift, the variance in GC content between individuals should be lower.
- **Fitness Cost**: Individuals that drift away from 42% GC (due to random mutations) will have lower fitness and produce fewer offspring.
