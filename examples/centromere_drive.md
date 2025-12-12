# Scenario: Centromere Drive (Length Selection)

"Centromere Drive" is a hypothesis where centromeres compete during female meiosis to be transmitted to the egg. Larger or "stronger" centromeres often have a transmission advantage.

In this simulation, we model this "drive" by applying positive selection on **chromosome length**. We set a target length that is larger than the starting length, encouraging the population to evolve longer arrays.

## 1. Initialize Simulation

- **Target Length (`--fit-len-opt`)**: 250,000 bp (Starting length is ~205,000 bp).
- **Tolerance (`--fit-len-std`)**: 10,000 bp (Standard deviation of the fitness function).

```bash
centrevo init \
  -N drive \
  -n 50 \
  -g 200 \
  --fit-len-opt 250000 \
  --fit-len-std 10000 \
  -o drive.db
```

## 2. Run Simulation

```bash
centrevo run -N drive -d drive.db
```

## 3. Expected Outcome

- **Array Expansion**: The average chromosome length in the population should increase over time as longer haplotypes are favored.
- **Fitness Increase**: Mean population fitness should rise as the population adapts to the driven length optimum.
- **Copy Number Variation**: Mechanistically, this expansion is driven by recombination events (like unequal crossing over or gene conversion) that duplicate repeat units.
