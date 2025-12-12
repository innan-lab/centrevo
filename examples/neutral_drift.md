# Scenario: Neutral Drift

This scenario simulates the evolution of centromeric sequences under "neutral drift," meaning there is no selection pressure favoring specific sequence properties. Evolution is driven solely by random mutations and recombination events.

In the absence of selection or "drive," centromeric arrays often maintain a relatively stable structure or decay slowly over very long timescales, depending on the balance between recombination (homogenization) and mutation (divergence).

## 1. Initialize Simulation

We set up a small population (N=50) for a short duration (200 generations) to observe the baseline dynamics.

```bash
centrevo init \
  -N neutral_drift \
  -n 50 \
  -g 200 \
  -o neutral_drift.db
```

## 2. Run Simulation

Execute the simulation. Since this is a small run, it should finish in a few seconds.

```bash
centrevo run -N neutral_drift -d neutral_drift.db
```

## 3. Expected Outcome

- **Diversity**: Check `nucleotide_diversity` (pi). It should remain relatively low or increase slowly due to mutation drift.
- **Fitness**: Since there is no selection, fitness will remain constant (1.0) for all individuals.
- **Structure**: The repeat structure should remain largely intact, with random point mutations accumulating over time.
