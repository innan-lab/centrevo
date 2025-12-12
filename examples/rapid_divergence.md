# Scenario: Rapid Divergence (High Mutation)

This scenario demonstrates what happens when the mutation rate is significantly higher than normal. In biological terms, this could represent a hyper-mutable region or a defect in DNA repair mechanisms.

Without strong purifying selection to weed out mutations, the repetitive structure of the centromere can decay rapidly, leading to a loss of the higher-order repeat (HOR) structure.

## 1. Initialize Simulation

We increase the mutation rate to `1e-4` (an order of magnitude higher than the default `1e-5`).

```bash
centrevo init \
  -N high_mut \
  -n 50 \
  -g 200 \
  --mutation-rate 1e-4 \
  -o high_mut.db
```

## 2. Run Simulation

```bash
centrevo run -N high_mut -d high_mut.db
```

## 3. Expected Outcome

- **Diversity**: Nucleotide diversity (pi) should increase much faster than in the neutral drift scenario.
- **Sequence Identity**: If you inspect the sequences, you will see more "noise" or differences between repeat units.
- **Decay**: Over longer timescales (e.g., thousands of generations), the distinct HOR pattern may become blurred as mutations accumulate faster than recombination can homogenize them.
