# bdtree

Birth-death tree models for [BEAST 3](https://github.com/CompEvol/beast3).

Implements the birth-death-sequential-sampling (BDSS) model for tree likelihood and simulation (Stadler & Yang, 2013).

**Contributors:** Fabio K. Mendes, Rong Zhang

## Models

| Class | Description |
|-------|-------------|
| `BirthDeathSequentialSampling` | BDSS tree likelihood with optional fossil sampling |
| `BirthDeathSerialSamplingTree` | Tree simulator under the BDSS model |

## Building

BEAST 3 dependencies are resolved from [Maven Central](https://central.sonatype.com/namespace/io.github.compevol) — no extra configuration needed.

```bash
mvn compile
mvn test
```

To develop against an unreleased SNAPSHOT, install BEAST 3 from source:

```bash
cd ~/Git/beast3
mvn install -DskipTests
```

## Running

```bash
# Run an analysis
mvn exec:exec -Dbeast.args="examples/testing/BDSSLikelihood.xml"
```

## Examples

- `examples/testing/BDSSLikelihood.xml` — MCMC analysis with BDSS tree prior (10 taxa, fossils)
- `examples/testing/BDSSTreeSimulator.xml` — Simulate trees from the BDSS model
- `examples/testing/Shankarappa.xml` — HIV sequence data analysis

## References

Stadler, T., & Yang, Z. (2013). Dating phylogenies with sequentially sampled tips. *Systematic Biology*, 62(5), 674-688. https://doi.org/10.1093/sysbio/syt030
