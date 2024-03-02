# Measuring patterns separability in complex data and network community embedding through the traveling salesman path

This repository offers a framework for the definition and measure of the geometric separability (linear and nonlinear) of mesoscale patterns in complex data visualization by solving the traveling salesman problem.

## Before starting

This software requires [Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) to execute calculations. Check the [requirements](REQUIREMENTS.md) and follow the instructions to install this dependency in your operating system.

## Usage

```matlab
[measures, metadata] = CommunitySeparability(embedding, communities, variant, 'positives', positiveCommunities, 'permutations', numberOfPermutations)
```

### Inputs

- `embedding` (required): Numeric matrix of at least two dimensions. Eventually, this input is the output obtained after applying a dimension reduction or network embedding method, where the rows represent the samples and the columns the dimensions.
- `communities` (required): Numeric array or cell array of strings containing the sample labels (i.e., classes or groups).
- `variant` (required): String representing the desired community separability index to be calculated. This value should be one of: `tsps` for a traleving salesman problem separability computation, `cps` for a centroid projection separability computation, or `ldps` for a linear discriminant projection separability computation.
- `positives` (optional): Numeric or cell array of strings containing the positive sample labels (i.e., positive classes or groups).
- `permutations` (optional): Number of desired permutations to compute a null-model based on the selected variant.

### Outputs

- `measures`: Structure whose elements are the statistical-based measures to assess (i.e., quantify) the separability of the provided inputs. If `permutations` are inputted, then, this structure will contain a substructure per measure with the results of the null-model evaluation.
- `metadata`: Structure whose elements are the different pairwise comparisons used to evaluate the separability of the provided inputs.

For more details about the described inputs and outputs, see [Example1](example1.m) and [Example2](example2.m).

## Reporting an issue

For reporting problems, questions, and suggestions; please, use the [Issues](https://github.com/biomedical-cybernetics/travelling-salesman-path/issues) section.
