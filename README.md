[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# CacheTest

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Symmetric separable convex resource allocation problems with structured disjoint interval bound constraints](https://doi.org/10.1287/ijoc.2023.0263) by M. H. H. Schoot Uiterkamp. 


## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0263

https://doi.org/10.1287/ijoc.2023.0263.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{CacheTest,
  author =        {M. H. H. Schoot Uiterkamp},
  publisher =     {INFORMS Journal on Computing},
  title =         {Symmetric separable convex resource allocation problems with structured disjoint interval bound constraints},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0263.cd},
  url =           {https://github.com/INFORMSJoC/2023.0263},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0263},
}  
```

## Description

The software in this repository is developed to solve symmetric separable convex resource allocation problems with structured disjoint interval bound constraints using an exact algorithm.


## Building / Replicating

The exact algorithm is contained in the Python module Alg_disjoint.py.

To reproduce the numerical results of the paper, the script Comparison_RAP_DIBC.py can be run. This requires the Gurobi Optimizer for Python to be installed beforehand. The script calls both Gurobi and the exact algorithm in Alg_disjoint.py.

Warning: running the second part of the script (scalability analysis) has a high total running time.




## Results

The raw output and figures in the paper are in the subdirectory "Results". This contains:

- For the evaluation on EV instances: one file for each considered charging requirement (25, 50, or 100 percent);
- For the scalability evaluation: one file for each combination of considered number of disjoint intervals (m = 2, 3, or 4) and algorithm (exact or (g)urobi); 
- Figure 1: Boxplots of the ratios of the execution times between Gurobi and the exact algorithm for each charging requirement;
- Figure 2: Execution times of the exact algorithm (circles, black) and Gurobi (triangles, gray) where the number of disjoint intervals is 2 (Figure 2a), 3 (Figure 2b), or 4 (Figure 2c).