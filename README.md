# Ray-tracing-in-cosmology
In this project, we fired 100,000 light rays in random directions for three different cases to study the effect of inhomogeneities in cosmology. The published paper that explains the codes and its results can be found on https://arxiv.org/pdf/1705.02328.pdf - Ray tracing and Hubble diagrams in post-Newtonian cosmology. We ran each code in parallel over about 80 cores. As each light ray trajectory run was independent, we split up each code to run for 1250 light rays to give us our total of 100,000.

The first code is Final_lambda_check.c. This code was used to test the speed and accuracy of our code. We considered the simple model of a universe with just a cosmological constant (lambda), which has exact known solutions that we can compare our results to.

The second code is Final_point_source_check.c. In this code we considered a universe with a cosmological constant and point masses of about 3 kiloparsecs radius to represent opaque galaxy buldges.

The second code is Final_homogeneous_halo_check.c. In this code we considered a universe with a cosmological constant and homogenous haloes of about 30 kiloparsecs radius to represent transparent dark matter haloes.
