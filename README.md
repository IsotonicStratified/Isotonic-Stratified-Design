# Isotonic design for a single-arm trial with a biomarker   
In single-arm trials with a predefined subgroup based on baseline biomarkers, it is often assumed that a biomarker defined subgroup, the biomarker positive subgroup, has the same or higher response to treatment compared to its complement, the biomarker negative subgroup. The goal is to determine if the treatment is effective in both groups or in the biomarker positive group only or not effective at all. We propose an isotonic stratified Simon’s design for this problem. The design has a joint set of decision rules for biomarker positive and negative subjects and utilizes joint estimation of response to treatment using assumed monotonicity of response probabilities between the biomarker negative and positive groups. 

## The materials
*Supplementary material.pdf* describes how to calculate the probability of rejecting the null hypothesis through the three routes in the paper.

*exact_pava.R* implements the exact calculation method described in the above file.

*pava simulation.R* simulates the isotonic design and the parallel Simon's design.
