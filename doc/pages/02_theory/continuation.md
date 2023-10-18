---
title: Continuation method
ordered_subpage: continuation.md
---

All the phase envelopes in this software are calculated with the continuation
method, as described by [Allgower and Georg](https://epubs.siam.org/doi/book/10.1137/1.9780898719154).

It is used to calculate lines and hyper-lines defined by a system of equations.
Using information from the last calculated point, the next point is
initializated with a reasonable estimation of its solution, making convergence
both faster and easier in hard to converge problems.

This method can be simplified in single steps

1. Specify the value of a single variable.
2. Solve a system of equations at the specification value.
3. Determine how all variables change with respect that specified variable.
4. Select a new specification with that information.

The utilization of the contiuation method in phase-equilibria diagrams was 
shown by [Michelsen](https://www.sciencedirect.com/science/article/pii/037838128080001X)
with the calculation of multicomponent phase-envelopes.