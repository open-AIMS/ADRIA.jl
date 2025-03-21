# Dynamic Multi-Criteria Decision Analysis

## Multi-Criteria Decision Making

Multi-criteria decision analysis (MCDA) encompasses a series of methods for conducting
decision making in a way that is formalised, structured and transparent. It evaluates the
decision objective according to numerical measures of decision criteria assembled by the
decision maker. Constructing explicit criteria over which to evaluate alternatives allows a
transparent evaluation of benefits, negatives and trade-offs in coming to a decision
solution.

Typical approaches to MCDA require the construction of a "decision matrix", which takes the
form of a $X^{n \cdot m}$, where $n$ is the number of alternate options available, and $m$
is the number of criteria.

In the context of ADRIA, the alternate options relate to the locations being assessed.
The criteria then relate to the common set of attributes on which the locations
are being judged, such as the projected heat stress (in terms of DHW), depth, and level of
incoming or outgoing connectivity.

MCDA methods provide a ranking according to the set of assessed criteria, in the form of:

$$r = g(X, w, d)$$

where:

- ``g()`` refers to a given MCDA algorithm (see [JMcDM.jl](https://jbytecode.github.io)).
- ``X`` is the decision matrix
- ``w`` is the weights afforded to each criteria, indicating their relative importance
- ``d`` is the desired directionality for each criterion (to minimize the criteria value, or to maximize)
- ``r`` is the ranking determined by ``g()``

When applied in conjunction with scenario analyses, locations are assessed at each decision
point.

When conducting location selection, the analyses are applied to the initial conditions
represented in the Domain.
