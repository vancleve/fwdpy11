Simulating quantitative traits, I.
==========================================

fwdpy11 allows the simulation of quantitative traits in a fairly general way.  These simulations differ from standard
population genetic simulations in that a "trait value" is calculated and that a mapping from trait value to fitness is
required. Optionally, one may want to add "noise" to trait values, perhaps representing non-genetic contributions to
phenotype.

This page covers simulations of quantitative traits in a single genomic region (*e.g.* using
:class:`fwdpy11.fwdpy11_types.Spop`).

Let's call our trait value :math:`P`, which is composed of a genetic component, :math:`G` and a random component
:math:`E`.

Trait values
-----------------------------

The following trait value functions are implemented:

* :class:`fwdpy11.trait_values.SpopAdditiveTrait`
* :class:`fwdpy11.trait_values.SpopMultTrait`
* :class:`fwdpy11.trait_values.SpopGBRTrait`



Adding noise to trait values
----------------------------------------------------------

Mapping trait values to fitness
----------------------------------------------------------


