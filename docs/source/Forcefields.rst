.. _forcefields:

Forcefields
===============


Force fields play an important role in membrane simulations, as they define the naming
conventions and structural details for lipids and other membrane components.
``twodanalysis`` currently supports the **CHARMM** and **AMBER** force fields.
However, the package was designed to be as general as possible, so many features will
still work with other force fields if the guidelines in :ref:`When problems arise
<when_problems_arise>` are followed.

We are actively working to extend full support to additional force fields in future
releases. In the meantime, we encourage users to open issues on our GitHub page if
their force field is not recognized or if specific functionality is needed.



We have introduced a simple way to specify the forcefield. You should add the flag forcefield with either "amber" or "charmm"
when initializing the ``Cumulative2D``, ``Voronoi2D`` or ``PackingDefects`` objects, as follows.

For amber:
----------
.. code-block:: python

    memb1 = Cumulative2D(forcefield = "amber")
    memb2 = Voronoi2D(forcefield = "amber")
    memb3 = PackingDefects(forcefield = "amber")

If you have any problems using the amber forcefield, please refer to :ref:`When problems arise <when_problems_arise>`.

For charmm:
-----------
.. code-block:: python

    memb1 = Cumulative2D(forcefield = "charmm")
    memb2 = Voronoi2D(forcefield = "charmm")
    memb3 = PackingDefects(forcefield = "charmm")

If you have any problems using the amber forcefield, please refer to :ref:`When problems arise <when_problems_arise>`.

.. note:: The default forcefield is charmm, if you do not specify the forcefield, charmm will be used automatically.

Other forcefields:
------------------

If you are using a different forcefield, please refer to :ref:`When problems arise <when_problems_arise>`.


.. rubric:: References

- `Ramirez, R. X., Bosch, A. M., Pérez, R., Guzman, H. V., & Monje, V. (2025). 2Danalysis: A toolbox for analysis of lipid membranes and biopolymers in two-dimensional space. Biophysical Journal. <https://www.cell.com/biophysj/fulltext/S0006-3495(25)00321-2>`_
