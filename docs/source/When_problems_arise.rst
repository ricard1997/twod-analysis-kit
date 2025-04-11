When problems arise
====================


Customizable lipids
+++++++++++++++++++


Customizing lipid heads
#######################

Our code will work perfectly and with zero to none user intervention as long as it is composed by cholesterol,
DODMA, and phospholipids that contains a P in their headgroup. For others lipids, including customized ones, we provide
an easy way of using them. At initialize :code:`Voronoi2D`, :code:`Cumulative2D` or :code:`PackingDefects`,
the user should provide the headgroup they want to be taken into account for each lipid by passing :code:`working_lip`.
For example, for a membrane composed by "DSPC", "DODMA", "CHL1", and "POPE", we should provide






.. code-block:: python

    import MDAnalysis as mda
    from twodanalysis import Cumulative2D
    from twodanalysis import Voronoi2D

    # Dictionary containing the lipid heads to be used
    working_lip = {
        "DSPC": {"head" : "P"},
        "DODMA" : {"head" : "N"},
        "CHL1" : {"head": "O3"},
        "POPE" : {"head": "P"}
    }


    # Start Cumulative 2D or Voronoi2D
    tpr = "membrane.tpr"
    xtc = "membrane.xtc"
    univ = mda.Universe(tpr, xtc)
    membrane = Cumulative2D(univ, working_lip = working_lip)

    # or

    membrane = Voronoi2D(univ, working_lip = working_lip)


.. note::
    This method can be used to include lipids that give errors, user defined lipids, and Coarse grain lipids. Furthermore,
    The user can define any atom here, which allows us to use any atom from the lipids and for the corresponding 2D projections.
    Lately, the user can also define multiple atoms. In such case, the algorithm take the center of mass of all the atoms as the head
    (This is an experimental feature which could not work in all analysis).


Customizing lipid tails
#######################

In cases such as splay angle, order parameters and packing defects, in addition to head information we need information of the tail.
Our code is able to guess this information for most standard lipids. However, there may be some cases where our code does not work properly.
For instance, when dealing with sterols, lipids, with one tail or three tails. In such case, the user should provide the bond that connects
the tail with the head. Follows an image rendered with VMD showing the bonds for CHL1, DODMA, and DSPC.

 .. image:: connection.png

The image above shows the names for the atoms that belong to the bond that connects the lipid tails with the headgroup for different lipids.
These bonds (atoms names) should be provided to the :code:`connection_chain` attribute so 2Danalysis process them correctly. Follows an example
for the mentioned lipids.


.. code-block:: python

    import MDAnalysis as mda
    from twodanalysis import Cumulative2D
    from twodanalysis import Voronoi2D

    # Dictionary containing the connection of the chains
    connection_chains = {

            "CHL1" : [("O3", "C3")],
            "DODMA" : [("C21", "C22"), ("C31", "C32")],
            "DSPC" : [("C21", "C22"), ("C31", "C32")],

        }

    # Start Cumulative 2D or Voronoi2D
    tpr = "membrane.tpr"
    xtc = "membrane.xtc"
    univ = mda.Universe(tpr, xtc)
    membrane = Cumulative2D(univ, connection_chains = connection_chains)

    # or

    membrane = Voronoi2D(univ, connection_chains = connection_chains)


.. note::
    This method can be used to include any lipid and we recommend to add the lipid chains in the order sn1, sn2 for those with two tails. This method
    works for lipids with any number of lipid tails, in case of three lipid tails one should provide a list with 3 bonds. Using this method,
    potentially any lipid can be added to the 2Danalysis framework, including customized lipids and MARTINI lids.

We also offer a nice way to check if the lipids tails are being assigned correctly by plotting them with :code:`visualize_polarity()` which would
output an image as follows:


.. code-block:: python

    membrane.visualize_polarity()
    plt.show()

.. image:: polarity.png

Periodicity
+++++++++++++++++++


All our code include periodicity handling by replicating a percentage of the data, by default 10% in each border. Sometimes,
specially for Voronoi2D, the default percentage of data replication sould not be enough and APL would show big values or look weird. In
such case, the user should increase the periodicity value as follows.

.. code-block:: python

    import MDAnalysis as mda
    from twodanalysis import Cumulative2D
    from twodanalysis import Voronoi2D
    from twodanalysis import PackingDefects


    tpr = "membrane.tpr"
    xtc = "membrane.xtc"
    univ = mda.Universe(tpr, xtc)



    # For Cumulative2D
    membrane = Cumulative2D(univ)
    membrane.periodicity = 0.2 # Increase replication of data to 20%
    # For Voronoi2D
    membrane = Voronoi2D(univ)
    membrane.periodicity = 0.3 # Increase replication of data to 30%

    # For PackingDefects
    membrane = PackingDefects(univ)
    membrane.periodicity = 0.5 # Increase replication of data to 50%

.. note::
    This is not needed unless your images are showing errors. This is definitely not needed if your edges are significantly
    smaller than the periodic box size.