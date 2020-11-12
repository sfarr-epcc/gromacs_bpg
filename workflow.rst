==========================================
Suggested workflow for GROMACS simulations
==========================================

This part of the GROMACS best practice guide walks you through a suggested
workflow for preparing, running, and analysing the results of your GROMACS
simulations. There are already many very good tutorials available to teach 
new users how to do this. Where applicable, we will refer to these. Here, 
we have collated and compiled the suggestions of these tutorials into what 
we believe to be the best way to prepare and run a GROMACS simulation.

This section is divided as follows: first we will describe how to choose and
prepare your system; we will then discuss the correct way to run your GROMACS 
simulations; and we will finish by discussing some of the post-analysis tools 
that are available.


--------------------
Preparing the system
--------------------

GMX2PDB
=======

The GROMACS ``pdb2gmx`` command is used to convert a coordinate file into a 
set of GROMACS topology files (in these examples, we will assume that the 
file is a ``.pdb`` file, but this is not a necessity). To run this:

.. code-block:: bash

    gmx pdb2gmx -f ${INPUT_FILE}.pdb

where ``${INPUT_FILE}.pdb`` is replaced with the input file name. You will be 
prompted to select the forcefield you would like to use. GROMACS comes with 
a number of AMBER and GROMOS forcefields, as well as a CHARMM and an OPLS-AA
option. You will also need to specify your water model (choices included are 
TIP models, and the SPC and SPC/E models). Specifying the water model here 
results in ``pdb2gmx`` to write a complete topology and will ensure that all
forcefields are consistent if the system needs to be hydrated. The code above 
produces three outputs: a system topology ``topol.top``, an "included 
topology" ``posre.itp``, and a forcefield coordinate file ``conf.gro``. 
Further to these files, ``pdb2gmx`` will output a number of interesting 
details to screen, such as the total mass of the system given the coordinates 
and topology being used as well as the net charge of the system. The charge 
is particularly important to note down and will be used in the `Solvating and 
ionise a system`_ step of system preparation.

More information about the flags and options of this program can be found in 
the GROMACS `PDB2GMX manual
<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html>`_.

.. note::

  There are some limitations to ``pdb2gmx``. This usually gives less accurate 
  results for branched/non-linear bonds.  **There's another one but it slips 
  my mind**

Generating your own forcefield file
-----------------------------------

.. note::

  It is advised to use one of the pre-existing GROMACS forcefields where 
  possible. Only consider generating your own forcefield if the ones 
  available on GROMACS do not fulfill your requirements.

GROMACS comes with a number of forcefields available out of the box. These 
forcefields are stored within the main GROMACS directory in 
``share/gromacs/top``. If the forcefield you want to use is not present, you
will need to generate your own forcefield files. To do this, create a 
directory ``<forcefield>.ff`` in your working directory, with ``<forcefield>``
replaced by a sensible name. Within this directory, create a 
``forcefield.doc`` file and write a simple one-sentence description of your 
forcefield -- this description is what will come up in ``pdb2gmx`` when you 
choose a forcefield. Next, generate a ``forcefield.itp`` included topology 
file. This file is a topology file where you can define the properties of 
atoms, bond, angles, dihedrals, *etc.*. You can find more information about 
generating topology files from scratch in the GROMACS manual `file format page
<http://manual.gromacs.org/documentation/2019.1/reference-manual/file-formats.html#top>`_.

Using a ``<forcefield>.ff`` directory has a number of advantages over writing 
out your system topologies directly. For one, this allows for better 
reproducibility in the even that you want to simulate a new system with this 
forcefield. It also has a number of functionalities that can be useful. For 
instance, adding a ``watermodels.dat`` file into the forcefield directory 
makes it easy to keep track of water models available. A line and descripting 
can be added in this file for each water-model included topology file. This 
file is what prompts the choice of water model in ``pdb2gmx``.

Once it is populated, running ``pdb2gmx`` in the directory containing your ``<forcefield>.ff`` directory will result in your new forcefield being included 
at the top of the list of selectable forcefields. If you are happy with your 
``<forcefield>.ff`` directory and you will use it a lot (and if you have the 
correct permissions to edit parts of the GROMACS directory), you can copy it to 
the ``share/gromacs/top`` of the GROMACS directory (or to ``$GMXLIB`` which 
should be the same directory). In doing so, your forcefield will become a 
permanent part of the forcefields that ``pdb2gmx`` can use.

For more information on generating your own forcefield, please see the GROMACS
manual pages about `adding a residue<http://manual.gromacs.org/documentation/2019.1/how-to/topology.html>`_
and `force field organisations<http://manual.gromacs.org/documentation/2019/reference-manual/topologies/force-field-organization.html>`_.

Preparing and solvating your simulation box
===========================================

Now that a topology has been generated, the next step is to generate a 
simulation box into which to place this topology. For this, use the 
``editconf`` command. This tool has a number of functionalities, including 
generating and orienting a simulation box, and filing it with pre-generated 
topologies. To create a simulation box with ``editconf``, run:

.. code-block:: bash

  gmx editconf -f ${INPUT}.gro -o ${OUTPUT}.gro -center -d ${SEPARATION} -bt ${BOX_TYPE}
  
where ``${INPUT}.gro`` is the input forcefield-compliant coordinate file, 
``${OUTPUT}.gro`` is the chosen output name (the default is ``out.gro``), 
the ``-c`` flag will place the system described in ``${INPUT}.gro`` into the 
centre of the simulation box, ``-d ${SEPARATION}`` defines the minimum 
separation between the input and the edge of the box (units are in nm), and 
``-bt ${BOX_TYPE}`` defines the type of box for the simulation (triclinic is 
the default, but other options are cubic, octohedral, or dodecahedral). There 
are a number of other ``editconf`` options, predominantly to have more 
control over defining the simulation box. These can be found in the GROMACS 
manual `gmx editconf page
<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-editconf.html>`_.

Generating a system from a GROMACS topology
===========================================

Here, I'll talk about ``insert-molecule`` (for generating liquids) and how to generate a bilayer.

Solvating and ionise a system
=============================

The aptly-named ``solvate`` tool can be used to create a box of solvent or 
to solvate a pre-existing box. To use it, run:

.. code-block:: bash

  gmx solvate -cp ${SOLUTE}.gro -cs ${SOLVENT}.gro -o ${OUTPUT}.gro -p ${TOPOLOGY}.top
  
where ``${SOLUTE}.gro`` is the simulation box configured using the steps 
described above, ``${SOLVENT}.gro`` is the solvent configuration file (node 
that GROMACS has a number of pre-defined solvent configuration files but that 
you can also prepare and use your own), and ``${TOPOLOGY}.top`` is the 
topology obtained when running `GMX2PDB`_. If using a GROMACS-provided 
solvent, the addition of this solvent should not alter the net charge of the 
system.

If the net charge of your system is already 0, you do not need to add ions 
to neutralise your system (and can therefore skip this passage). If, on the 
other hand, your system has a non-zero net charge, you may wish to consider 
adding ions to neutralise your system. This is done in two parts: first, you 
need to use the ``grompp`` tool to generate a ``.tpr`` file to which ions can 
be added, and then you must replace some of the recently-added solvent 
molecules with the necessary counterions using ``genion``.

The GROMACS preprocessor tool reads in coordinate and topology files to 
generate an atomic-level input file (with a ``.tpr`` extension). This ``.tpr`` 
file contains all of the parameters needed for all atoms in the system.

For further information, please see the GROMACS manual `gmx solvate
<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html>`_, 
`gmx grompp<http://manual.gromacs.org/current/onlinehelp/gmx-grompp.html>`_, 
and `gmx genion
<http://manual.gromacs.org/documentation/2020.4/onlinehelp/gmx-genion.html>`_ 
pages.

--------------------
Running a simulation
--------------------

----------------------------------
Post-processing and analysis tools
----------------------------------
