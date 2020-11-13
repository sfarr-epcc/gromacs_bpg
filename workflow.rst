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
the GROMACS 
`PDB2GMX manual<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html>`_.

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
generating topology files from scratch in the GROMACS manual 
`file format page<http://manual.gromacs.org/documentation/2020.4/reference-manual/file-formats.html#top>`_.

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
manual pages about 
`adding a residue<http://manual.gromacs.org/documentation/2020.4/how-to/topology.html>`_
and `force field organisations<http://manual.gromacs.org/documentation/2020/reference-manual/topologies/force-field-organization.html>`_.

Preparing and solvating your simulation box
===========================================

Now that a topology has been generated, the next step is to generate a 
simulation box into which to place this topology. For this, use the 
``editconf`` command. This tool has a number of functionalities, including 
generating and orienting a simulation box, and filing it with pre-generated 
topologies. To create a simulation box with ``editconf``, run:

.. code-block:: bash

  gmx editconf -f ${INPUT}.gro -center -d ${SEPARATION} -bt ${BOX_TYPE} \
               -o ${OUTPUT}.gro
  
where ``${INPUT}.gro`` is the input forcefield-compliant coordinate file, 
``${OUTPUT}.gro`` is the chosen output name (the default is ``out.gro``), 
the ``-c`` flag will place the system described in ``${INPUT}.gro`` into the 
centre of the simulation box, ``-d ${SEPARATION}`` defines the minimum 
separation between the input and the edge of the box (units are in nm), and 
``-bt ${BOX_TYPE}`` defines the type of box for the simulation (triclinic is 
the default, but other options are cubic, octohedral, or dodecahedral). There 
are a number of other ``editconf`` options, predominantly to have more 
control over defining the simulation box. These can be found in the GROMACS 
manual 
`gmx editconf page<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-editconf.html>`_.

.. Generating a system from a GROMACS topology
.. ===========================================
.. 
.. Here, I'll talk about ``insert-molecule`` (for generating liquids) and how to generate a bilayer.

Solvating a system
===============================

The aptly-named ``solvate`` tool can be used to create a box of solvent or 
to solvate a pre-existing box. To use it, run:

.. code-block:: bash

  gmx solvate -cp ${SOLUTE}.gro -cs ${SOLVENT}.gro -p ${TOPOLOGY}.top \
              -o ${OUTPUT}.gro
  
where ``${SOLUTE}.gro`` is the simulation box configured using the steps 
described above, ``${SOLVENT}.gro`` is the solvent configuration file (node 
that GROMACS has a number of pre-defined solvent configuration files but that 
you can also prepare and use your own), and ``${TOPOLOGY}.top`` is the 
topology obtained when running `GMX2PDB`_. If using a GROMACS-provided 
solvent, the addition of this solvent should not alter the net charge of the 
system.

For further information, please see the GROMACS manual 
`gmx solvate<http://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html>`_

Creating a charge-neutral system
================================

If the net charge of your system is already 0, you do not need to add ions 
to neutralise your system (and can therefore skip this passage). If, on the 
other hand, your system has a non-zero net charge, you may wish to consider 
adding ions to neutralise your system. This is done in two parts: first, you 
need to use the ``grompp`` tool to generate a ``.tpr`` file to be used when 
adding ions, and then you must replace some of the recently-added solvent 
molecules with the necessary counterions using ``genion``.

The GROMACS preprocessor tool ``grompp`` reads in coordinate and topology 
files to generate an atomic-level input file (with a ``.tpr`` extension). 
This ``.tpr`` file contains all of the parameters needed for all atoms in 
the system. We will go into more details about the ``grompp`` tool in the 
`Running a simulation`_ section. For now, the important part is that, to 
generate a run input ``.tpr`` file, ``grompp`` needs a structure (``.gro``) 
file, a topology (``.top``) file, and a file defining the instructions for 
the simulation run (this is kept in an ``.mdp`` file). This ``.mdp`` file can 
be kept empty when ionising the system as no actual simulation is to be run. 
To generate the ``,tpr`` file, run:

.. code-block:: bash

  gmx grompp -f ${RUN_FILE}.mdp -c ${COORDINATES}.gro -p ${TOPOLOGY}.top \
             -o ${OUTPUT}.tpr
  
It is likely that ``grompp`` will output a number of notes to screen (one of 
which should be reminding you of the net non-zero charge of your system). In 
this case, these can be ignored (this is an exception and is not usually true).

Now that the ``.tpr`` has been generated, ``genion`` can be used to make the 
charge of the system neutral. The system charge is decreased by replacing a 
number of parts of the system with anions and cations. This is done by 
running the following (note that the ``${INPUT}.tpr`` named below is likely 
to be the ``${OUTPUT.tpr}`` generated in the ``grompp`` step above): 

.. code-block:: bash

  gmx genion -s {INPUT}.tpr -p ${TOPOLOGY}.top -neutral -o ${OUTPUT}.gro
             
You will be prompted to choose the group within your system (solvents, 
solutes, protein backbones, *etc.*) that you would like ions to replace, with 
the frequency of occurrence of each group also shown. Note that some groups 
may have overlap completely and be different names for the same group. In 
general, it is best to replace solvent molecules with ions (the group named 
``SOL``). Once a group is chosen, ``genion`` will replace a number of that 
group with anions and cations until the system is charge neutral. The default 
anion name is ``CL``, though this name can be changed with the ``-nname`` 
flag, and the default cation name is ``NA`, but this name can be changed with 
the ``nname`` flag. By default, the cation and anion charges are 1 and -1 
respectively, but this can be changed with the ``-pq`` flag for the cation and 
the ``-nq`` flag for the anion.

For further information, please see the GROMACS manual  
`gmx grompp<http://manual.gromacs.org/current/onlinehelp/gmx-grompp.html>`_, 
and `gmx genion<http://manual.gromacs.org/documentation/2020.4/onlinehelp/gmx-genion.html>`_ 
pages.

--------------------
Running a simulation
--------------------

This section describes how to set the GROMACS simulation parameters, how to 
generate a run input file from a GROMACS topology and parameter file, how to 
run a simulation in GROMACS, and how to analyse the results produced. It is 
assumed that you already have a system topology ready to use (by following 
the steps in the `Preparing the system`_ section) -- if this is not the case, 
and if you are unsure how to create this topology, please read through that 
section.

Creating a run parameter file
=============================


A GROMACS molecular dynamics parameter (``.mdp``) file defines the simulation 
parameters to be used during a simulation. A number of options can be set in 
this script, including: defining the simulation integrator that will define 
the method used to solve Newton's equations to propagate the system forward in 
time; setting the size of the simulation timestep and total simulation time; 
setting the restrictions within which the system will be simulated (such as 
setting a system pressure/temperature through a thermostat or barostat); 
setting or adjusting the way the simulation forcefield is interpreted (by 
*e.g.* defining the way short- and long-ranged interactions are calculated 
and at what distance they are truncated); to define which simulation 
properties to output (and the output frequency); and many more options. Given 
the number of options and variables that can be included, not included, or 
kept as default, we will not go over all of the options here and will instead 
look at and explain an example molecular dynamics parameter file. You can find 
a list of all available options in the GROMACS manual
`molecular dynamics parameters page<http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html>`_.

Example molecular dynamics parameter file
-----------------------------------------

The GROMACS manual has the following 
`example script<http://manual.gromacs.org/documentation/2018/user-guide/file-formats.html#mdp>`_:

.. code-block:: bash

  ; Intergrator, timestep, and total run time
  integrator               = md
  dt                       = 0.002
  nsteps                   = 500000
  
  ; Logs and outputs
  nstlog                   = 5000
  nstenergy                = 5000
  
  ; Bond constraints
  constraints              = all-bonds
  constraint-algorithm     = lincs
  
  ; Van der Waals interactions
  vdwtype                  = Cut-off
  rvdw                     = 1.0
  cutoff-scheme            = Verlet
  DispCorr                 = EnerPres
  
  ; Coulombic interactions
  coulombtype              = PME
  rcoulomb                 = 1.0
  
  ; Thermostat
  tcoupl                   = V-rescale
  tc-grps                  = Protein  SOL
  ref-t                    = 300      300
  tau-t                    = 0.1      0.1
  
  ; Barostat
  pcoupl                   = Parrinello-Rahman
  ref-p                    = 1.0
  tau-p                    = 2.0
  compressibility          = 4.5e-5

First note that, while the the example above is ordered in a sensible way, 
with commands grouped by what they are defining (*e.g.* temperature, pressure, 
van der Waals interactions, *etc.*), the order in which the individual 
commands are written should not matter. Having said that, we would recommend 
grouping commands affecting similar simulation aspects together to help 
future readability. Also, if the same command appears twice in a 
``.mdp`` file, the second appearance will override the first.

The first block of the example script defines the molecular dynamics 
integrator as a Verlet leap-frog algorithm (``integrator = md``), declares 
that the simulation timestep will be 2 fs (``dt = 0.002``, where the default 
unit is ps), and that the simulation will run for a total of 500,000 *dt*
timesteps (``nstep = 500000``) or 1 ns.

The next block defines the simulation outputs. ``nstlog`` sets the time 
interval between each output to log (``md.log``) of the energy components and 
physical properties of the system at 5,000 *dt*. ``nstenergy`` sets the time 
interval between each output to the energy file (``ener.edr``) of the energy 
components of the system at 5,000 *dt* -- note that this file is written in 
binary.

In this example, all bonds are constrained and set to be rigid. This is done 
with the ``constraints = all-bonds`` command. Furthermore, the constraint 
algorithm is set to be the linear constraint solver algorithm with the 
``constraint-algortihm = LINCS`` command.

The van der Waals interactions are set as truncated (``vdwtype = cutoff``), 
with a cutoff distance of 1 nm (``rvdw = 1.0``). This means that no van der 
Waals interactions will be computed for pairs of particles whose 
centre-of-mass separation greater than 1 nm. To save in simulation time, a 
neighbour-list cutoff scheme is used. The ``cutoff-scheme = Verlet`` command 
specifies how this list is generated. A long-ranged dispersion correction to 
the energy and pressure is considered here with the ``DispCorr = EnerPress`` 
command. The Coulombic interactions will be calculated using the smooth 
particle-mesh Ewald (SPME) method (``coulombtype = PME``), with an interaction 
cutoff of 1 nm (``rcoulomb = 1.0``). 

The thermostat used for this simulation is defined by the 
``tcoupl = v-rescale`` -- in this case, the velocity rescaling algorithm is 
used. The ``tc-grps`` is there to specify that the protein and solvent 
(``SOL``) should have separate heat baths for this simulations. The reference 
temperature (or desired temperature) is set by ``ref-t``. In this case, the 
reference temperature for both the protein and the solvent have been set to 
300 K. Note that the reference temperature must be set for every group defined 
in ``tc-grps`` and that these temperatures do not need to be the same.
The temperature coupling time constant is defined by ``tau-t`` as being 0.1 
ps. This coupling time constant controls how much, and over what timescales, 
the thermostat is allowed to fluctuate -- a high coupling time usually leads 
to long, low-amplitude oscillations about the mean the system temperature, 
whereas low coupling  often leads to short, large-amplitude oscillations. 
Note that, like with the system temperature, this must be defined for each 
group separately (hence the two values inthe example script).

In this example script, the barostat is defined with the ``pcoul`` parameter 
as the Parinello-Rahman barostat. The reference (or desired) pressure is set 
at 1 atm with the ``ref-p`` command, and the coupling time constant ``tau-p`` 
is set to 2 ps. Much like the temperature coupling time constant set for the 
thermostat, the pressure coupling time constant is used to dictate the 
frequency and amplitude of fluctuations during a simulation. Finally, the 
``compressibility`` parameter is used to define the compressibility of the 
system (how the volume of the system changes as pressure is changed). In this 
case, it is set as 4.5e-5 bar^-1.

Generating your simulation input file
-------------------------------------

Once you have prepared your ``.mdp`` file, you are ready to combine it with 
the topology you've prepared to create a run input ``.tpr`` file. For this, we 
will use the GROMACS pre-processing tool ``grompp``. This is very similar to 
the step described in the `Creating a charge-neutral system`_ section, but 
with more care regarding the warnings that are output. Like before, this is 
done by running:

.. code-block:: bash

  gmx grompp -f ${RUN_FILE}.mdp -c ${COORDINATES}.gro -p ${TOPOLOGY}.top \
             -o ${OUTPUT}.tpr
             
where ``${RUN_FILE}.mdp`` is file discussed in 
`Creating a run parameter file`_, and ``${COORDINATES}.gro`` and 
``${TOPOLOGY}.top`` were generated following the instructions in the 
`Preparing the system`_ section. The ``${OUTPUT}.tpr`` file generated here 
is the only file needed to proceed with running a GROMACS molecular dynamics 
simulation.

Running your simulation
-----------------------

With the run input ``.tpr`` file now generated, we are ready to run a GROMACS 
simulation. For this, we will use the ``mdrun`` command:

.. code-block:: bash

  gmx mdrun -s ${INPUT}.tpr
  
This command will run the simulation with the topology that you've prepared 
and the molecular dynamics parameters that you've chosen.

It is possible to add flags to ``mdrun`` to alter some of the parameters that 
had been set in the molecular dynamics parameter file. For instance, the 
``-nsteps`` flag can be used to override the number of timesteps that the 
simulation should run for. Also, there are a number of useful options for 
defining input files (and input file types), output files, and parameters 
related to the computational system on which you are running (such as the 
``-nt`` option to set the number of MPI threads that the simulation should 
use). More information on these and other options can be found on the GROMACS 
`gmx mdrun<http://manual.gromacs.org/documentation/5.1/onlinehelp/gmx-mdrun.html>`_
page.

.. Good simulation practices
.. =========================
.. 
.. Energy minimisation
.. -------------------
.. 
.. System equilibration
.. --------------------
.. 
.. Production runs
.. ---------------
.. 

.. Post-processing and analysis tools
.. ==================================

.. With the simulation complete, we can analyse the simulation trajectory and 
.. understand what the simulation has demonstrated. GROMACS offers a number of 
.. post-simulation analysis tools. In this section, we will discuss tools that 
.. can be used to: generate the thermodynamic properties of interest; obtain 
.. radial distribution functions and correlation functions; 
.. 
.. Thermodynamic properties of the system
.. --------------------------------------
.. 
.. The GROMACS ``energy`` tool can be used to extract energy components from an 
.. energy (``.edr``) file. By default, this tool will generate an XMGrace file. 
.. To use this, run:
.. 
.. .. code-block:: bash
.. 
  .. gmx energy -f ${INPUT_ENERGY_FILE}.edr -o ${OUTPUT_XMGRACE_FILE}.xvg
..   
.. When running this, you will get a prompt asking which property you would like 
.. output (*e.g.* potential energy, kinetic energy, pressure, temperature, 
.. *etc.*). Enter the correct number to generate an XMGrace file that, when 
.. plotted, will show you how that property varied over the simulation run. 
.. There are a number of other options for the ``energy`` command, and these 
.. can be found in the GROMACS manual 
.. `gmx energy<http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-energy.html#gmx-.. energy>`_
.. page.
.. 
.. Radial distribution function and correlation functions
.. ------------------------------------------------------
.. 
