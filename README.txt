  ___  ____  __  __
 / __)(  _ \(  \/  )  
 \__ \ )___/ )    (   Statistical Parametric Mapping
 (___/(__)  (_/\/\_)  SPM -  http://www.fil.ion.ucl.ac.uk/spm

                              R E A D M E

________________________________________________________________________

This README gives a brief introduction to the installation and use of
the SPM package. Full details can be found on the SPMweb site:
                 http://www.fil.ion.ucl.ac.uk/spm

A manifest for this release is contained in the file Contents.m
The release is described in the file spm.man
 
________________________________________________________________________
                                                                     SPM

   Statistical Parametric Mapping refers to the construction and
   assessment of spatially extended statistical process used to test
   hypotheses about [neuro]imaging data from SPECT/PET & fMRI. These
   ideas have been instantiated in software that is called SPM.

________________________________________________________________________
                                                            Installation

The SPM software is a suite of MATLAB functions, scripts and data
files, with some externally compiled C routines, implementing
Statistical Parametric Mapping. MATLAB, a commercial engineering
mathematics package, is required to use SPM. MATLAB is produced by The
MathWorks, Inc.  Natick, MA, USA. http://www.mathworks.com/
eMail:info@mathworks.com. SPM requires only core MATLAB to run (no
special toolboxes are required).

SPM99 is written for Matlab 5.2 under UNIX. (SPM will not work with
versions of Matlab 5 prior to 5.2, including Matlab 4.) Binaries of the
external C-mex routines are provided for Solaris, Linux and Windows
only, users of other UNIX platforms need an ANSI C compiler to compile
the supplied C source (Makefile provided: spm_MAKE). See
http://www.fil.ion.ucl.ac.uk/spm/distrib.html for details.

( Whilst the majority of the code is implemented as MatLab functions    )
( & scripts containing standard MatLab commands, a number of features   )
( specific to the UNIX version have been used. SPM also uses external   )
( C programs, linked to MatLab as C-mex files, to perform some of the   )
( more computationally intensive operations.  Some of these latter C    )
( programs use UNIX system calls to implement SPMs "memory mapping",    )
( mapping disk resident image volumes into memory. (Specifically, the   )
( calls are mmap (mman.h) in spm_map_vol.c; munmap (mman.h) in          )
( spm_unmap_vol.c; & readdir (dirent.h) in spm_list_files.c. Not all    )
( UNIX flavours support mman.h.                                         )

With the compiled c-mex files in place, simply prepend the SPM
directory to your MatLab path to complete the installation. (Type `help
path` in matlab for information on the MatLab path.)

________________________________________________________________________
                                                         Getting started
                                                         
SPM is invoked with the command `spm` at the MatLab prompt. We
recommend you start by reviewing the help system, by selecting "About
SPM" from the splash screen. This initially displays the "spm.man"
topic, detailing this release. Press the "Menu" button to display a
representation of the SPM menu window, with buttons linked to
appropriate help pages.

Before attempting to analyze data using SPM, we recommend you spend
some time reading. It is essential to understand the concepts of
Statistical Parametric Mapping in order to effectively use the software
as a research tool. You should begin with the SPMweb pages,
particularly the "Documentation" page. Of the resources listed there
perhaps the most useful starting point are the SPM course notes, which
explain the concepts and theories implemented in SPM at a lower level
than the articles in the peer reviewed literature. There is no manual.

Note that SPM uses Analyze format images as standard, although it can
also read MINC & ECAT-7 images. You will either need to convert your
image files to one of these formats (preferably Analyze), or construct
an additional module for the SPM memory mapping subsystem to read your
file format. Image conversion utilities for your image file format may
be available in other packages, or may have been specially written by
other SPM users. (Consult the SPM email discussion list, described
below, by first searching the archives, and posting a query if
necessary.) Unfortunately we have no resources to provide image
conversion software, although we will collaborate in developing SPM
memory mapping read-modules for popular image formats for inclusion in
SPM.

________________________________________________________________________
                                                               Resources

The SPMweb site is the central repository for SPM resources:
                 http://www.fil.ion.ucl.ac.uk/spm
Introductory material, installation details, documentation, course
details and patches are published on the site.

There is an SPM eMail discussion list, hosted at <spm@mailbase.ac.uk>.
The list is monitored by the authors, and discusses theoretical,
methodological and practical issues of Statistical Parametric Mapping
and SPM. Subscribe by sending the one line message: "join spm firstname
lastname" to <mailbase@mailbase.ac.uk>. (Users at NIH or UC-Davis
should join their local SPM feeds.) The SPMweb site has further
details:
                 http://www.fil.ion.ucl.ac.uk/spm/help

Please report bugs to the authors at <spm-authors@fil.ion.ucl.ac.uk>.
Peculiarities may actually be features, and should be raised on the SPM
eMail discussion list, <spm@mailbase.ac.uk>.

________________________________________________________________________
                                       Disclaimer, copyright & licencing

SPM (being the collection of files given in the manifest in the
Contents.m file) is free but copyright software, distributed under the
terms of the GNU General Public Licence as published by the Free
Software Foundation (either version 2, as given in file spm_LICENCE.man, or
at your option, any later version). Further details on "copyleft" can
be found at http://www.gnu.org/copyleft/. In particular, SPM is
supplied as is.  No formal support or maintenance is provided or
implied.

________________________________________________________________________
SPM is developed by members and collaborators of the
                              Wellcome Department of Cognitive Neurology

@(#)README.txt 2.5 00/01/25
