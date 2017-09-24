Python version of PPS+

This directory contains the python version of PPS+. 

Scripts in directory “python_scripts” were developed by Jessika Fidler. 

These scripts replace the original Ruby source code, therefore this version does not need any Ruby dependency.

In the PPS+ VM (Version 1.4), replace the "tools" folder by the decompressed "tools/tools_py.tar.xz" folder.
(The python scripts are there contained in directory: “tools/PhyloPythiaS/1_3/python_scripts”.)

Modify the configuration file, s.t. it contains the following entries:

# set the maximum number of processors to be used
processors=4
#
# this settings means that the python version of PPS will be used instead of the Ruby version of PPS
ruby=False

See "example_config_files", which contains several example configuration files. 
