## Setting up necessary drivers
For now, iOptron mount does not need any special drivers. Communications can  be done via the rotctl command.

To use the AtikSDK python wrapper, the drivers and pythono library must be install. First, install the prereqs for the AtikSDK drivers as detailed in chapter 4 of the SDK Developer's Guide. To use the Python library on Windows (yuck) you will have to run the python prereq .exe installer. Once all prerequisites are installed, navigate to the directory containing the .whl file and run...

pip install -m Atik_Python_SDK-1.5.1-py3-none-any.whl

You should be able to use the AtikSDK python library from here. For example usage, see the scripts in the examples folder.