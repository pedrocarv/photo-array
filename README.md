Windows 10
Python 3.6

Drivers:
 - ASCOM Platform 64
 - ASCOM Platform 64 Developer
   - Use "PlatformDeveloperHelp.chm" as ASCOM reference
 - iOptron ASCOM driver
 - ATIK DLL
 - RS232 converter driver (it's got a part number on it)
 
-------------------------------------------------------
What is going on here?

This is the software database for SCIP. It includes code to run the mount, the 
cameras, and the temperature/humidity sensor. It also includes code for
processing the data. 

-------------------------------------------------------
How is the Python script communicating with the hardware?

The cameras are using functions defined in the DLL file. The ctypes library is 
required to access the funtions in the DLL since they are written in C.

The mount is using an ASCOM driver. ASCOM is a layer between the Python and
the driver. See "PlatformDeveloperHelp.chm" for documentation of the 
functionality. Go to ASCOM Namespaces > ASCOM.DriverAccess > Telescope Class.

The temperature sensor is connected to an Arduino, and Python communicates with
it using a COM port. The Arduino code is in SCIP_operation > si7021.

-------------------------------------------------------
Folders:

scip_operation: All the operation code for running SCIP autonomously or 
otherwise. Anything that communicates with hardware is here. 

data_processing: Code for processing images and temperature/humidity data

snr_simulation: Code for simulating the SNR of SCIP

ref_stars: reference star information

find_sunrise_sunset_shadowheight: Scripts for finding when sunrises and sunsets
occur at a given location, and for plotting the shadowheights at a given 
location. Can also plot when the local shadow height is above 500 km (no local
photoelectrons) AND the geomagnetic conjugate shadowheight is below 500 km 
(photoelectrons will come from the geo conj to the local). Useful for planning
when to observe.

8446_storms: a script that looks at what data is available from Arecibo from
years of previous observations. Could be useful if you want to analyze some 
older data. Talk to Prof. Waldrop.

radiative_transfer: the radiative transfer (RT) model (Fortran, called from
Python) and the exospheric H density retrieval. See the next section to build
and run it in parallel.

-------------------------------------------------------
Radiative transfer: building and running in parallel

The RT forward model (`forward`) and line-of-sight model (`los_1x`) are Fortran
modules compiled with f2py. They are multithreaded with OpenMP and work on
macOS and Linux.

1) Build the modules (once, and again after any change to the Fortran code):

```bash
conda env create -f photom.environment.yaml   # first time only
conda activate photom
cd radiative_transfer/rt_ft
./build_1x.sh
```

`build_1x.sh` builds both modules multithreaded and places them in
`radiative_transfer/`. The other `build_*.sh` scripts build the forward model
multithreaded but their LOS model single-threaded (see
`radiative_transfer/rt_ft/README.md`). Without conda (Debian/Ubuntu), install
`gfortran liblapack-dev python3-dev` with apt and `numpy meson ninja` with pip,
then run `./build_1x.sh`.

2) Run as usual. By default every forward and LOS call uses all CPU cores; no
code changes are needed. Set `OMP_NUM_THREADS` to use fewer threads:

```bash
python osse_test.py                    # all cores
OMP_NUM_THREADS=4 python osse_test.py  # 4 threads
```

- `rt_inversion.inversion_init(..., cores=N)` inverts the wedges in N worker
  processes and gives each worker cores/N threads automatically (unless
  `OMP_NUM_THREADS` is already set). Its workers use the `spawn` start method,
  so scripts that call it need an `if __name__ == "__main__":` guard.
- If you start your own worker processes, set `OMP_NUM_THREADS` to roughly
  cores / number of workers, and use the `spawn` start method on Linux: GNU
  OpenMP hangs in workers created with the default `fork`.
- The `ncore` argument of `rt_common.run_forward_los` / `run_los` has no
  effect; the parallelism comes from the OpenMP threads.

On a 12-core M2 Max, one forward model takes about 0.9 s (34 s with the
original single-threaded code) and the retrieval in `osse_test.py` takes about
15 s (about 10 min before).

-------------------------------------------------------
In the "scip_operation" folder:

Note that data and log files are saved to:
C:/Users/photometer/Box Sync/SCIP/Synced_Raw_Data/

"main_run_autonomously.py" starts a night of observation. It will give you 
errors if it cannot connect to the mount, cameras, or temp sensor. You MUST 
align the mount and set it to the zero position first. This script checks 
the schedule.txt file to find out if it is time to observe and what the 
observing settings should be. You need a valid schedule file, or it won't work.
Add new obsering sessions to schedule.txt to schedule the sessions. Follow the
format in schedule.txt. 

"observation_functions.py" contains all the functions for deciding what actions
to take during an autonomous operation. It finds the sun status (daylight, 
twilight, or night) and observes the zenith airglow or the calibration stars as 
appropriate. 

"photometer_classes.py" has the classes for the camera, mount, and temp sensor.
You should think of these as a layer of abstraction between the rest of the 
Python control code and the access to the hardware. The ATIKCamera class is 
basically a wrapper for the DLL functions which are accessed through ctypes. 
The Binocular class makes two ATIKCamera objects (one for each camera) and 
commands them both at once. The CubeProMount class has functions which access 
the ASCOM commands for controlling the mount. TempHumidSensor uses the Python 
package "serial" to talk to the Arduino which is connected to the sensor. 
Any time you want the hardware to do something, you should do it through the 
appropriate object. If you need functionality that isn't there, add new 
functions to the appropriate class. For example, if you want both cameras to
take a picture at the same time, use the function in the Binocular class 
called "dual_exposure".

"util.py" has some extra functions like logging and saving data in FITs files.
The "log" file logs most everything that happens as the code runs. It should 
be helpful for diagnosing wierd behavior.

-------------------------------------------------------
Operation process:

Read the instructions on how to align the mount! Found in the Box folder at
SCIP > SCIP_documents. This process takes practice. It is not easy.

Setup to do during the day:
- Check that dessicant is still sufficiently dry
- Check that power and internet is working
- Check that heat and/or AC is working
- Level the mount head
- Attach the power block and telescopes
- Test computer connections
- Check that Box Sync is working
- Align the scopes in roll and pitch, and focus the telescopes
    - Take pictures of a very distant object to do this
    - The off-band does not focus very well
- Set up Teamviewer or your preferred system
- Turn on heat/AC as needed
        
Nighttime setup before observation session starts:
- Disconnect mount from computer
- Power on the mount
- Follow mount alignment instructions
    - Set the Zero Position
    - One or two star align
    - Return to Zero positon 
- Connect mount to computer
- Connect and power cameras, connect temp sensor
- Turn on heat or AC
- Turn on Teamviewer or your preferred system
- Start session

Daily during observation session:
- Remotely check temperature and humidity often
- Cover the dome with the sun shield and tie it down at sunrise
- Remove the sun sheild at sunset

Storage:
- AC and dome can stay indoors at field site
- Anything else with any kind of electronics should be in a controlled
environment (like ECEB). Filters will be damaged by temperature cycling.

-------------------------------------------------------
This codebase is a work in progress.
You can always reach out to dawnhaken@gmail.com if you are stuck.

 

