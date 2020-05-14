# Building 3-D Workshop Tutorial

## Notes

This was the example case used throughout the May 13th, 2020 blastFoam workshop (thank you to everyone who attended!). This case is designed to show a simple 3-D case setup using an STL file with snappyHexMesh and adaptive mesh refinement. The a c-4 charge placed inside a L-Shaped building. The c-4 is activated using the linear activation model.

The case takes approximately 5 min to run on a four-core desktop. Some additional notes:

- There are a few examples of how to add probes and/or surfaces and other runtime processing directives to an OpenFOAM/blastFoam case. An example ```probesAndSurfaces``` file is included in the system directory. One can cut/paste/modify these examples into the functions{}; section of the ```controlDict``` file. An alternative, and handy, way to acheive this is to use the ```foamGet``` command to generate a probes (or other) file and include the functions{}; section of the ```controlDict```.

Example: 

I want to generate a basic probe (gage) that I can then modify get pressure/time histories at discrete locations in my calculation. 

```bash
cd <root-case-directory>
foamGet probes # answer (1) to get the probes file and (2) to get a configuration file to set output times etc.
``` 
You will now have a file called "probes" in the system directory. Edit that file so that the probes are in the desired locations. Then, in the controlDict file, find the functions{}; section of the case. Add this line(s):

```c++
functions
{
	#includeFunc "probes"; // this will load the system/probes file that was created with the $ foamGet probes command
};
```

Alternatively, you can do this manually directly in the ```controlDict``` file:


```c++
functions
{

    pressureProbes
    {
        type            probes;
        functionObjectLibs ("libsampling.so");
        writeControl    timeStep;
        writeInterval   1;
        probeLocations
        (
			// probe locations in ( x y z ) format, unlimited number -- needs to be inside the domain!
			( 1.999 0.001 0.5)
			// ( 0 0 0 ) // etc... 
        );
        fields
        (
            p // fields you want at each probe location: p U T rho etc. 
        );
    }

};
```

During/after running the simulation, results will be in the ```postProcessing``` directory.


