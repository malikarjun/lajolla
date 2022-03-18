# CSE 272 - Final Project
UCSD CSE 272 Final Project - Photon Mapping and Thin-Film Interference.

Author: Ziyang Zhang and Mallikarjun Swamy

# Build
All the dependencies are included. Use CMake to build.
If you are on Unix systems, try
```
mkdir build
cd build
cmake ..
```

# Run the test scene
```
cd build
./lajolla ../scenes/\photonmap_test/pm_pool.xml
```

Note: Since the pool scene might take long time to render, the `../scenes/photonmap_test/pm_cbox_basic.xml`
can be rendered to sanity test the implementation.
