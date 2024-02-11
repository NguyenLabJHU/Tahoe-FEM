This is the progress I made on compressing the sub-volume data from the
synchrotron. I had a lot of difficulty getting things to simulate smoothly,
and only made it about halfway. Note that this is a simple ellipsoidal 
approximation we received, and isn't very good.

Each of the subdirectories here contains one step of the simulation I attempted.
The 'settle' one was just to resolve some problems with the particles
overlapping, and didn't do any compression. 'iso-0-1kpa' was a successful
compression of the particles from 0 to 1 kpa. I had difficulties with stability,
so the time-steps are very small and it took ~2 days to finish. I was not able
to get 'iso-0-30kpa' to work right, it just kept blowing up. I would suggest
playing with the scale, as these particles are quite a bit smaller than the ones
Beichuan tested with and his simulations didn't need such small time-steps.

The intended steps are as follows:
isotropic, 0-1 kpa initialization
isotropic, 1-30 kpa to establish a confining stress
triaxial, 30 kpa, the final simulation
