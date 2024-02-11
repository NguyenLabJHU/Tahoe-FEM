The images from the synchrotron are in IDL format, which means you need access
to IDL to read them. Lots of the Engineering computers have it; I just ssh'ed
into one and ran it remotely. Assuming you have the 'CARS_tomography' plugin, the following commands might be of use:

export IDL_PATH=CARS_tomography  # For Unix systems, sets the plugin directory

vol=read_tomo_volume('tenA.volume')  # Reads the file named 'tenA.volume'
window,0,xsize=170,ysize=170         # Opens a window to draw on
tvscl,vol[*,*,350]                   # Displays the 350th slice
write_png,'test.png',tvrd()          # Writes whatever is displayed to an image

for i=350,400 do tvscl,vol[600:750,600:750,i]  # Animates several slices


IDL makes grayscale images, but the data seems to be on an arbitrary scale. It
adjusts depending on the overall brightness level, so some slices might look
brighter or darker even when they aren't.


The 'slices' directory I have included has 6 small slices from tenA.volume,
named according to the z-depth.
