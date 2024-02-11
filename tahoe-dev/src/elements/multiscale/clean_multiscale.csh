#!/bin/csh
set DIRS = (. \
	continuum \
	materials \
	FEA \
	FEA/multiscale_elements \
	FEA/multiscale_elements/coarsescale \
	FEA/multiscale_elements/finescale \
	FEA/data_processing \
	FEA/shape_functions \
	FEA/integration \
	FEA/matrix)

foreach dir ($DIRS)
		echo $dir
		rm $dir/*.o $dir/*.d*
end
