FaceT::FaceT
(SurfaceT& surface,iArrayT& connectivity, dArrayT& coordinates);
{
        fNumNodes = Length.connectivity;
        fSurface = surface ;
        fConnectivity.Allocate(fNumNodes);
        for (i = 0; i<fNumNodes ; i++) {
                fConnectivity[i] = connectivity[i];
        }
        fCoordinates = coordinates;

}

