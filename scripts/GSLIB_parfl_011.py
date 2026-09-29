# -*- coding: utf-8 -*-
"""
Created on Tue May  9 17:39:29 2023

@author: Witty
"""
def GSLIB_GAM(gamexe,gampar,data,outfl,GSLIB_GAM_params):
    (nvar,ivar1,tmin,tmax,nsim,grid_params,ndir,nlag,
     xdir,standardize,nvarg,ivtail,ivhead,ivtype) = GSLIB_GAM_params.values()
    (nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn) = grid_params.values()
    """
    GSLIB Variogram Program for the gridded data (regularly spaced). 
    
    Parameters
    ----------
    gamexe : string
        GSLIB executable for creation of a correlogram considering regularly 
        spaced data
    gampar : string
        parameter file for gam exe.
    data : string
        the input data in a simplified Geo-EAS formatted file. The data are 
        ordered rowwise (X cycles fastest, then Y, then Z).
    nvar : int
        the number of variables.
    ivar1 : int
        column number of the variables in the data file.
    tmin : float
        lower trimmimg limit - all values less than this value are ignored.
    tmax : float
        upper trimmimg limit - all values greater than this value are ignored.
    outfl : string
        the output variograms are written to a single output file named outfl. 
        The output file contains the variograms ordered by direction and then 
        variogram type specifi ed in the parameter file (the directions cycle 
        fastest then the variogram number). For each variogram there is a one-
        line description and then nlag lines each with the following:
            1. lag number (increasing from 1 to nlag).
            2. average separation distance for the lag. 
            3. the semivariogram value (whatever type was specified).
            4. number of pairs for the lag. 
            5. mean of the data contributing to the tail.
            6. mean of the data contributing to the head.
            7. the tail and head variances (for the correlogram).
            
            The vargplt program documented in Section VI.1.8 of Deutsch and Journel 
            may be used to create PostScript displays of multiple variograms.
      
    isim : int
        the grid or realization number. Recall that realizations or grids
        are written one after another; therefore, if igrid=2 the input file must 
        contain at least 2 nx*ny*nz values and the second set of nx*ny*nz values 
        will be taken as the second grid.
    nx : int
        definition of the grid system (x-axis)
        number of cells in x direction.
    xmn : float
        definition of the grid system (x-axis)
        coordinates at the center of the first blocks in x direction.
    xsiz : float
        definition of the grid system (x-axis) 
        constant block size in x direction.
    ny : int
        definition of the grid system (y-axis)
        number of cells in y direction.
    ymn : float
        definition of the grid system (y-axis)
        coordinates at the center of the first blocks in y direction.
    ysiz : float
        definition of the grid system (y-axis) 
        constant block size in y direction.
    nz : int
        definition of the grid system (z-axis)
        number of cells in z direction.
    zmn : float
        definition of the grid system (z-axis)
        coordinates at the center of the first blocks in z direction.
    zsiz : float
        definition of the grid system (z-axis) 
        constant block size in z direction.
    ndir : int
        the number of directions to consider. All directions are considered for 
        all of the nvarg variograms specified below
    nlag : int
        the number of lags to consider. The same number of lags are considered 
        for all directions
    ixd : int
        the node offset in x coordinate direction as used to define a direction.
    iyd : int
        the node offset in y coordinate direction as used to define a direction.
    izd : int
        the node offset in z coordinate direction as used to define a direction.
    standardize : int
        if set to 1, the semivariogram values will be divided by the variance.
        If not used set to 0.
    nvarg : int
        the number of variograms to compute.
    ivtail : int
        specificaton of variables which are to be used for the tail.
    ivhead : int
        specificaton of variables which are to be used for the tail. For direct
        variograms the ivtail array is identical to the ivhead array. Cross 
        variograms are computed by having the tail variable different from the 
        head variable.
    ivtype : int
        The ivtype variable corresponds to the integer code in the list given 
        in Section III.1 of Deutsch and Journel.
        1 = Semivariogramm; -> predef. indicator variable
        2 = Cross semivariogram
        3 = Covariance
        4 = Correlogram
        5 = General relative semivariogramm
        6 = Pairwise relative semivariogram
        7 = Semivariogram of logarithms
        8 = Semimadogram
        
        9 = Indicator Semivariogram(continuous variable)
        10 = Indicator Semivariogram(categorical variable)
            cut: whenever the ivtype is set to 9 or 10, i.e., asking for an 
            indicator variogram, then a cutoff must be specified immediately 
            after the ivtype parameter on the same line in the input file 
            (see Figure III.3, Deutsch and Journel). Note that if an indicator 
            variogram is being computed, then the cutoff/category applies to 
            variable ivtail(i) in the input file [although the ivhead(i) variable 
            is not used, it must be present in the file to maintain consistency
            with the other variogram types].     
        
    Returns
    -------
    None.

    """    
#  Variogram of regularly Spaced 3-D Data, wrapper for GAM.exe from GSLIB (.exe must be in working directory)



    file = open(gampar, "w")
    file.write("                      Parameters for GAM                                  \n")
    file.write("                      *******************                                  \n")
    file.write("                                                                           \n")
    file.write("START OF PARAMETERS:                                                       \n")
    file.write(str(data) + "                     - file with data                          \n")
    file.write(str(nvar) + " " + str(ivar1) + "  - number of varables,column numbers       \n")
    file.write(str(tmin) + " " + str(tmax) + "   - trimming limits                         \n")
    file.write(str(outfl) + "                    - file for variogram output               \n")
    file.write(str(nsim) + "                     - the grid or realisation number          \n")
    file.write(str(nx) + " " + str(xmn) + " " +
               str(xsiz) + "                     - nx, xmn, xsiz                           \n")
    file.write(str(ny) + " " + str(ymn) + " " +
               str(ysiz) + "                     - ny, ymn, ysiz                           \n")
    file.write(str(nz) + " " + str(zmn) + " " +
               str(zsiz) + "                     - nz, zmn, zsiz                           \n")
    file.write(str(ndir) + " " + str(nlag) + "   - number of directions, number of lags    \n")
    for i in range(ndir):
        ixd,iyd,izd = xdir[i]
        file.write(str(ixd) + " " + str(iyd) + " " +
                   str(izd) + "                      - ixd(1), iyd(1), izd(1)                 \n")
    file.write(str(standardize) + "              - standardize sills? (0=no, 1=yes)        \n")
    file.write(str(nvarg) + "                    - number of variograms                    \n")
    file.write(str(ivtail) + " " +
               str(ivhead) + " " +
               str(ivtype) + "                   - tail var., head var., variogram type    \n")
    file.close()


#%%
# ----------------------------------------------------------------------------
# GAMV
# ----------------------------------------------------------------------------
#%%
def GSLIB_GAMV(gamvexe,gamvpar,data,outfl,GSLIB_GAMV_params):
    
    (icolx,icoly,icolz,nvar,ivar1,tmin,tmax,nlag,xlag,lagtol,
     ndir,xdir,standardize,nvarg,
     ivtail,ivhead,ivtype) = GSLIB_GAMV_params.values()
    
    """GSLIB variogram program for irregularly spaced data
    
        
    Parameters
    ----------
    gamvexe : string
        GSLIB executable for creation of a variogram considering irregularly 
        spaced data.
    gamvpar : string
        parameter file for GAMV.
    data : string
        input data in a simplied Geo-EAS formatted file.
    outfl : string
        the output variograms are written to a single output file named outfl.
        The output file contains the variograms ordered by direction and then 
        variogram type specified in the parameter file (the directions cycle 
        fastest then the variogram number). For each variogram there is a one-
        line description and then nlag lines each with the following.
            1. lag number (increasing from 1 to nlag).
            2. average separation distance for the lag.
            3. the semivariogram value (whatever type was specified).
            4. number of pairs for the lag.
            5. mean of the data contributing to the tail.
            6. mean of the data contributing to the head.
            7. the tail and head variances (for the correlogram).
    icolx : int
        the columns for the x coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data, e.g.,
        2D data may be handled by setting icolz to 0).
    icoly : int
        the columns for the y coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data - see
        iclox.
    icolz : int
        the columns for the y coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data - see
        iclox.
    nvar : int
        the number of variables.
    ivar1 : int
        variable column order in the data file.
    tmin : float
        lower trimmimg limit - all values less than this value are ignored.
    tmax : float
        upper trimmimg limit - all values greater than this value are ignored.
    nlag : int
        the number of lags to compute (same for all directions).
    xlag : int
        the unit lag separation distance.
    lagtol : float
        the lag tolerance. This could be one-half of xlag or smaller to allow 
        for data on a pseudoregular grid. If lagtol is entered as negative or 
        zero, it will be reset to xlag/2. A pair will report to multiple lags 
        if lagtol is greater than one-half of xlag.
    ndir : int
        the number of directions to consider. All these directions are 
        considered for all the nvarg variograms specified below.
    azm : float
        the azimuth angle. 
        e.g.    azm=0 is north
                azm=90 is east
                azm=135 is southeast
    atol : float
        the half-window azimuth tolerance. It is restricted once the deviation 
        from the direction vector exceeds the bandwidth
    bandwh : float
        the azimuth bandwidth. It is the horizontal bandwidth or maximum 
        acceptable horizontal deviation from the direction vector
    dip : float
        The dip angle is measured in negative degrees down from horizontal.
        i.e. dip=0 is horizontal, dip=-90 is vertical downward
    dtol : float
        angular tolerance from dip.
    bandwd : float
        bandwd is the vertical \bandwidth" or maximum acceptable deviation 
        perpendicular to the dip direction in the vertical plane.
    standardize : int
        if set to 1, the semivariogram values will be divided by the variance.
    nvarg : int
        the number of variograms to compute.
    ivtail : int
        specificaton of variables which are to be used for the tail.
    ivhead : int
        specificaton of variables which are to be used for the tail. For direct
        variograms the ivtail array is identical to the ivhead array. Cross 
        variograms are computed by having the tail variable different from the 
        head variable.
    ivtype : int
        The ivtype variable corresponds to the integer code in the list given 
        in Section III.1 of Deutsch and Journel.
        1 = Semivariogramm; -> predef. indicator variable
        2 = Cross semivariogram
        3 = Covariance
        4 = Correlogram
        5 = General relative semivariogramm
        6 = Pairwise relative semivariogram
        7 = Semivariogram of logarithms
        8 = Semimadogram
        
        9 = Indicator Semivariogram(continuous variable)
        10 = Indicator Semivariogram(categorical variable)
            cut: whenever the ivtype is set to 9 or 10, i.e., asking for an 
            indicator variogram, then a cutoff must be specified immediately 
            after the ivtype parameter on the same line in the input file 
            (see Figure III.3, Deutsch and Journel). Note that if an indicator 
            variogram is being computed, then the cutoff/category applies to 
            variable ivtail(i) in the input file [although the ivhead(i) variable 
            is not used, it must be present in the file to maintain consistency
            with the other variogram types].
            
    Returns
    -------
    None.

    """
#  Variogram of Irregularly Spaced 3-D Data, wrapper for GAMV.exe from GSLIB (.exe must be in working directory)



    file = open(gamvpar, "w")
    file.write("                      Parameters for GAMV                                   \n")
    file.write("                      *******************                                   \n")
    file.write("                                                                            \n")
    file.write("START OF PARAMETERS:                                                        \n")
    file.write(str(data) + "                     - file with data                           \n")
    file.write(str(icolx) + " " + str(icoly) +
        " " + str(icolz) + " - columns for X, Y, Z coordinates                              \n")
    file.write(str(nvar) + " " + str(ivar1) + " - number of varables,column numbers         \n")
    file.write(str(tmin) + " " + str(tmax) + "  - trimming limits                           \n")
    file.write(str(outfl) + "                    - file for variogram output                \n")
    file.write(str(nlag) + "                     - number of lags                           \n")
    file.write(str(xlag) + "                - lag separation distance                       \n")
    file.write(str(lagtol) + "                    - lag tolerance                           \n")
    file.write(str(ndir) + "                      - number of directions                    \n")
    for i in range(ndir):
        azm,atol,bandwh,dip,dtol,bandwd = xdir[i]
        file.write(str(azm) + " " + str(atol) + " " +
                   str(bandwh) + " " + str(dip) + " " +
                   str(dtol) + " " + str(bandwd) + "   - azm,atol,bandwh,dip,dtol,bandwd    \n")
    file.write(str(standardize) + "                   - standardize sills? (0=no, 1=yes)    \n")
    file.write(str(nvarg) + "                      - number of variograms                   \n")
    file.write(str(ivtail) + " " +
                   str(ivhead) + " " +
                   str(ivtype) + " - tail var., head var., variogram type                   \n")

    file.close()



 

    # lag         # 1. lag number (increasing from 1 to nlag).
    # avsepdist   # 2. average separation distance for the lag.
    # svargval    # 3. the semivariogram value (whatever type was specied).
    # npairs      # 4. number of pairs for the lag.
    # meantail    # 5. mean of the data contributing to the tail.
    # meanhead    # 6. mean of the data contributing to the head.
    # tailheadvar # 7. the tail and head variances (for the correlogram).

#%%
# ----------------------------------------------------------------------------
# VARMAP
# ----------------------------------------------------------------------------
#%%                    
def GSLIB_VARMAP(varmapexe,parfl,datafl,outfl,gslib_params):       

    (nx,ny,nz,xsiz,ysiz,zsiz,icolx,icoly,icolz,nvar,ivar1,tmin,tmax,igrid,
     nxlag,nylag,nzlag,dxlag,dylag,dzlag,minpairs,standardize,nvarg,
     ivtail,ivhead,ivtype,cut) = gslib_params.values()       

   
    """GSLIB variogram program for irregularly spaced data
    
    Variograms are traditionally presented as 1D curves: g(h) as a function of the
    distance h along a particular direction. It is often useful to have a global view
    of the variogram values in all directions. The variogram map [28, 89] is a 2D
    plot of g(h1; h2) of the sample semivariogram for all experimentally available
    separation vectors h = (h1; h2)); see Figure III.5. The value g(0) = 0 plots at
    the center of the figures. The value g(h1; h2) plots as a grayscale or colorscale
    at offset (h1; h2) from that center. Any pixel not filled by an experimental
    value is left uninformed (it could be filled in with some type of interpolation).
    The variogram volume is the 3D plot of g(h1; h2; h3).
    Directions of anisotropy are usually evident from a variogram map; see
    Figure III.5. Visualization of a variogram volume is best done using a graph-
    ical routine that allows continuous sectioning of the volume along any of the
    three axes h1; h2; h3, [162].
        
    Parameters
    ----------
    gamvexe : string
        GSLIB executable for creation of a variogram considering irregularly 
        spaced data.
    parfl: string
        parameter file for VARMAP.
    varmappar : string
        parameter file for VARMAP.
    datafl : string
        input data in a simplied Geo-EAS formatted file.
    nvar : int
        nvar and ivar(1) . . . ivar(nvar): the number of variables and their
        column order in the data file.
    tmin : float
        lower trimmimg limit - all values less than this value are ignored.
    tmax : float
        upper trimmimg limit - all values greater than this value are ignored.
     igrid : integer
        set to 1 if the data are on a regular grid and set to 0 if the data
        are irregularly spaced.       
     nx, ny, nz: integer
        the number of nodes in the x, y, and z coordinate
        directions (used if igrid is set to 1). Any of the numbers may be set to
        one if that coordinate is not present in the data; e.g., 2D data may be
        handled by setting nz to 1.  
    xsiz, ysiz, zsiz: integer
        the separation distance between the nodes in the
        x, y, and z coordinate directions (used if igrid is set to 1).
    icolx : int
        the columns for the x coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data, e.g.,
        2D data may be handled by setting icolz to 0).
    icoly : int
        the columns for the y coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data - see
        iclox.
    icolz : int
        the columns for the y coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data - see
        iclox.
    outfl : string
        the output file for the variogram maps/volumes are written to a single 
        output file named out. This file contains each variogram
        map/volume written sequentially as a grid with the x direction cycling
        fastest, then y, then z.        
    nxlag, nylag, nzlag : integer
        the number of lags to compute in the X, Y and Z directions.
    dxlag, dylag, dzlag : float
        the lag tolerances or \cell sizes" in the X, Y and Z directions.
        minpairs: the minimum number of pairs needed to define a variogram
        value (set to missing if fewer than minpairs is found).
    standardize : integer
        if set to 1, the semivariogram values will be divided by
        the variance.
    nvarg : int
        the number of variograms to compute.
    ivtail : int
        specificaton of variables which are to be used for the tail.
    ivhead : int
        specificaton of variables which are to be used for the tail. For direct
        variograms the ivtail array is identical to the ivhead array. Cross 
        variograms are computed by having the tail variable different from the 
        head variable.
    ivtype : int
        The ivtype variable corresponds to the integer code in the list given 
        in Section III.1 of Deutsch and Journel.
        1 = Semivariogramm; -> predef. indicator variable
        2 = Cross semivariogram
        3 = Covariance
        4 = Correlogram
        5 = General relative semivariogramm
        6 = Pairwise relative semivariogram
        7 = Semivariogram of logarithms
        8 = Semimadogram
        
        9 = Indicator Semivariogram(continuous variable)
        10 = Indicator Semivariogram(categorical variable)
            cut: whenever the ivtype is set to 9 or 10, i.e., asking for an 
            indicator variogram, then a cutoff must be specified immediately 
            after the ivtype parameter on the same line in the input file 
            (see Figure III.3, Deutsch and Journel). Note that if an indicator 
            variogram is being computed, then the cutoff/category applies to 
            variable ivtail(i) in the input file [although the ivhead(i) variable 
            is not used, it must be present in the file to maintain consistency
            with the other variogram types].
    cut : float
        whenever the ivtype is set to 9 or 10, i.e., asking for an indicator
        variogram, then a cutoff must be specified immediately after the ivtype
        parameter on the same line in the input file (see Figure III.3). Note that
        if an indicator variogram is being computed, then the cutoff/category
        applies to variable ivtail(i) in the input file [although the ivhead(i)
        variable is not used, it must be present in the file to maintain consistency
        with the other variogram types].
            
    Returns
    -------
    None.

    """

    """
    #  Variogram of Irregularly Spaced 3-D Data, wrapper for GAMV.exe from GSLIB (.exe must be in working directory)
    Parameters for VARMAP
    *********************
    START OF PARAMETERS:
    ../data/cluster.dat \file with data
    1 3 \ number of variables: column numbers
    -1.0e21 1.0e21 \ trimming limits
    0 \1=regular grid, 0=scattered values
    50 50 1 \if =1: nx, ny, nz
    1.0 1.0 1.0 \ xsiz, ysiz, zsiz
    1 2 0 \if =0: columns for x,y, z coordinates
    varmap.out \file for variogram output
    10 10 0 \nxlag, nylag, nzlag
    5.0 5.0 1.0 \dxlag, dylag, dzlag
    5 \minimum number of pairs
    1 \standardize sill? (0=no, 1=yes)
    1 \number of variograms
    1 1 1 \tail, head, variogram type"""
    

    file = open(parfl, "w")
    file.write("                      Parameters for VARMAP                                 \n")
    file.write("                      *********************                                 \n")
    file.write("                                                                            \n")
    file.write("START OF PARAMETERS:                                                        \n")
    file.write(str(datafl) + "                     - file with data                           \n")
    file.write(str(nvar) + " " + str(ivar1) + " - number of varables,column numbers         \n")
    file.write(str(tmin) + " " + str(tmax) + "  - trimming limits                           \n")    
    file.write(str(igrid) + "                    - 1=regular grid, 0=scattered values        \n")

    file.write(str(nx) + " " + str(ny) +
        " " + str(nz) + " - if =1: nx, ny, nz                                               \n")
    file.write(str(xsiz) + " " + str(ysiz) +
        " " + str(zsiz) + " - xsiz, ysiz, zsiz                                              \n")
    file.write(str(icolx) + " " + str(icoly) +
        " " + str(icolz) + " - if =0: columns for x,y, z coordinates                        \n")    
    file.write(str(outfl) + "                    - file for variogram output           \n")
    file.write(str(nxlag) + " " + str(nylag) +
        " " + str(nzlag) + " - nxlag, nylag, nzlag                                          \n")        
    file.write(str(dxlag) + " " + str(dylag) +
        " " + str(dzlag) + " - dxlag, dylag, dzlag                                          \n") 
    file.write(str(minpairs) + "                    - minimum number of pairs               \n") 
    file.write(str(standardize) + "                   - standardize sills? (0=no, 1=yes)    \n")
    file.write(str(nvarg) + "                      - number of variograms                   \n")
    file.write(str(ivtail) + " " +
                   str(ivhead) + " " +
                   str(ivtype) + " - tail var., head var., variogram type                   \n")

    file.close()

def GSLIB_SGSIM(sgsimexe,datfl,parfl,outfl,GSLIB_SGSIM_params):
    """
    Sequential Gaussian Simulation program. It is used to model continuous 
    variables probabilistically by taking advantage of its remarkable features. 

    Parameters
    ----------
    sgsimexe : string
        program to perform a Sequential Gaussian Simulation.
    datfl : string
        the input continuous data in a simplified Geo-EAS formatted file.
    parfl : string
        SGSIM parameter file.
    icolx : int
        the columns for the x coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data, e.g.,
        2D data may be handled by setting icolz to 0).
    icoly : int
        the columns for the y coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data - see
        iclox.
    icolz : int
        the columns for the y coordinates (any of the column numbers 
        may be set to zero if that coordinate is not present in the data - see
        iclox.
    icolvr : int
        column number for variable to be simulated.
    icolwt : int
        column number for the declustering weight.
    icolsec : int
        column numbers for the secondary variable (e.g., for external drift 
        if used).
    tmin : float
        lower trimmimg limit - all values less than this value are ignored.
    tmax : float
        upper trimmimg limit - all values greater than this value are ignored.
    itrans : int
        if set to 0, then no transformation will be performed; the variable is 
        assumed already standard normal (the simulation results will also be 
        left unchanged). If itrans=1, transformations are performed.
    transfl : string
        output file for the transformation table if transformation is required 
        (itrans=0).
    ismooth : int
        if set to 0, then the data histogram, possibly with declustering 
        weights is used for transformation; if set to 1, then the data are 
        transformed according to the values in another file (perhaps from 
        histogram smoothing).
    smthfl : string
        file with the values to use for transformation to normal scores
        if ismooth is set to 1).
    icolvrsmthfl : int
        column in smth for the variable (set to 1 if smth is the output 
                                         from histsmth).
    icolwtsmthfl : int
        column in smth for the declustering weight (set to 2 if smth is the 
                                                    output from histsmth).
    zmin : float
        the minimum allowable data value. These are used in the 
        back-transformation procedure.
    zmax : float
        the maximum allowable data value. These are used in the 
        back-transformation procedure.
    ltail : int
        specify the back-transformation implementation in the lower tail of 
        the distribution: ltail=1 implements linear interpolation to the lower 
        limit zmin, and ltail=2 implements power model interpolation, with  
        ω = ltpar, to the lower limit zmin.
        The middle class interpolation is linear.
    ltpar : float
        see ltail description.
    utail : int
        utail and utpar specify the back-transformation implementation in the
        upper tail of the distribution: utail=1 implements linear interpolation
        to the upper limit zmax, utail=2 implements power model interpolation, 
        with ω = utpar, to the upper limit zmax, and utail=4 implements 
        hyperbolic model extrapolation with ω = utpar. The hyperbolic tail 
        extrapolation is limited by zmax.
    utpar : float
        see utail description.
    idbg : int
        debugging level between 0 and 3. The larger the debugging level, the 
        more information written out.
    dbgfl : string
        the file for the debugging output.
    outfl : string
        the output grid is written to this file. The output file will contain
        the results, cycling fastest on x, then y, then z, then simulation by
        simulation.
    nsim : int
        number of simulaitons to generate.
    nx : int
        definition of the grid system (x-axis)
        number of cells in x direction.
    ny : int
        definition of the grid system (y-axis)
        number of cells in y direction.
    nz : int
        definition of the grid system (z-axis)
        number of cells in z direction.
    xmn : float
        definition of the grid system (x-axis)
        coordinates at the center of the first blocks in x direction.
    ymn : float
        definition of the grid system (y-axis)
        coordinates at the center of the first blocks in y direction.
    zmn : float
        definition of the grid system (z-axis)
        coordinates at the center of the first blocks in z direction.
    xsiz : float
        definition of the grid system (x-axis) 
        constant block size in x direction.
    ysiz : float
        definition of the grid system (y-axis) 
        constant block size in y direction.
    zsiz : float
        definition of the grid system (y-axis) 
        constant block size in y direction.
    seed : int
        random number seed (a large odd integer e.g. 69069).
    ndmin : int
        the minimum number of original data that should be used to simulate a 
        grid node. If there are fewer than ndmin data points, the node is not 
        simulated.
    ndmax : int
        the maximum number of original data that should be used to simulate a 
        grid node. If there are fewer than ndmin data points, the node is not 
        simulated.
    ncnode : int
        the maximum number of previously simulated nodes to use for the 
        simulation of another node.
    sstrat : int
        if set to 0, the data and previously simulated grid nodes are searched 
        separately: The data are searched with a super block search, and the 
        previously simulated nodes are searched with a spiral search (see 
        Section II.4, Deutsch and Journel). If set to 1, the data are relocated 
        to grid nodes and  a spiral search is used and the parameters ndmin and 
        ndmax are not considered.
    multgrid : int
        multiple grid simulation will be performed if this is set to 1 
        (otherwise a standard spiral search for previously simulated nodes is 
        considered).
    nmult : int
        the number of multiple grid refinements to consider (used only if 
        multgrid is set to 1).
    noct : int
        the number of original data to use per octant. If this parameter is set 
        ≤ 0, then it is not used; otherwise, it overrides the ndmax parameter 
        and the data are partitioned into octants and the closest noct data in 
        each octant are retained for the simulation of a grid node.
    radiushmax : float
        the search radii in the maximum horizontal direction.
    radiushmin : float
        the search radii in the minimum horizontal direction.
    radiusvert : float
        the search radii in the vertical direction.
    sang1 : float
        the angle parameters that describe the orientation of the search 
        ellipsoid. sang1 rotates the original Y axis (principal direction) 
        in the horizontal plane.
    sang2 : float
        the angle parameters that describe the orientation of the search 
        ellipsoid. sang2 rotates the rotates the principal direction from 
        the horizontal.
    sang3 : float
        the angle parameters that describe the orientation of the search 
        ellipsoid. sang 3 leaves the principal direction, defined by ang1 and 
        ang2, unchanged. The two directions orthogonal to that principal 
        direction are rotated clockwise relative to the principal direction 
        when looking toward the origin.
    covtab1 : int
        Maximum x point in covariance table (odd number).
    covtab2 : int
        Maximum y point in covariance table (odd number).
    covtab3 : int
        Maximum z point in covariance table (odd number).
    ktype : int
        the kriging type (0 = simple kriging, 1 = ordinary kriging, 2 = simple 
        kriging with a locally varying mean, 3 = kriging with an external drift, 
        or 4 = collocated cokriging with one secondary variable) used 
        throughout the loop over all nodes. SK is required by theory; only in 
        cases where the number of original data found in the neighborhood is 
        large enough can OK be used without the risk of spreading data values 
        beyond their range of influence.
    rho : float
        correlation coefficient to use for collocated cokriging (used only 
        if ktype = 4).
    varred : float
        variance reduction factor to use for collocated cokriging (used only 
        if ktype = 4).
    secfl : TYPE
        the  le for the locally varying mean, the external drift variable, or 
        the secondary variable for collocated cokriging (the secondary variable 
        must be gridded at the same resolution as the model being constructed
        by sgsim).
    icolsecfl : int
        column for secondary variable.
    nst : int
        the number of semivariogram structures.
    c0 : float
        the isotropic nugget constant.
    variogram : list
        it: int
         the type of structure (the power model is not allowed)
             it = 1: spherical model,
             it = 2: exponential model
             it = 3: Gaussian model
             it = 5: hole effect model
        cc: float
         the c parameter - sill value
        ang1: float
         rotates the original Y axis (principal direction) in the horizontal 
         plane: This angle is measured in degrees clockwise.
        ang2: float
         rotates the principal direction from the horizontal: This angle is 
         measured in negative degrees down from horizontal.
        ang3: float
         leaves the principal direction, defined by ang1 and ang2, unchanged. 
         The two directions orthogonal to that principal direction are rotated 
         clockwise relative to the principal direction when looking toward 
         the origin.
        a_hmax: float
         the maximum horizontal range
        a_hmin: float
         the minimum horizontal range
        a_vert: float
         the vertical range

    Returns
    -------
    sim_array : TYPE
        DESCRIPTION.

    """
    (icolx,icoly,icolz,icolvr,icolwt,icolsec,tmin,tmax,itrans,transfl,ismooth,
     smthfl,icolvrsmthfl,icolwtsmthfl,zmin,zmax,ltail,ltpar,utail,utpar,idbg,
     dbgfl,nsim,grid_params,seed,ndmin,ndmax,ncnode,
     sstrat,multgrid,nmult,noct,radiushmax,radiushmin,radiusvert,sang1,sang2,
     sang3,covtab1,covtab2,covtab3,ktype,rho,varred,secfl,icolsecfl,nst,c0,
     variogram)  = GSLIB_SGSIM_params.values()
    nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn = grid_params.values()
    f = open(parfl, "w")
    f.write("              Parameters for SGSIM                                        \n")
    f.write("              ********************                                        \n")
    f.write("                                                                          \n")
    f.write("START OF PARAMETERS:                                                      \n")
    f.write(str(datfl) + "               -file with data                               \n")
    f.write(
        str(icolx) + " " + str(icoly) + " " + str(icolz) + " " + str(icolvr) +
        " " + str(icolwt) + " " + str(icolsec) +
        "                                 -  columns for X,Y,Z,vr,wt,sec.var.          \n")
    f.write(str(tmin)+" "+str(tmax)+"     -  trimming limits                           \n")
    f.write(str(itrans) + "               -transform the data (0=no, 1=yes)            \n")
    f.write(str(transfl) + "              -  file for output trans table               \n")
    f.write(str(ismooth) + "              -  consider ref. dist (0=no, 1=yes)          \n")
    f.write(str(smthfl) + "               -  file with ref. dist distribution          \n")
    f.write(str(icolvr) + " " +
            str(icolwt) + "               -  columns for vr and wt                     \n")
    f.write(str(zmin) + " " + str(zmax) + "   - zmin,zmax(tail extrapolation)          \n")
    f.write(str(ltail) + " " +
            str(ltpar) + "                -  lower tail option, parameter              \n")
    f.write(str(utail) + " " + str(utpar) + "   - upper tail option, parameter         \n")
    f.write(str(idbg) + "                 -debugging level: 0,1,2,3                    \n")
    f.write(str(dbgfl) + "                -file for debugging output                   \n")
    f.write(str(outfl) + "                -file for simulation output                  \n")
    f.write(str(nsim) + "                 -number of realizations to generate          \n")
    f.write(str(nx) + " " + str(xmn) + " " + str(xsiz) + "                             \n")
    f.write(str(ny) + " " + str(ymn) + " " + str(ysiz) + "                             \n")
    f.write(str(nz) + " " + str(zmn) + " " + str(zsiz) + "                             \n")
    f.write(str(seed) + "                 -random number seed                          \n")
    f.write(str(ndmin) + " " + str(ndmax) + "   -min and max original data for sim     \n")
    f.write(str(ncnode) + "               -number of simulated nodes to use            \n")
    f.write(str(sstrat) + "               -assign data to nodes (0=no, 1=yes)          \n")
    f.write(str(multgrid) + " " + 
            str(nmult) + "                -multiple grid search (0=no, 1=yes),num      \n")
    f.write(str(noct) + "                 -maximum data per octant (0=not used)        \n")
    f.write(str(radiushmax) + " " +
            str(radiushmin) + " " +
            str(radiusvert) + "           -maximum search  (hmax,hmin,vert) \n")
    f.write(str(sang1) + " " + 
            str(sang2) + " " +
            str(sang3) + "                -angles for search ellipsoid                 \n")
    f.write(str(covtab1)+" "+str(covtab2)+" "+
           str(covtab3)+"                 -size of covariance lookup table             \n")
    f.write(str(ktype) + " " + str(rho) + " " +
            str(varred) + "               - ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC       \n")
    f.write(str(secfl) + "                -  file with LVM, EXDR, or COLC variable     \n")
    f.write(str(icolsecfl) + "            -  column for secondary variable             \n")
    f.write(str(nst)+" "+str(c0)+ "    -  nst, nugget effect                        \n")
    for i in range(0,nst):
        it,cc,ang1,ang2,ang3,a_hmax,a_hmin,a_vert = variogram[i]

        f.write(str(it)+" "+str(cc)+" "+str(ang1)+
                   " "+str(ang2)+" "+
                   str(ang3)+"            -  "+str(i+1)+"  it,cc,ang1,ang2,ang3        \n")
        f.write("      "+str(a_hmax)+" "+
                   str(a_hmin)+" "+
                   str(a_vert) + "        -  "+str(i+1)+" - a_hmax, a_hmin, a_vert     \n")

    f.close()


def GSLIB_SISIM_LM(sisim_lmexe,parfl,datafl,priormfl,
                   dbgfl,output,GSLIB_SISIM_LM_params):
    
    """
    Sequential Indicator Simulation program. A technique that combines the 
    indicator theory with the sequential concept to simulate non-parametric 
    categorical distributions.
    
    Parameters
    ----------
    sisimexe : string
        program to perform a Sequential Gaussian Simulation.
    parfl : string
        parameter file?.
    vartype : int
        the variable type (1=continuous, 0=categorical).
    ncat : int
        DESCRIPTION.
    catx : float
        the threshold values or category codes (there should be ncat values on 
                                                this line of input).
    gpdfx : int
        the global cdf or pdf values (there should be ncat values on this line 
                                      of input).
    datafl : string
        the input data in a simplified Geo-EAS file.
    icolx : int
        column number for the x coordinates.
    icoly : int
        column number for the y coordinates.
    icolz : int
        column number for the z coordinates.
    icolvr : int
        column number for variable to be simulated.
    directik : string
        already transformed indicator values are read from this file. Missing 
        values are fied as less than tmin, which would correspond to a 
        constraint interval. Otherwise, the cdf data should steadily increase 
        from 0 to 1 and soft categorical probabilities must be between 0 to 1 
        and sum to 1.0.
    icolsx : int
        column for the x coordinates.
    icolsy : int
        column for the y coordinates.
    icolsz : int
        column for the z coordinates.
    icol : int
        columns for the indicator variables.
    imbsim : int
        set to 1 if considering Markov-Bayes option for cokriging with soft 
        indicator data, otherwise, set to 0.
    bz : float
        if imbsim is set to 1, then the B(z) calibration values are needed.
    tmin : float
        lower trimmimg limit - all values less than this value are ignored.
    tmax : float
        upper trimmimg limit - all values greater than this value are ignored.
    zmin : float
        minimum attribute values when considering a continuous variable.
    zmax : float
        minimum attribute values when considering a continuous variable.
    ltail : int
        ltail and ltpar specify the extrapolation in the lower tail: 
            ltail=1 implements linear interpolation to the lower limit zmin; 
            ltail=2 power model interpolation, with ω = ltpar, to the lower 
            limit zmin; 
            ltail = 3 implements linear interpolation between tabulated 
                quantiles (only for continuous variables).
    ltpar : float
        see ltail.
    middle : int
        middle and midpar specify the interpolation within the middle of the 
        distribution:
            middle = 1 implements linear interpolation; 
            middle = 2 implements power model interpolation, with ω = midpar; 
            middle = 3 allows for linear interpolation between tabulated 
            quantile values (only for continuous variables)            
    midpar : float
        see middle.
    utail : int
        utail and utpar specify the extrapolation in the upper tail of the
        distribution:
           utail=1 implements linear interpolation to the upper limit zmax, 
           utail=2 implements power model interpolation, with ω = utpar, to the 
           upper limit zmax,
           utail=3 implements linear interpolation between tabulated quantiles, 
           utail=4 implements hyperbolic model extrapolation with ω = utpar. 
           The hyperbolic tail extrapolation is limited by zmax (only for 
           continuous variables)
    utpar : float
        see utpar.
    tabfl : string
        If linear interpolation between tabulated values is the option selected 
        for any of the three regions, then this simplified Geo-EAS format file 
        is opened to read in the values. One legitimate choice is exactly the
        same file as the conditioning data, i.e., datafl. Note that tab specifies
        the tabulated values for all classes.
    icolvrt : int
        the column numbers for the values in tabfl.
    icolwtt : int
        the column numbers for the declustering weights in tabfl. Note that 
        declustering weights can be used but are not required - in this case 
        just set the column number less than or equal to zero.
    idbg : int
        an integer debugging level between 0 and 3. The larger the debugging 
        level, the more information written out.
    dbgfile : string
        the file for the debugging output.
    output : string
        the output grid is written to this file. The output file will contain
        the results, cycling fastest on x, then y, then z, then simulation by 
        simulation.
    nsim : int
        the number of simulations to generate.
    nx : int
        definition of the grid system (x-axis)
        number of cells in x direction.
    ny : int
        definition of the grid system (y-axis)
        number of cells in y direction.
    nz : int
        definition of the grid system (z-axis)
        number of cells in z direction.
    xmn : float
        definition of the grid system (x-axis)
        coordinates at the center of the first blocks in x direction.
    ymn : float
        definition of the grid system (y-axis)
        coordinates at the center of the first blocks in y direction.
    zmn : float
        definition of the grid system (z-axis)
        coordinates at the center of the first blocks in z direction.
    xsiz : float
        definition of the grid system (x-axis) 
        constant block size in x direction.
    ysiz : float
        definition of the grid system (y-axis) 
        constant block size in y direction.
    zsiz : float
        definition of the grid system (y-axis) 
        constant block size in y direction.
    seed : int
        random number seed (a large odd integer).
    ndmax : int
        the maximum number of original data that should be used to simulate a 
        grid node.
    ncnode : int
        the maximum number of previously simulated nodes to use for the 
        simulation of another node.
    maxsec : int
        the maximum number of soft data (at node locations) that will be used 
        for the simulation of a node. This is particularly useful to restrict 
        the number of soft data when an exhaustive secondary variable informs 
        all grid nodes.
    sstrat : int
        assign data to nodes. (0=no, 1=yes)    
        if set to 0, the data and previously simulated grid nodes are searched 
        separately: The data are searched with a super block search and the 
        previously simulated nodes are searched with a spiral search (see Section 
        II.4, Deutsch and Journel). If set to 1, the data are relocated to grid 
        nodes and a spiral search is used; the parameters ndmin and ndmax are 
        not considered.
    multgrid : int
        a multiple grid simulation will be performed if this is set to 1 
        (otherwise a standard spiral search will be considered).
    nmult : int
        the target number of multiple grid refinements to consider (used only 
        if multgrid is set to 1).
    noct : int
        the number of original data to use per octant. If this parameter is set 
        ≤ 0, then it is not used; otherwise, the closest noct data in each 
        octant are retained for the simulation of a grid node.
    radiushmax : float
       the search radii in the maximum horizontal direction.
    radiushmin : float
       the search radii in the minimum horizontal direction.
    radiusvert : float
       the search radii in the vertical direction.
    sang1 : float
       the angle parameters that describe the orientation of the search 
       ellipsoid. sang1 rotates the original Y axis (principal direction) 
       in the horizontal plane.
    sang2 : float
       the angle parameters that describe the orientation of the search 
       ellipsoid. sang2 rotates the rotates the principal direction from 
       the horizontal.
    sang3 : float
       the angle parameters that describe the orientation of the search 
       ellipsoid. sang 3 leaves the principal direction, defined by ang1 and 
       ang2, unchanged. The two directions orthogonal to that principal 
       direction are rotated clockwise relative to the principal direction 
       when looking toward the origin.
    covtab1 : int
        Maximum x point in covariance table (odd number).
    covtab2 : int
        Maximum y point in covariance table (odd number).
    covtab3 : int
        Maximum z point in covariance table (odd number).
    mik : int
        if mik is set to 0, then a full indicator kriging is performed at each 
        grid node location to establish the conditional distribution. If mik is 
        set to 1, then the median approximation is used; i.e., a single 
        variogram is used for all categories; therefore, only one kriging
        system needs to be solved and the computer time is significantly reduced. 
        The variogram corresponding to category mikcat will be used.
    mikcat : float
        threshold num IK.
    ktype : int
        the kriging type (0 = simple kriging, 1 = ordinary kriging) used 
        throughout the loop over all nodes. SK is required by theory, only in 
        cases where the number of original data found in the neighborhood is 
        large enough can OK be used without the risk of spreading data values 
        beyond their range of influence [87]. The global pdf values (specified
        with each category) are used for simple kriging.
    
    variogram : list
    
        The following set of parameters are required for each of the ncat categories:
        
        nst : int
            the number of semivariogram structures.
        c0 : float
            the isotropic nugget constant.
        it: int
            the type of structure (the power model is not allowed)
                 it = 1: spherical model,
                 it = 2: exponential model
                 it = 3: Gaussian model
                 it = 5: hole effect model
        cc: float
         the c parameter - sill value
        ang1: float
         rotates the original Y axis (principal direction) in the horizontal 
         plane: This angle is measured in degrees clockwise.
        ang2: float
         rotates the principal direction from the horizontal: This angle is 
         measured in negative degrees down from horizontal.
        ang3: float
         leaves the principal direction, defined by ang1 and ang2, unchanged. 
         The two directions orthogonal to that principal direction are rotated 
         clockwise relative to the principal direction when looking toward 
         the origin.
        a_hmax: float
         the maximum horizontal range
        a_hmin: float
         the minimum horizontal range
        a_vert: float
         the vertical range
         
         Each semivariogram model refers to the corresponding indicator transform.
         A Gaussian variogram with a small nugget constant is not a legitimate
         variogram model for a discontinuous indicator function. There is no need 
         to standardize the parameters to a sill of one since only the relative 
         shape affects the kriging weights.

    Returns
    -------
    None.

    """
    (vartype,ncat,catx,gpdfx,icolx,icoly,icolz,
     icolvr,tmin,tmax,zmin,zmax,ltail,ltpar,middle,midpar,
     utail,utpar,tabfl,icolvrt,icolwtt,idbg,nsim,
     grid_params,seed,ndmax,ncnode,sstrat,multgrid,
     nmult,noct,radiushmax,radiushmin,radiusvert,sang1,sang2,sang3,
     covtab1,covtab2,covtab3,mik,mikcat,variogram) = GSLIB_SISIM_LM_params.values()
    

    nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn = grid_params.values()
#  Wrapper for SISIM_LM.exe from GSLIB (.exe must be in working directory)

    file = open(parfl, "w")
    file.write("                      Parameters for SISIM_LM                                \n")
    file.write("                      *******************                                    \n")
    file.write("                                                                             \n")
    file.write("START OF PARAMETERS:                                                         \n")
    file.write(str(vartype)+"                      - 1=continuous(cdf), 0=categorical(pdf)   \n")
    file.write(str(ncat)+"                      - number thresholds/categories            \n")
    cat1, cat2, cat3, cat4, cat5 = catx
    file.write(str(cat1)+" "+str(cat2)+" "+str(cat3)+" "+
               str(cat4)+" "+str(cat5)+"                      - thresholds / categories               \n")
    gpdf1, gpdf2, gpdf3, gpdf4, gpdf5 = gpdfx
    file.write(str(gpdf1)+" "+str(gpdf2)+" "+str(gpdf3)+
               " "+str(gpdf4)+" "+str(gpdf5)+"                      - global cdf / pdf                      \n")
    file.write(str(datafl)+"                      - file with data                          \n")
    file.write(str(icolx)+" "+str(icoly)+" "+
               str(icolz)+" "+str(icolvr)+"                      - columns for X,Y,Z, and variable       \n")
    file.write(str(priormfl)+"                      - file with gridded indicator prior mean  \n")
    file.write(str(tmin)+" "+str(tmax)+"                      - trimming limits                         \n")
    file.write(str(zmin)+" "+str(zmax)+"                      - minimum and maximum data value          \n")
    file.write(str(ltail)+" "+str(ltpar)+"                      - lower tail option and parameter       \n")
    file.write(str(middle)+" "+str(midpar)+"                      - middle option and parameter       \n")
    file.write(str(utail)+" "+str(utpar)+"                      - upper tail option and parameter       \n")
    file.write(str(tabfl)+"                      -   file with tabulated values            \n")
    file.write(str(icolvrt)+" "+str(icolwtt)+"                      - columns for variable, weight       \n")
    file.write(str(idbg)+"                      - debugging level: 0,1,2,3                \n")
    file.write(str(dbgfl)+"                      - file for debugging output               \n")
    file.write(str(output)+"                      - file for simulation output              \n")
    file.write(str(nsim)+"                      - number of realizations                  \n")
    file.write(str(nx)+" "+str(xmn)+" "+str(xsiz)+"                      - nx,xmn,xsiz                           \n")
    file.write(str(ny)+" "+str(ymn)+" "+str(ysiz)+"                      - ny,ymn,ysiz                           \n")
    file.write(str(nz)+" "+str(zmn)+" "+str(zsiz)+"                      - nz,zmn,zsiz                           \n")
    file.write(str(seed)+"                      - random number seed                      \n")
    file.write(str(ndmax)+"                      - maximum original data  for each kriging \n")
    file.write(str(ncnode)+"                      - maximum previous nodes for each kriging \n")
    file.write(str(sstrat)+"                      - assign data to nodes? (0=no,1=yes)      \n")
    file.write(str(multgrid)+" "+str(nmult)+"                      - multiple grid search? (0=no,1=yes),num  \n")
    file.write(str(noct)+"                      - maximum per octant    (0=not used)      \n")
    file.write(str(radiushmax)+" "+str(radiushmin)+
               " "+str(radiusvert)+"                      - maximum search radii                    \n")
    file.write(str(sang1)+" "+str(sang2)+" "+
               str(sang3)+"                      - angles for search ellipsoid             \n")
    file.write(str(covtab1)+" "+str(covtab2)+" "+
               str(covtab3)+"                      - -size of covariance lookup table        \n")
    file.write(str(mik)+" "+str(mikcat)+"                      - 0=full IK, 1=median approx. (cutoff)    \n")
    for i in range(0,ncat):
        nst_var1, c0, nest_var = variogram[i]
        
        file.write(str(nst_var1)+" "+str(c0)+ "          \ "+str(i+1)+" nst, nugget effect         \n")
        for nest in nest_var:
            it,cc,ang1,ang2,ang3,a_hmax,a_hmin,a_vert = nest
            file.write(str(it)+" "+str(cc)+" "+str(ang1)+
                       " "+str(ang2)+" "+str(ang3)+"   \      it,cc,ang1,ang2,ang3       \n")
            file.write("      "+str(a_hmax)+" "+
                       str(a_hmin)+" "+str(a_vert)+"   \      a_hmax, a_hmin, a_vert             \n")
    file.close()

#sisimexe = str(WD) + "gslib90/sisim_lm.exe"  # sisim executable
#parfl = str(WD) + "parfiles/sisim.par"  # parameter file
#datafl = str(WD) + "save/tertiary.gslib"  # file with data
#dbgfl = str(WD) + "save/sisim.dbg"  # file for debugging output
#output = str(WD) + "save/sisim.out"  # file for kriging output
#GSLIB_SISIM_params = {
#    "vartype": 0, "ncat": 2, "catx": [0, 1, "", "", ""], "gpdfx": [0.5, 0.5, "", "", ""],
#    "icolx": 1, "icoly": 2, "icolz": 3, "icolvr": 4, "directik": direct.ik, "icolsx": 1,
#    "icolsy": 2, "icolsz": 3, "icol": 4, "imbsim": 5, "bz": 0,
#    "tmin": 0, "tmax": 1.0e21, "zmin": 0, "zmax": 30, "ltail": 1, "ltbar": 0,
#    "middle": 1, "midpar": 1, "utail": 1, "utpar": 1, "tabfl": "cluster.dat", "icolvrt": 3,
#    "icolwtt": 0, "idbg": 1, "nsim": 1, "grid_params": grid_params, "seed": 12345,
#    "ndmax": 40, "ncnode": 40, "maxsec": 0, "sstrat": 1, "multgrid": 0, "nmult": 0,
#    "noct": 0, "radiushmax": 100, "radiushmin": 100, "radiusvert": 10, "sang1": 90,
#    "sang2": 0, "sang3": 0, "covtab1": 50, "covtab2": 50, "covtab3": 50, "mik": 0,
#    "mikcat": 1, "ktype": 1, "variogram": [[1, 0.3, 2, 1, 90, 0.0, 0.0, 70, 70, 7],
#                                              [1, 0.3, 2, 1, 90, 0.0, 0.0, 70, 70, 7]]
#}

def GSLIB_SISIM(sisimexe,parfl,datafl,dbgfile,output,GSLIB_SISIM_params):
    """
    Sequential Indicator Simulation program. A technique that combines the 
    indicator theory with the sequential concept to simulate non-parametric 
    categorical distributions.
    
    Parameters
    ----------
    sisimexe : string
        program to perform a Sequential Gaussian Simulation.
    parfl : string
        parameter file?.
    vartype : int
        the variable type (1=continuous, 0=categorical).
    ncat : int
        DESCRIPTION.
    catx : float
        the threshold values or category codes (there should be ncat values on 
                                                this line of input).
    gpdfx : int
        the global cdf or pdf values (there should be ncat values on this line 
                                      of input).
    datafl : string
        the input data in a simplified Geo-EAS file.
    icolx : int
        column number for the x coordinates.
    icoly : int
        column number for the y coordinates.
    icolz : int
        column number for the z coordinates.
    icolvr : int
        column number for variable to be simulated.
    directik : string
        already transformed indicator values are read from this file. Missing 
        values are fied as less than tmin, which would correspond to a 
        constraint interval. Otherwise, the cdf data should steadily increase 
        from 0 to 1 and soft categorical probabilities must be between 0 to 1 
        and sum to 1.0.
    icolsx : int
        column for the x coordinates.
    icolsy : int
        column for the y coordinates.
    icolsz : int
        column for the z coordinates.
    icol : int
        columns for the indicator variables.
    imbsim : int
        set to 1 if considering Markov-Bayes option for cokriging with soft 
        indicator data, otherwise, set to 0.
    bz : float
        if imbsim is set to 1, then the B(z) calibration values are needed.
    tmin : float
        lower trimmimg limit - all values less than this value are ignored.
    tmax : float
        upper trimmimg limit - all values greater than this value are ignored.
    zmin : float
        minimum attribute values when considering a continuous variable.
    zmax : float
        minimum attribute values when considering a continuous variable.
    ltail : int
        ltail and ltpar specify the extrapolation in the lower tail: 
            ltail=1 implements linear interpolation to the lower limit zmin; 
            ltail=2 power model interpolation, with ω = ltpar, to the lower 
            limit zmin; 
            ltail = 3 implements linear interpolation between tabulated 
                quantiles (only for continuous variables).
    ltpar : float
        see ltail.
    middle : int
        middle and midpar specify the interpolation within the middle of the 
        distribution:
            middle = 1 implements linear interpolation; 
            middle = 2 implements power model interpolation, with ω = midpar; 
            middle = 3 allows for linear interpolation between tabulated 
            quantile values (only for continuous variables)            
    midpar : float
        see middle.
    utail : int
        utail and utpar specify the extrapolation in the upper tail of the
        distribution:
           utail=1 implements linear interpolation to the upper limit zmax, 
           utail=2 implements power model interpolation, with ω = utpar, to the 
           upper limit zmax,
           utail=3 implements linear interpolation between tabulated quantiles, 
           utail=4 implements hyperbolic model extrapolation with ω = utpar. 
           The hyperbolic tail extrapolation is limited by zmax (only for 
           continuous variables)
    utpar : float
        see utpar.
    tabfl : string
        If linear interpolation between tabulated values is the option selected 
        for any of the three regions, then this simplified Geo-EAS format file 
        is opened to read in the values. One legitimate choice is exactly the
        same file as the conditioning data, i.e., datafl. Note that tab specifies
        the tabulated values for all classes.
    icolvrt : int
        the column numbers for the values in tabfl.
    icolwtt : int
        the column numbers for the declustering weights in tabfl. Note that 
        declustering weights can be used but are not required - in this case 
        just set the column number less than or equal to zero.
    idbg : int
        an integer debugging level between 0 and 3. The larger the debugging 
        level, the more information written out.
    dbgfile : string
        the file for the debugging output.
    output : string
        the output grid is written to this file. The output file will contain
        the results, cycling fastest on x, then y, then z, then simulation by 
        simulation.
    nsim : int
        the number of simulations to generate.
    nx : int
        definition of the grid system (x-axis)
        number of cells in x direction.
    ny : int
        definition of the grid system (y-axis)
        number of cells in y direction.
    nz : int
        definition of the grid system (z-axis)
        number of cells in z direction.
    xmn : float
        definition of the grid system (x-axis)
        coordinates at the center of the first blocks in x direction.
    ymn : float
        definition of the grid system (y-axis)
        coordinates at the center of the first blocks in y direction.
    zmn : float
        definition of the grid system (z-axis)
        coordinates at the center of the first blocks in z direction.
    xsiz : float
        definition of the grid system (x-axis) 
        constant block size in x direction.
    ysiz : float
        definition of the grid system (y-axis) 
        constant block size in y direction.
    zsiz : float
        definition of the grid system (y-axis) 
        constant block size in y direction.
    seed : int
        random number seed (a large odd integer).
    ndmax : int
        the maximum number of original data that should be used to simulate a 
        grid node.
    ncnode : int
        the maximum number of previously simulated nodes to use for the 
        simulation of another node.
    maxsec : int
        the maximum number of soft data (at node locations) that will be used 
        for the simulation of a node. This is particularly useful to restrict 
        the number of soft data when an exhaustive secondary variable informs 
        all grid nodes.
    sstrat : int
        assign data to nodes. (0=no, 1=yes)    
        if set to 0, the data and previously simulated grid nodes are searched 
        separately: The data are searched with a super block search and the 
        previously simulated nodes are searched with a spiral search (see Section 
        II.4, Deutsch and Journel). If set to 1, the data are relocated to grid 
        nodes and a spiral search is used; the parameters ndmin and ndmax are 
        not considered.
    multgrid : int
        a multiple grid simulation will be performed if this is set to 1 
        (otherwise a standard spiral search will be considered).
    nmult : int
        the target number of multiple grid refinements to consider (used only 
        if multgrid is set to 1).
    noct : int
        the number of original data to use per octant. If this parameter is set 
        ≤ 0, then it is not used; otherwise, the closest noct data in each 
        octant are retained for the simulation of a grid node.
    radiushmax : float
       the search radii in the maximum horizontal direction.
    radiushmin : float
       the search radii in the minimum horizontal direction.
    radiusvert : float
       the search radii in the vertical direction.
    sang1 : float
       the angle parameters that describe the orientation of the search 
       ellipsoid. sang1 rotates the original Y axis (principal direction) 
       in the horizontal plane.
    sang2 : float
       the angle parameters that describe the orientation of the search 
       ellipsoid. sang2 rotates the rotates the principal direction from 
       the horizontal.
    sang3 : float
       the angle parameters that describe the orientation of the search 
       ellipsoid. sang 3 leaves the principal direction, defined by ang1 and 
       ang2, unchanged. The two directions orthogonal to that principal 
       direction are rotated clockwise relative to the principal direction 
       when looking toward the origin.
    covtab1 : int
        Maximum x point in covariance table (odd number).
    covtab2 : int
        Maximum y point in covariance table (odd number).
    covtab3 : int
        Maximum z point in covariance table (odd number).
    mik : int
        if mik is set to 0, then a full indicator kriging is performed at each 
        grid node location to establish the conditional distribution. If mik is 
        set to 1, then the median approximation is used; i.e., a single 
        variogram is used for all categories; therefore, only one kriging
        system needs to be solved and the computer time is significantly reduced. 
        The variogram corresponding to category mikcat will be used.
    mikcat : float
        threshold num IK.
    ktype : int
        the kriging type (0 = simple kriging, 1 = ordinary kriging) used 
        throughout the loop over all nodes. SK is required by theory, only in 
        cases where the number of original data found in the neighborhood is 
        large enough can OK be used without the risk of spreading data values 
        beyond their range of influence [87]. The global pdf values (specified
        with each category) are used for simple kriging.
    
    variogram : list
    
        The following set of parameters are required for each of the ncat categories:
        
        nst : int
            the number of semivariogram structures.
        c0 : float
            the isotropic nugget constant.
        it: int
            the type of structure (the power model is not allowed)
                 it = 1: spherical model,
                 it = 2: exponential model
                 it = 3: Gaussian model
                 it = 5: hole effect model
        cc: float
         the c parameter - sill value
        ang1: float
         rotates the original Y axis (principal direction) in the horizontal 
         plane: This angle is measured in degrees clockwise.
        ang2: float
         rotates the principal direction from the horizontal: This angle is 
         measured in negative degrees down from horizontal.
        ang3: float
         leaves the principal direction, defined by ang1 and ang2, unchanged. 
         The two directions orthogonal to that principal direction are rotated 
         clockwise relative to the principal direction when looking toward 
         the origin.
        a_hmax: float
         the maximum horizontal range
        a_hmin: float
         the minimum horizontal range
        a_vert: float
         the vertical range
         
         Each semivariogram model refers to the corresponding indicator transform.
         A Gaussian variogram with a small nugget constant is not a legitimate
         variogram model for a discontinuous indicator function. There is no need 
         to standardize the parameters to a sill of one since only the relative 
         shape affects the kriging weights.

    Returns
    -------
    None.

    """
               
    (vartype,ncat,catx,gpdfx,icolx,icoly,icolz,icolvr,directik,icolsx,
    icolsy,icolsz,icol,imbsim,bz,tmin,tmax,zmin,zmax,ltail,ltbar,middle,midpar,utail,utpar,
    tabfl,icolvrt,icolwtt,idbg,nsim,grid_params,
    seed,ndmax,ncnode,maxsec,sstrat,multgrid,nmult,noct,radiushmax,radiushmin,radiusvert,
    sang1,sang2,sang3,covtab1,covtab2,covtab3,mik,mikcat,ktype,variogram) = GSLIB_SISIM_params.values()
    
    nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn = grid_params.values()

 
    file = open(parfl, "w")
    file.write("                      Parameters for SISIM                                   \n")
    file.write("                      *******************                                    \n")
    file.write("                                                                             \n")
    file.write("START OF PARAMETERS:                                                         \n")
    file.write(str(vartype)+"                      \ 1=continuous(cdf), 0=categorical(pdf)   \n")
    file.write(str(ncat)+"                         \ number thresholds/categories            \n")
    cat1, cat2, cat3, cat4, cat5 = catx
    file.write(str(cat1)+" "+str(cat2)+" "+str(cat3)+" "+
               str(cat4)+" "+str(cat5)+"           \   thresholds / categories               \n")
    gpdf1, gpdf2, gpdf3, gpdf4, gpdf5 = gpdfx
    file.write(str(gpdf1)+" "+str(gpdf2)+" "+str(gpdf3)+
               " "+str(gpdf4)+" "+str(gpdf5)+"     \   global cdf / pdf                      \n")
    file.write(str(datafl)+"                       \ file with data                          \n")
    file.write(str(icolx)+" "+str(icoly)+" "+
               str(icolz)+" "+str(icolvr)+"        \   columns for X,Y,Z, and variable       \n")
    file.write(str(directik)+"                     \ file with soft indicator input          \n")
    icol1, icol2, icol3, icol4 = icol
    file.write(str(icolsx)+" "+str(icolsy)+" "+
               
               str(icolsz)+" "+
               str(icol1)+" "+str(icol2)+" "+
               str(icol3)+" "+str(icol4)+" \   columns for X,Y,Z, and indicators             \n")
    file.write(str(imbsim)+"                       \   Markov-Bayes simulation (0=no,1=yes)  \n")
    bz1, bz2, bz3, bz4, bz5 = bz
    file.write(str(bz1)+" "+str(bz2)+" "+str(bz3)+
               " "+str(bz4)+" "+str(bz5)+"         \      calibration B(z) values            \n")
    file.write(str(tmin)+" "+str(tmax)+"           \ trimming limits                         \n")
    file.write(str(zmin)+" "+str(zmax)+"           \ minimum and maximum data value          \n")
    file.write(str(ltail)+" "+str(ltbar)+"         \   lower tail option and parameter       \n")
    file.write(str(middle)+" "+str(midpar)+"       \   middle     option and parameter       \n")
    file.write(str(utail)+" "+str(utpar)+"         \   upper tail option and parameter       \n")
    file.write(str(tabfl)+"                        \   file with tabulated values            \n")
    file.write(str(icolvrt)+" "+str(icolwtt)+"     \      columns for variable, weight       \n")
    file.write(str(idbg)+"                         \ debugging level: 0,1,2,3                \n")
    file.write(str(dbgfile)+"                      \ file for debugging output               \n")
    file.write(str(output)+"                       \ file for simulation output              \n")
    file.write(str(nsim)+"                         \ number of realizations                  \n")
    file.write(str(nx)+" "+str(xmn)+" "+str(xsiz)+"  \ nx,xmn,xsiz                           \n")
    file.write(str(ny)+" "+str(ymn)+" "+str(ysiz)+"  \ ny,ymn,ysiz                           \n")
    file.write(str(nz)+" "+str(zmn)+" "+str(zsiz)+"  \ nz,zmn,zsiz                           \n")
    file.write(str(seed)+"                         \ random number seed                      \n")
    file.write(str(ndmax)+"                        \ maximum original data  for each kriging \n")
    file.write(str(ncnode)+"                       \ maximum previous nodes for each kriging \n")
    file.write(str(maxsec)+"                       \ maximum soft indicator nodes for kriging\n")
    file.write(str(sstrat)+"                       \ assign data to nodes? (0=no,1=yes)      \n")
    file.write(str(multgrid)+" "+str(nmult)+"      \ multiple grid search? (0=no,1=yes),num  \n")
    file.write(str(noct)+"                         \ maximum per octant    (0=not used)      \n")
    file.write(str(radiushmax)+" "+str(radiushmin)+
               " "+str(radiusvert)+"               \ maximum search radii                    \n")
    file.write(str(sang1)+" "+str(sang2)+" "+
               str(sang3)+"                        \ angles for search ellipsoid             \n")
    file.write(str(covtab1)+" "+str(covtab2)+" "+
               str(covtab3)+"                        \ -size of covariance lookup table             \n")
    file.write(str(mik)+" "+str(mikcat)+"          \ 0=full IK, 1=median approx. (cutoff)    \n")
    file.write(str(ktype)+"                        \ 0=SK, 1=OK                              \n")
    for i in range(0,ncat):
        nst_var1, c0, nest_var = variogram[i]
        
        file.write(str(nst_var1)+" "+str(c0)+ "          \ "+str(i+1)+" nst, nugget effect         \n")
        for nest in nest_var:
            it,cc,ang1,ang2,ang3,a_hmax,a_hmin,a_vert = nest
            file.write(str(it)+" "+str(cc)+" "+str(ang1)+
                       " "+str(ang2)+" "+str(ang3)+"   \      it,cc,ang1,ang2,ang3       \n")
            file.write("      "+str(a_hmax)+" "+
                       str(a_hmin)+" "+str(a_vert)+"   \      a_hmax, a_hmin, a_vert             \n")

    file.close()