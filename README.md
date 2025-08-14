# Automating-Nickel-1m

Scott Hakoda
Utah Tech University 

Site: University of California Observatories, Santa Cruz, California 
Mentors: Kyle Westfall, Will Deich

Automating Observing Procedures for the Nickel 1-m Telescope at Lick Observatory

Lick Observatory’s Nickel 1-m Telescope, located atop Mount Hamilton, is utilized by a wide variety of observers, including faculty, staff, and students both within and outside of the University of California Observatories. At the beginning of each observation session, observers must perform a set of startup and calibration procedures. These include focusing, adjusting the telescope’s pointing, and collecting bias frames and sky flats. Currently, the focusing process is inefficient, requiring observers to select a focus star and manually image it at a range of different focus values. This is a process that takes approximately 10 minutes and must be repeated multiple times per night, making it a tedious and error-prone task. This project aims to address this inefficiency by automating the focusing process using Python scripting. The Keck Task Library (KTL) API was used for communication with the Nickel 1-m Telescope and manipulation of its instruments.  The new automation system selects and slews to the nearest known focus star based on the Nickel’s current position and takes a series of exposures of the star at a range of focus values. After each exposure, the quality of the focus is assessed by measuring the full width at half maximum (FWHM) of the star at each focus value. Based on the collected measurements, a curve is fit to the data to determine the optimal focus value. The new focusing procedure should significantly reduce telescope downtime spent on focusing by up to 90% and improve scientific return on data by standardizing the image quality. 


move_to_target.py:
    automatically moves the Nickel telescope to the nearest suitable star for focus or pointing calibration based on a catalog of known stars.

    Uses the Nickel's current RA and DEC to locate the closest pointing/focus (p/f) star from the star list (point_focus.txt)
    Creates a SkyCoord object for the telescope and each of the p/f stars and uses the SkyCoord method .seperation to get the distance between the telescope's position and each of the p/f stars. The p/f star with the smallest seperation is returned and the telescope is moved to its RA and DEC.

    Argument	    Type	Default	    Description
    --star_type     str	    'focus'	    Type of star: 'focus' or 'pointing'
    --dry_run	    flag	False	    Test mode - no actual telescope movement

    Examples

    1. Move to focus star:
        kpython move_to_target.py --star_type focus
    
    2. Test pointing star slelction:
        kpython move_to_target.py --star_type pointing --dry_run


photometry.py:
    analyzes FITS images to measure star properties, particularly the Full Width at Half Maximum (FWHM) for focus optimization.

    1. Iterative Background Subtraction
        - uses SigmaClip, detect_threshold, and detect_sources to determine sources and grow a grow a mask around them to be used when determining the background. 
        - median background is subtracted from the data 
        - process repeated until max iterations or threshold change is insignificant 
        - background subtracted data and detcted sources are returned
    
    2. Source Shape Analysis. 
        - use moments to characterize source propertis 
        - calculate M0, M1, and M2
        - σ_x = √(M2_x - M1_x²)
        - σ_y = √(M2_y - M1_y²)
        - FWHM_x = 2√(2ln(2)) * σ_x ≈ 2.355 * σ_x
        - FWHM_y = 2√(2ln(2)) * σ_y ≈ 2.355 * σ_y

    3. Source Selection
        a. Coordinate-based (if focus_coords provided):
            - calculate distance to each detected source 
            - select closest source to specified coordinates
        b. brightness-based:
            - select brightest source (highest M0)
            - used when coordinates not specified
    
    returns:    
        focus_star = {
            'Label': source_id,
            'M0': total_flux,
            'Centroid': (x_position, y_position),
            'FWHM': average_fwhm (from moment analysis),
            'Focus': focus_value,
            'ObsNum': observation_number
        }


automate.py:
    focus finding system for the Nickel telescope that automatically determines the optimal focus position by taking a series of exposures at different focus values and analyzing the FWHM of stars. 

    Keyword Class:
        manages all KTL keywords

    Event Class:
        executes telescope operation and image acquisition 
        a sequence consists of:
            1. change focus to target value
            2. take exposure
            3. perform photometery on resulting image 

    Focus Finding:
        - automatic
            1. intial sampling: takes two images at starting focus and starting focus + step_size
            2. analyze each images fwhm
            3. curve following: uses curve_finder to traverse the curve in the direction that minimizes the fwhm
            4. range expansion: uses curve_helper to fill in additional points aroudn the minimum 
            5. fit quadratic curve to collected data
        - manual
            1. takes images at regular intervals from start to end focus
            2. analyze each images fwhm 
            3. fit quadratic curve to collected data

    Program Flow:
        1. Parse command line arguments
        2. Create Keyword object (KTL interface)
        3. Perform safety checks (controller ready, motion enabled)
        4. Execute focus finding:
            - Auto mode: curve_finder algorithm
            - Manual mode: evalute focus over range
            - Refit/Reevaluate: data reprocessing
        5. Analyze results:
            - Detect outliers
            - Fit quadratic curve
            - Find optimal focus
        6. Save data and display results

    Argument	            Type    Default	    Description
    -fs, --focus_start	    int	    350	        Starting focus value
    -fe, --focus_end	    int	    None	    Ending focus value 
    -s, --step_size	        int	    5	        Step size for focus increments
    -el, --length_exposure	float	1.0	        Exposure length in seconds
    -es, --exposure_speed	str	    'Fast'	    Exposure speed: Slow, Medium, Fast
    --focus_coords	        float   None	    Coordinates(x, y) of star to use for focus
    --refit	                flag	False	    Refit curve using last focus session data
    --omit	                int 	None	    Observation numbers to omit from fitting
    --reevaluate	        flag	False	    Reprocess last focus session data
    --verbose	            flag	False	    Enable detailed debug output

    Examples

    1. automatic focus finding: starting at 360 incrementing focus by 2
        kpython automate.py -fs 360 -s 2

    2. specified range focus finding: starting at 360 and stopping at 370 incrementing focus by 2 (6 exposures)
        kpython automate.py -fs 360 -fe 370 -s 2

    3. refit curve: excluding bad measurements:
        kpython automate.py --refit --omit 1001 1002 1003

    4. reevalute curve: using specific coordinate of focus star
        kpython automate.py --reevaluate --focus_coords 480 580 