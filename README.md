# Automating-Nickel-1m
Automating Observing Procedures for Lick Observatory's Nickel 1-m telescope 

Scott Hakoda
Utah Tech University 

Site: University of California Observatories, Santa Cruz, California 
Mentors: Kyle Westfall, Will Deich

Automating Observing Procedures for the Nickel 1-m Telescope at Lick Observatory

Lick Observatory’s Nickel 1-m Telescope, located atop Mount Hamilton, is utilized by a wide variety of observers, including faculty, staff, and students both within and outside of the University of California Observatories. At the beginning of each observation session, observers must perform a set of startup and calibration procedures. These include focusing, adjusting the telescope’s pointing, and collecting bias frames and sky flats. Currently, the focusing process is inefficient, requiring observers to select a focus star and manually image it at a range of different focus values. This is a process that takes approximately 10 minutes and must be repeated multiple times per night, making it a tedious and error-prone task. This project aims to address this inefficiency by automating the focusing process using Python scripting. The Keck Task Library (KTL) API will be used for communication with the Nickel 1-m Telescope and manipulation of its instruments.  The automation system will select and slew to the nearest known focus star based on the Nickel’s current position, and take a series of exposures of the star at a range of focus values. After each exposure, the quality of the focus is assessed by measuring the full width at half maximum (FWHM) of the star at each focus value. Based on the collected measurements, a curve is fit to the data to determine the optimal focus value. This should significantly reduce telescope downtime spent on focusing by up to 90% and improve scientific return on data by standardizing the image quality. 

