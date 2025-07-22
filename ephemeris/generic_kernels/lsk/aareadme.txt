

    LEAPSECONDS KERNEL VERSION NAIF0011

The file naif0011.tls is a unix-style text file. It is suitable for use
on all unix boxes, including PC/linux and MAC/OSX machines.

For PCs running Windows, use naif0011.tls.pc.

Use of one of these files is required for all SPICE computations
dealing with times on or after July 01, 2015 00:00:00 UTC if you want
accurate conversions between UTC and Ephemeris Time (a.k.a. TDB).
Failure to use naif0011 for times after this epoch will result
in a one second time conversion error.

You may begin using naif0011 in place of naif0010 right now; there
is no need to wait until July 01, 2015. Time conversions for epochs
prior to July 01, 2015 will be the same using either naif0010 or
naif0011.


