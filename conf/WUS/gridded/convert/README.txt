DO NOT USE THESE FILES FOR FORTRAN CALCULATIONS
THEY EXIST ONLY TO FACILITATE CONVERSION

The files in this directory have been consolidated from different focal-mech
specific flavors of WUS grid inputs.

Comments '1! have been added where necessary.

These files have also had the default mMax adjusted to 7.5 for any GR MFDs.
The default mMax was set to 7.0 and was always being overidden by the value
in any supplied mMax file. This ultimately yielded the correct MFDs but
mMax was specified for almost every node.