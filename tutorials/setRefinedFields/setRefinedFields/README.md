# setRefinedFields tutorial

## Notes

This case is mean to illustrate many of the features of setRefinedFields including refinement selection (hex/polyhedral), set and zone creation, setting boundary values, and non-uniform region values. The debug flag is used so that intermediate steps are written along with additional fields. By default a serial case using hexehedral refinement will be run, but the '-t/-tet' flag can be used to use a  tet mesh (created with GMSH). The '-j/-parallel' flag can be used to run the case in parallel.

Because 'system/setRefinedFields' is modified, the 'system/setFieldsDict_mod' is used to actually run the case