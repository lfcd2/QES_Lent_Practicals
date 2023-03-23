# Changes

Models are placed in separate files and imported for readability.

Models now take a list of dicts rather than a whole range of args.  **TODO: Use \*args to make this optional**

lfcd2OceanTools serves to add more helper functions (and a .py version of working.pyc):

1. add_emissions - adds emissions to the atmosphere dict when given a list of dicts.
2. add_fancy_labels - formats the ylabels for common variables
3. Modifier Class - allows variables in a list of dicts to be easily modified
4. modify_single_dict and modify_dicts - apply modifiers to a list of dicts

Some instances of applying fluxes to boxes have been rewritten to reduce the number of if statements used for speed.

Additionally, extra comments have been added throughout.