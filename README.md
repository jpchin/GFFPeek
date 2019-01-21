# GFFPeek

A python tool for taking quick peeks into GFF files.  Takes a GFF file and one or two strings as input.  Searches the "attributes" property of each gene in the GFF file for the query string (not case sensitive).  It then draws a diagram with proportional gene lengths and spacings of the gene of interest plus 5 genes up and down stream.

#Current state:

Interesting, but not reads for serious use.

# Dependencies

1) Requires the Python Pillow library

# Known issues

1) Very short genes lead to malformed arrows.
2) Needs a better way to display gene annotations (currently those at the right hadn side get cut off by the edge of the image).
3) Something funny is going on with the "title" text in the upper left corner.
