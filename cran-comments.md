## Test environments
* local OS X install(x86_64, darwin17.0), R 4.0.5
* win-builder (release, old release, development)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs. 

## Downstream dependencies
There are currently no downstream dependencies for this package

## Resubmission
This is a resubmission. In this version I have addressed comments by Gregor Seyer:

* Please do not start the description with "This package", package name,
title or similar. \
**Done**


* Please always write package names, software names and API (application
programming interface) names in single quotes in title and description.
e.g: --> 'CINmetrics' \
**Done**

* If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <[https:...]https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title") \
**Done**

* Please add small executable examples in your Rd-files to illustrate the
use of the exported function but also enable automatic testing. \
**Done**
