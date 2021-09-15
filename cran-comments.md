## Test environments
* local R installation, R 4.1.1
* Linux (Ubuntu 20.04), macOS (10.15) and Windows (Server 2019), R devel and release (through GitHub Actions)
* Windows Server 2008, Ubuntu 20.04, Fedora Linux (devtools::check_rhub())
* Windows (devtools::check_win_devel())

## R CMD check results

0 errors | 0 warnings | 1 note

- New submission

- Possibly mis-spelled words in DESCRIPTION:
  Ferte (21:85)
  GBJ (18:69)
  Hejblum (22:20)
  sGBJ (19:20, 19:42)
  Thiebaut (22:5)

We believe the spelling to be correct

## Response to CRAN comments

Thanks for the comments, we adressed all of them.

- Please reduce the length of the title to less than 65 characters.

Ok, we did as requested

- Please rather use the Authors@R field and declare Maintainer, Authors and Contributors with their appropriate roles with person() calls.

Ok, we did as requested

- Please always explain all acronyms in the description text.

Ok, we did as requested

- Please write references in the description of the DESCRIPTION file in the form authors (year) <doi:...> (If you want to add a title as well please put it in quotes: "Title")

Ok, we did as requested
