## R CMD check results

0 errors | 0 warnings | 2 notes

* The first note indicates that this is a new release
* the second note only occurs on rhub with the windows platform. I believe it is due to [R-hub issue #503](https://github.com/r-hub/rhub/issues/503).

## Test environments

with `devtools::check()`

* Windows 10 Pro (build 19045), R-release, 64 bit

with `devtools::check_mac_release()`

* Mac OS X 11.5.2, R-release, 64 bit

with `rhub::check_for_cran()`

* Windows Server 2022, R-devel, 64 bit
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran

## Notes

The note from `rhub::check_for_cran()` is copied here

```
❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
```
