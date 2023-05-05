## R CMD check results

0 errors | 0 warnings | 4 notes

* The first note indicates that this is a new release.
* The second note only occurs on rhub with the Windows platform. I believe it is due to [R-hub issue #503](https://github.com/r-hub/rhub/issues/503).
* The third note only occurs on rhub with Fedora Linux. I believe this is due to a missing dependency ("tidy") on rhub's test machine. 
* The fourth note occurs on rhub with (both) Linux platforms. I believe CRAN's servers are faster and will not produce this note.

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

The notes from `rhub::check_for_cran()` are copied below

```
❯ On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Dean Koch <dkoch@ualberta.ca>'
  
  New submission
```

```
❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
```

```
❯ On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
```

```
❯ On fedora-clang-devel (r-devel)
  checking examples ... [25s/54s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
             user system elapsed
  sk_export 7.628  0.164  16.732
```
