# ObjectiveCephes

This is an API wrapper around the [Cephes library](http://www.netlib.org/cephes/index.html).

The inspiration for doing this came from Evan Miller's blog post, [Statistical Shortcomings in Standard Math Libraries (And How To Fix Them)](http://www.evanmiller.org/statistical-shortcomings-in-standard-math-libraries.html).

I've used the sources from the Perl library, [Math-Cephes](https://github.com/shlomif/Math-Cephes), to get the Cephes library to build on macOS 10.12. Not sure how to go about modifying the code to have it buildable on ARM, so if you have any ideas, a PR is welcome.

Next step is to port the Math-Cephes test cases to XCTest-based Unit Tests. Here's the list of test files from Math-Cephes:

- ~~bessels.t~~
- ~~betas.t~~
- ~~cmplx.t~~
- ~~dists.t~~
- ~~elliptics.t~~
- ~~explog.t~~
- fract.t
- gammas.t
- ~~hypergeometrics.t~~
- hypers.t
- mat.t
- misc.t
- new_cmplx-2.t
- new_cmplx.t
- ~~poly.t~~ (Can't figure out how to test this)
- ~~style-trailing-space.t~~ (Can't test in Perl. Not going to attempt in Objective-C)
- trig.t
- utils.t

I've also taken the liberty to put the Cephes Math C functions into their own "namespace". Currently the vast majority of the functions are prefixed by the `cfs_` prefix. The intention here is that when the Objective-C / Swift API wrappers are built, all possible collisions with the `math.h` library will be avoided. **NOTE:** Most, but not all functions have been namespaced.

PRs are welcome.
