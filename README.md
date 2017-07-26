# ObjectiveCephes

This is an API wrapper around the [Cephes library](http://www.netlib.org/cephes/index.html).

The inspiration for doing this came from Evan Miller's blog post, [Statistical Shortcomings in Standard Math Libraries (And How To Fix Them)](http://www.evanmiller.org/statistical-shortcomings-in-standard-math-libraries.html).

I've used the sources from the Perl library, [Math-Cephes](https://github.com/shlomif/Math-Cephes) to get the Cephes library to build on macOS 10.12. Not sure how to go about modifying the code to have it buildable on ARM, so if you have any ideas, a PR is welcome.

