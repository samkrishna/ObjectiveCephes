# ObjectiveCephes

This is an API wrapper around the [Cephes library](http://www.netlib.org/cephes/index.html).

The inspiration for doing this came from Evan Miller's blog post, [Statistical Shortcomings in Standard Math Libraries (And How To Fix Them)](http://www.evanmiller.org/statistical-shortcomings-in-standard-math-libraries.html).

Trying to figure out how to build Cephes in such a way that it's importable or available to [Homebrew](https://brew.sh) seemed to be too heavy a lift for me, esp. with the Data Science revolution going on in the Mobile era. The Cephes library's last **version update** was back in 1992.

Since Evan made such a compelling argument, I thought I'd take a crack at this. In the mean time, there's a lot of warnings to clean up and API wrappers to write, so that's what I'm doing here.
