# LunaSRIM documentation

lunasrim is a pure Python 3 re-implementation of SRIM-2013 that reproduces `TRIM.exe` and `SRModule.exe` digit for digit. The [README](../README.md) gives an overview and a quick start.

- [User guide](user-guide.md): installation, the Python API (materials, stopping tables, energy loss, TRIM calculations, reading results, pysrim-compatible classes, SRIM input files), the command line, examples, performance and limits.
- [Physics](physics.md): the models of SRIM-2013 as lunasrim implements them, from electronic and nuclear stopping to the TRIM Monte Carlo, energy straggling and implantation profiles, with their approximations.
- [Validation](validation.md): how lunasrim is compared with the original programs, the 46 TRIM and 24 SR reference cases, the processor dependence of the x87 arithmetic and the known limits.
- [Web application](web.md): the tools of the web application, the online version, running it locally or with Docker, and the deployment of the static site.
- [Development](development.md): repository layout, tests, new reference cases, documentation and the release process.

The web application runs online at https://skowrons94.github.io/LunaSRIM/.
