# Evidence for the software submission checklist

This checklist describes the current public repository, including explicit
FAI/FCI outputs, the updated core workflows and HBV sensitivity experiments.
The verification records below were generated on 20 October 2026.

| Checklist item | Repository evidence |
| --- | --- |
| Source code | `code/01_flood_extraction/` through `code/05_hbv_sensitivity/` |
| Small demonstration dataset | `examples/case_data/`; HBV forcings are generated internally |
| Dependencies and operating-system versions | Root README, `requirements-core.txt`, HBV requirements and measured run reports |
| Tested software versions | Root README, `docs/core_verification.json` and the HBV complete-run report |
| Non-standard hardware | None for the supplied CPU demonstrations |
| Installation instructions | Root README and `docs/HBV_README.md` |
| Typical installation time | Root README gives a labelled 5-15-minute estimate; installation was not timed |
| Demo instructions | Root README, `examples/example_four_methods.py` and HBV `run_all.py` |
| Expected outputs | `expected_outputs/core/` and `code/05_hbv_sensitivity/expected_outputs/` |
| Expected demo runtime | Measured core and HBV runtimes in the README and verification records |
| Use on other data | Root README commands, input contract, unit conventions and HBV parameter instructions |
| Optional reproduction of every quantitative result | Not met: complete global datasets and every model/figure pipeline are outside this repository |

The four core demos completed in 5.69 seconds and the five core-definition
tests present on 20 October passed in 1.20 seconds. The current test suite also
includes output/window tests. The complete HBV workflow, including the figure, took 1.948
seconds; all 14 labelled scenario results matched the final experiment's
reference values. These are single-run measurements on Windows 11. Installation
time, other operating systems and independent-colleague installation testing
were not measured or tested.

- [Installation, dependencies and usage](../README.md)
- [Core verification record](core_verification.json)
- [Algorithm definitions and version history](ALGORITHMS_AND_VERSION_NOTES.md)
- [HBV methods and validation](HBV_README.md)
- [HBV complete-run report](../code/05_hbv_sensitivity/expected_outputs/run_report.json)
- [MIT License for the code](../LICENSE)

Repository: https://github.com/Bright3098/global-flood-scaling-organization

Algorithm descriptions appear in the manuscript Methods and Supplementary
Notes 1 and 3-6. HBV-specific locations are Supplementary Note 6 and Extended
Data Fig. 5. Data-provider terms apply to the original observational datasets.
An archival DOI has not been assigned.
